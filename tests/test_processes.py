"""Tests for SimbioProcess — the real simbio CRN bridge."""

import numpy as np
import pytest
from process_bigraph import Composite, allocate_core, gather_emitter_results

from pbg_simbio import SimbioProcess, build_crn_model

REVERSIBLE = {
    "species": {"A": 1.0, "B": 2.0, "AB": 0.0},
    "reactions": [
        {"name": "bind", "reactants": ["A", "B"], "products": ["AB"], "rate": 1.0},
        {"name": "unbind", "reactants": ["AB"], "products": ["A", "B"], "rate": 0.3},
    ],
    "volume": 1.0,
}


def test_build_crn_model():
    model, reaction_names = build_crn_model(
        REVERSIBLE["species"], REVERSIBLE["reactions"], REVERSIBLE["volume"]
    )
    assert reaction_names == ["bind", "unbind"]
    # species are real attributes on the dynamically built Compartment
    for name in ("A", "B", "AB"):
        assert hasattr(model, name)


def test_ports_are_dicts():
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    assert isinstance(proc.inputs(), dict)
    assert isinstance(proc.outputs(), dict)
    assert set(proc.inputs()) == {"concentrations", "rates"}
    assert set(proc.outputs()) == {"concentrations"}


def test_initial_state():
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    assert proc.initial_state()["concentrations"] == {"A": 1.0, "B": 2.0, "AB": 0.0}


def test_update_returns_deltas():
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    out = proc.update({"concentrations": {"A": 1.0, "B": 2.0, "AB": 0.0}}, interval=1.0)
    d = out["concentrations"]
    assert set(d) == {"A", "B", "AB"}
    # A + B -> AB consumes A and B, produces AB
    assert d["AB"] > 0
    assert d["A"] < 0 and d["B"] < 0


def test_mass_conservation():
    """A and B are conserved in (free + bound) form by the closed CRN."""
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    out = proc.update({"concentrations": {"A": 1.0, "B": 2.0, "AB": 0.0}}, interval=2.0)
    d = out["concentrations"]
    # d[AB] of AB formed must equal A consumed and B consumed
    assert abs(d["AB"] + d["A"]) < 1e-6
    assert abs(d["AB"] + d["B"]) < 1e-6


def test_rate_override_input():
    """A higher 'bind' rate fed via the rates port produces more product."""
    proc = SimbioProcess(config=REVERSIBLE, core=allocate_core())
    state = {"concentrations": {"A": 1.0, "B": 2.0, "AB": 0.0}}
    slow = proc.update({**state, "rates": {"bind": 0.5}}, interval=1.0)
    fast = proc.update({**state, "rates": {"bind": 5.0}}, interval=1.0)
    assert fast["concentrations"]["AB"] > slow["concentrations"]["AB"]


def test_composite_run_accumulates_in_store():
    """Deltas accumulate in the shared store across a real Composite run."""
    core = allocate_core()
    doc = {
        "simbio": {
            "_type": "process",
            "address": "local:SimbioProcess",
            "config": REVERSIBLE,
            "interval": 0.5,
            "inputs": {"concentrations": ["stores", "concentrations"]},
            "outputs": {"concentrations": ["stores", "concentrations"]},
        },
        "stores": {"concentrations": dict(REVERSIBLE["species"])},
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {"emit": {"concentrations": "map[string,float]", "time": "float"}},
            "inputs": {
                "concentrations": ["stores", "concentrations"],
                "time": ["global_time"],
            },
        },
    }
    sim = Composite({"state": doc}, core=core)
    sim.run(5.0)
    rows = gather_emitter_results(sim)[("emitter",)]
    assert len(rows) > 1
    final = rows[-1]["concentrations"]
    # AB grows toward steady state; A+AB and B+AB stay conserved
    assert final["AB"] > 0.5
    assert abs((final["A"] + final["AB"]) - 1.0) < 1e-6
    assert abs((final["B"] + final["AB"]) - 2.0) < 1e-6


def test_stoichiometry():
    """2A -> A2 consumes A twice as fast as it makes A2."""
    cfg = {
        "species": {"A": 2.0, "A2": 0.0},
        "reactions": [
            {"name": "dimerize", "reactants": [["A", 2]], "products": ["A2"], "rate": 1.0}
        ],
        "volume": 1.0,
    }
    proc = SimbioProcess(config=cfg, core=allocate_core())
    out = proc.update({"concentrations": {"A": 2.0, "A2": 0.0}}, interval=1.0)["concentrations"]
    assert abs(out["A"] + 2 * out["A2"]) < 1e-6  # 2 A consumed per A2 formed
