"""Tests for the simbio composite generators and their registration."""

from process_bigraph import Composite, allocate_core, gather_emitter_results

from pbg_simbio.composites.crn import michaelis_menten, reversible_binding


def test_reversible_binding_generator_registered():
    from pbg_superpowers.composite_generator import _REGISTRY

    matches = [eid for eid in _REGISTRY if eid.endswith(".simbio_reversible_binding")]
    assert matches, f"generator missing; have {list(_REGISTRY)[:5]}"


def test_michaelis_menten_generator_registered():
    from pbg_superpowers.composite_generator import _REGISTRY

    matches = [eid for eid in _REGISTRY if eid.endswith(".simbio_michaelis_menten")]
    assert matches, f"generator missing; have {list(_REGISTRY)[:5]}"


def test_reversible_binding_runs():
    core = allocate_core()
    doc = reversible_binding(kf=1.0, kr=0.2, interval=0.5)
    sim = Composite({"state": doc}, core=core)
    sim.run(3.0)
    rows = gather_emitter_results(sim)[("emitter",)]
    assert rows[-1]["concentrations"]["AB"] > 0


def test_michaelis_menten_produces_product():
    core = allocate_core()
    doc = michaelis_menten(kon=1.0, koff=0.5, kcat=2.0, e0=0.2, s0=5.0, interval=0.5)
    sim = Composite({"state": doc}, core=core)
    sim.run(10.0)
    rows = gather_emitter_results(sim)[("emitter",)]
    final = rows[-1]["concentrations"]
    assert final["P"] > 0  # substrate is turned into product
    # enzyme is conserved (free + bound)
    assert abs((final["E"] + final["ES"]) - 0.2) < 1e-6
