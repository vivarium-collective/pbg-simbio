"""Composite generators for simbio Chemical Reaction Networks.

Each generator returns a process-bigraph document wiring a :class:`SimbioProcess`
to a shared ``concentrations`` store plus a RAM emitter, so the dashboard's
Composites tab can run and sweep them.
"""

from __future__ import annotations

from pbg_superpowers.composite_generator import composite_generator


def _crn_document(*, species, reactions, volume, interval):
    """Wire a SimbioProcess + emitter around a shared concentrations store."""
    return {
        "simbio": {
            "_type": "process",
            "address": "local:SimbioProcess",
            "config": {
                "species": species,
                "reactions": reactions,
                "volume": volume,
            },
            "interval": interval,
            "inputs": {"concentrations": ["stores", "concentrations"]},
            "outputs": {"concentrations": ["stores", "concentrations"]},
        },
        "stores": {"concentrations": dict(species)},
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


@composite_generator(
    name="simbio_reversible_binding",
    description="A + B <-> AB reversible mass-action binding, integrated by simbio.",
    parameters={
        "kf": {"type": "float", "default": 1.0,
               "description": "Forward (association) rate constant"},
        "kr": {"type": "float", "default": 0.2,
               "description": "Reverse (dissociation) rate constant"},
        "a0": {"type": "float", "default": 1.0, "description": "Initial [A]"},
        "b0": {"type": "float", "default": 2.0, "description": "Initial [B]"},
        "interval": {"type": "float", "default": 1.0,
                     "description": "Composite tick size"},
    },
)
def reversible_binding(core=None, *, kf=1.0, kr=0.2, a0=1.0, b0=2.0, interval=1.0):
    species = {"A": a0, "B": b0, "AB": 0.0}
    reactions = [
        {"name": "bind", "reactants": ["A", "B"], "products": ["AB"], "rate": kf},
        {"name": "unbind", "reactants": ["AB"], "products": ["A", "B"], "rate": kr},
    ]
    return _crn_document(species=species, reactions=reactions, volume=1.0,
                         interval=interval)


@composite_generator(
    name="simbio_michaelis_menten",
    description="Mass-action enzyme kinetics E+S<->ES->E+P, integrated by simbio.",
    parameters={
        "kon": {"type": "float", "default": 1.0,
                "description": "E+S association rate"},
        "koff": {"type": "float", "default": 0.5,
                 "description": "ES dissociation rate"},
        "kcat": {"type": "float", "default": 2.0,
                 "description": "Catalytic turnover rate ES->E+P"},
        "e0": {"type": "float", "default": 0.2, "description": "Initial enzyme [E]"},
        "s0": {"type": "float", "default": 5.0, "description": "Initial substrate [S]"},
        "interval": {"type": "float", "default": 1.0,
                     "description": "Composite tick size"},
    },
)
def michaelis_menten(core=None, *, kon=1.0, koff=0.5, kcat=2.0, e0=0.2, s0=5.0,
                     interval=1.0):
    species = {"E": e0, "S": s0, "ES": 0.0, "P": 0.0}
    reactions = [
        {"name": "bind", "reactants": ["E", "S"], "products": ["ES"], "rate": kon},
        {"name": "unbind", "reactants": ["ES"], "products": ["E", "S"], "rate": koff},
        {"name": "cat", "reactants": ["ES"], "products": ["E", "P"], "rate": kcat},
    ]
    return _crn_document(species=species, reactions=reactions, volume=1.0,
                         interval=interval)
