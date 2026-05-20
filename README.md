# pbg-simbio

A [process-bigraph](https://github.com/vivarium-collective/process-bigraph)
wrapper for [**simbio**](https://github.com/dyscolab/simbio), a Python package
for simulating **Chemical Reaction Networks (CRNs)** built on top of
[`poincare`](https://github.com/dyscolab/poincare).

`SimbioProcess` is a **real bridge**: every `update()` assembles a genuine
`simbio.Compartment` from your reaction spec and integrates it with simbio's
LSODA solver. No kinetics are reimplemented — simbio builds and solves the
ODEs.

## What it does

You describe a CRN declaratively (species + mass-action reactions). The process
exposes the species concentrations as a `map[string,float]` store that the
surrounding bigraph owns, integrates the network forward by each `interval`,
and emits the **change** in concentration per species as an additive delta — so
a sibling process (influx, dilution, transport, another reaction module, a
controller modulating rate constants) composes naturally with it.

## Installation

```bash
# From PyPI (once published):
pip install pbg-simbio
# or with uv:
uv pip install pbg-simbio

# For development (editable):
uv venv --python 3.12 .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
```

simbio requires Python ≥ 3.12.

> Once installed, the process registers automatically via
> `bigraph_schema.package.discover` — no manual `register_link()` calls needed.

## Quick start

```python
from process_bigraph import Composite, allocate_core, gather_emitter_results
from pbg_simbio.composites.crn import reversible_binding

core = allocate_core()
doc = reversible_binding(kf=1.0, kr=0.2, a0=1.0, b0=2.0, interval=0.5)
sim = Composite({"state": doc}, core=core)
sim.run(10.0)

results = gather_emitter_results(sim)
print(results[("emitter",)][-1])   # final concentrations + time
```

Or drive the process directly with your own CRN:

```python
from pbg_simbio import SimbioProcess

proc = SimbioProcess(config={
    "species": {"A": 1.0, "B": 2.0, "AB": 0.0},
    "reactions": [
        {"name": "bind",   "reactants": ["A", "B"], "products": ["AB"], "rate": 1.0},
        {"name": "unbind", "reactants": ["AB"], "products": ["A", "B"], "rate": 0.3},
    ],
    "volume": 1.0,
}, core=core)

delta = proc.update({"concentrations": {"A": 1.0, "B": 2.0, "AB": 0.0}}, interval=1.0)
```

## API reference

### `SimbioProcess` (process-bigraph `Process`)

| Config | Type | Default | Meaning |
|---|---|---|---|
| `species` | `map[string,float]` | `{}` | Species name → initial concentration |
| `reactions` | `list` | `[]` | Mass-action reaction specs (see below) |
| `volume` | `float` | `1.0` | Compartment volume |

A reaction spec is a dict:

```python
{
  "name": "bind",                  # optional, defaults to r0, r1, ...
  "reactants": ["A", "B"],          # name | ["A", 2] | {"species": "A", "stoichiometry": 2}
  "products":  ["AB"],
  "rate": 1.0,                      # mass-action rate constant
}
```

| Port | Direction | Type | Meaning |
|---|---|---|---|
| `concentrations` | input | `map[string,float]` | Current absolute concentrations (integration initial condition) |
| `rates` | input | `map[string,float]` | Optional per-reaction rate-constant overrides, keyed by reaction name |
| `concentrations` | output | `map[string,float]` | Per-species concentration **delta** over the interval (composes additively) |

### Composite generators

Discoverable via `pbg_superpowers.composite_generator.discover_generators()`:

- `simbio_reversible_binding` — `A + B <-> AB` reversible binding.
- `simbio_michaelis_menten` — mass-action enzyme kinetics `E + S <-> ES -> E + P`.

## Architecture

```
shared store "concentrations" (map[string,float], absolute)
        │  read (initial condition)            ▲  delta (accumulates)
        ▼                                       │
   SimbioProcess.update(interval) ── builds simbio.Compartment, solves LSODA
```

simbio concept → PBG mapping:

| simbio | pbg-simbio |
|---|---|
| `Compartment` / `Species` / `MassAction` | built dynamically from `config["reactions"]` |
| `Simulator(Model).solve(values=…, t_span=…)` | called once per `update()` |
| absolute concentration trajectory | emitted as per-step **deltas** so the bigraph composes |

## Demo

```bash
source .venv/bin/activate
python demo/demo_report.py     # writes and opens demo/report.html
```

The report runs three configurations (reversible binding, Michaelis–Menten,
and a rate-stressed variant), with time-series charts, metrics, an architecture
diagram, and an interactive PBG document tree.

## Limitations and assumptions

- Reactions are **mass-action** (`MassAction`); custom rate laws are not yet
  surfaced through the config schema (simbio itself supports `RateLaw` /
  `AbsoluteRateLaw`).
- The process re-solves over `(0, interval)` each step from the store's current
  concentrations; very large intervals reduce the resolution of any in-step
  dynamics that a sibling never observes (the emitted delta is still exact for
  the closed reaction system over that interval).
- simbio requires Python ≥ 3.12.
