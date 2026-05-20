"""Process-bigraph wrapper for `simbio` Chemical Reaction Network simulation.

This is a **real bridge**: ``SimbioProcess.update`` builds a genuine
``simbio.Compartment`` from the reaction spec in its config and integrates it
with simbio's LSODA solver (via `poincare`). No reaction kinetics are
reimplemented here — the ODEs are assembled and solved by simbio itself.

simbio: https://github.com/dyscolab/simbio  (extends `poincare`)
"""

from __future__ import annotations

from typing import Any

import numpy as np
from process_bigraph import Process


def _normalize_reactant(entry: Any) -> tuple[str, int]:
    """Accept ``"A"``, ``["A", 2]`` or ``{"species": "A", "stoichiometry": 2}``.

    Returns ``(species_name, stoichiometry)``.
    """
    if isinstance(entry, str):
        return entry, 1
    if isinstance(entry, dict):
        return entry["species"], int(entry.get("stoichiometry", 1))
    if isinstance(entry, (list, tuple)):
        name = entry[0]
        stoich = int(entry[1]) if len(entry) > 1 else 1
        return name, stoich
    raise TypeError(f"Cannot interpret reactant {entry!r}")


def build_crn_model(
    species: dict[str, float],
    reactions: list[dict[str, Any]],
    volume: float = 1.0,
):
    """Assemble a `simbio.Compartment` subclass from a declarative CRN spec.

    Parameters
    ----------
    species
        Mapping of species name -> initial concentration.
    reactions
        Each reaction is a dict with ``reactants``, ``products`` (lists of
        species names, ``[name, stoich]`` pairs, or ``{"species", "stoichiometry"}``
        dicts), a numeric ``rate`` (mass-action rate constant), and an optional
        ``name`` (defaults to ``r0``, ``r1``, ...).
    volume
        Compartment volume.

    Returns
    -------
    (Model, reaction_names)
        The dynamically-created `Compartment` subclass and the ordered list of
        reaction names (whose rate parameters can be overridden at solve time).
    """
    from simbio import Compartment, MassAction, Parameter, Species, Volume, assign
    from simbio.core import concentration, volume as volume_field

    namespace: dict[str, Any] = {"__annotations__": {}}
    namespace["volume"] = volume_field(default=float(volume))
    namespace["__annotations__"]["volume"] = Volume

    species_objs: dict[str, Any] = {}
    for name, initial in species.items():
        obj = concentration(default=float(initial))
        namespace[name] = obj
        namespace["__annotations__"][name] = Species
        species_objs[name] = obj

    reaction_names: list[str] = []
    for index, rxn in enumerate(reactions):
        name = rxn.get("name", f"r{index}")
        reaction_names.append(name)

        rate_param = f"_rate_{name}"
        namespace[rate_param] = assign(default=float(rxn.get("rate", 1.0)))
        namespace["__annotations__"][rate_param] = Parameter

        def _terms(entries):
            terms = []
            for entry in entries:
                sp_name, stoich = _normalize_reactant(entry)
                if sp_name not in species_objs:
                    raise KeyError(
                        f"Reaction references unknown species {sp_name!r}; "
                        f"declare it in `species`."
                    )
                term = species_objs[sp_name]
                terms.append(stoich * term if stoich != 1 else term)
            return terms

        namespace[f"eq_{name}"] = MassAction(
            reactants=_terms(rxn.get("reactants", [])),
            products=_terms(rxn.get("products", [])),
            rate=namespace[rate_param],
        )

    model = type("SimbioCRNModel", (Compartment,), namespace)
    return model, reaction_names


class SimbioProcess(Process):
    """Time-stepped bridge to a simbio Chemical Reaction Network.

    Two equivalent ways to define the network:

    * **Antimony** — set ``config["antimony"]`` to an Antimony model string. It
      is compiled to SBML (libantimony), extracted (libSBML), and rebuilt with
      simbio's core (see :mod:`pbg_simbio.antimony_loader`). Arbitrary rate
      laws (mass-action, Hill, ...) are supported.
    * **Reaction spec** — set ``config["species"]`` + ``config["reactions"]``
      for a mass-action network built directly via :func:`build_crn_model`.

    The surrounding bigraph owns the absolute species concentrations in a
    shared store. Each step the process reads those concentrations, integrates
    the network forward by ``interval`` with simbio, and emits the **change** in
    each concentration as a ``map[string,float]`` delta — so a sibling process
    (an influx, dilution, transport term, or another reaction module) writing
    the same store composes additively.

    Inputs
    ------
    concentrations : map[string,float]
        Current absolute concentration of each species (used as the integration
        initial condition). Missing species fall back to their model initial.
    parameters : map[string,float]
        Optional overrides for the model's parameters (rate constants, Hill
        coefficients, ...), keyed by parameter name. A controller, a
        temperature model, or a slow regulatory process can write these.

    Outputs
    -------
    concentrations : map[string,float]
        Per-species concentration delta over the interval (accumulates
        additively in the shared store).
    """

    config_schema = {
        # Antimony model string (primary path). If empty, species+reactions are used.
        "antimony": {"_type": "string", "_default": ""},
        "species": {"_type": "map[float]", "_default": {}},
        # List of reaction dicts; kept permissive so arbitrary CRN topologies
        # pass through bigraph-schema unchanged.
        "reactions": {"_type": "list", "_default": []},
        "volume": {"_type": "float", "_default": 1.0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config=config, core=core)
        self._model = None
        self._simulator = None
        self._species_names: list[str] = []
        self._param_names: list[str] = []
        self._species_initials: dict[str, float] = {}
        # For the reaction-spec path, reaction names carry _rate_<name> params.
        self._reaction_names: list[str] = []
        # Pre-compute the initial state so initial_state() works before _build().
        if not self.config["antimony"]:
            self._species_initials = dict(self.config["species"])

    def _build(self):
        from simbio import Parameter, Simulator

        if self.config["antimony"]:
            from .antimony_loader import model_from_antimony

            self._model, self._species_names, self._param_names = model_from_antimony(
                self.config["antimony"]
            )
            self._species_initials = {
                name: float(getattr(self._model, name).initial or 0.0)
                for name in self._species_names
            }
        else:
            self._model, self._reaction_names = build_crn_model(
                species=self.config["species"],
                reactions=self.config["reactions"],
                volume=self.config["volume"],
            )
            self._species_names = list(self.config["species"])
            self._species_initials = dict(self.config["species"])
            # reaction-spec rate constants are exposed as _rate_<name> params
            self._param_names = list(self._reaction_names)

        self._simulator = Simulator(self._model)

    def _param_object(self, name):
        """Resolve a parameter name to its simbio Parameter on the model."""
        for attr in (name, f"_rate_{name}"):
            obj = getattr(self._model, attr, None)
            if obj is not None:
                return obj
        return None

    def inputs(self):
        return {
            "concentrations": "map[string,float]",
            "parameters": "map[string,float]",
        }

    def outputs(self):
        return {"concentrations": "map[string,float]"}

    def initial_state(self):
        if self._simulator is None:
            self._build()
        return {"concentrations": dict(self._species_initials)}

    def update(self, state, interval):
        if self._simulator is None:
            self._build()

        current = {**self._species_initials, **(state.get("concentrations") or {})}
        param_overrides = state.get("parameters") or {}

        values: dict[Any, float] = {}
        for name in self._species_names:
            values[getattr(self._model, name)] = float(current[name])
        for name, value in param_overrides.items():
            param = self._param_object(name)
            if param is not None:
                values[param] = float(value)

        result = self._simulator.solve(
            values=values,
            t_span=(0.0, float(interval)),
            save_at=[0.0, float(interval)],
        )

        delta = {}
        for name in self._species_names:
            final = float(np.asarray(result[name])[-1])
            delta[name] = final - float(current[name])
        return {"concentrations": delta}
