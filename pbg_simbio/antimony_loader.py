"""Load a genuine simbio model from an Antimony string.

Why not ``simbio.io.sbml.loads``? simbio 1.1.0's bundled SBML/MathML importer
is broken against the ``poincare``/``symbolite`` versions it resolves on PyPI
(``symbolite.core.as_function`` was removed; ``add_species`` builds a bare
``Variable`` instead of a ``Species``; ``get_species_reference`` calls the
keyword-only ``Species`` positionally as if it were a ``Reactant``; and the
``MathMLSymbol`` rate-law type fails poincare's unit translation). Repairing
all of that means maintaining a fork of simbio's importer.

Instead this module keeps the bridge **fully real** without touching that
importer: the model text is parsed by **libantimony** (real), the reaction
network is extracted with **libSBML** (real), and the model is rebuilt with
simbio's own working core constructors (`Species`, `Parameter`, `RateLaw`) —
so simbio still assembles and integrates the ODEs. Each reaction's kinetic law
is evaluated into a genuine simbio/symbolite expression, so arbitrary rate laws
(mass-action, Hill, Michaelis–Menten, ...) are supported, not just mass-action.
"""

from __future__ import annotations

import math
from typing import Any


def antimony_to_sbml(antimony: str, *, model_name: str | None = None) -> tuple[str, str]:
    """Compile an Antimony string to SBML via libantimony. Returns (sbml, name)."""
    import antimony as sb

    sb.clearPreviousLoads()
    code = sb.loadAntimonyString(antimony)
    if code < 0:
        raise ValueError(f"Antimony parse error: {sb.getLastError()}")
    name = model_name or sb.getMainModuleName()
    return sb.getSBMLString(name), name


def _eval_rate_law(formula: str, names: dict[str, Any]):
    """Evaluate a libSBML L3 infix formula into a symbolite expression.

    ``names`` maps every species / parameter / compartment id to its simbio
    object (or a float for a fixed compartment size). simbio's `Species` and
    `Parameter` overload the arithmetic operators, so ordinary Python
    evaluation of the formula builds the symbolic rate law.
    """
    python_expr = formula.replace("^", "**")
    return eval(python_expr, {"__builtins__": {}}, names)


def model_from_sbml(sbml: str, name: str = "model"):
    """Build a simbio `Compartment` subclass from an SBML string.

    Returns ``(Model, species_names, parameter_names)``.
    """
    import libsbml
    from poincare.reactions import Reactant, RateLaw
    from simbio import Compartment, Parameter, Volume, assign
    from simbio.core import amount, concentration, volume as volume_field

    document = libsbml.readSBMLFromString(sbml)
    model = document.getModel()
    if model is None:
        raise ValueError("libSBML could not parse the model")

    compartment_size = {
        c.getId(): (c.getSize() if c.isSetSize() else 1.0)
        for c in model.getListOfCompartments()
    }

    namespace: dict[str, Any] = {"__annotations__": {}}
    volume = next(iter(compartment_size.values()), 1.0)
    namespace["volume"] = volume_field(default=float(volume))
    namespace["__annotations__"]["volume"] = Volume

    eval_names: dict[str, Any] = {}
    species_objs: dict[str, Any] = {}
    species_names: list[str] = []
    for s in model.getListOfSpecies():
        sid = s.getId()
        if s.isSetInitialConcentration():
            initial, is_conc = s.getInitialConcentration(), True
        elif s.isSetInitialAmount():
            initial, is_conc = s.getInitialAmount(), False
        else:
            initial, is_conc = 0.0, True
        if initial is None or (isinstance(initial, float) and math.isnan(initial)):
            initial = 0.0
        obj = (concentration if is_conc else amount)(default=float(initial))
        namespace[sid] = obj
        namespace["__annotations__"][sid] = type(obj)
        species_objs[sid] = obj
        eval_names[sid] = obj
        species_names.append(sid)

    parameter_names: list[str] = []
    for p in model.getListOfParameters():
        pid = p.getId()
        value = p.getValue() if p.isSetValue() else 0.0
        param = assign(default=float(value))
        namespace[pid] = param
        namespace["__annotations__"][pid] = Parameter
        eval_names[pid] = param
        parameter_names.append(pid)

    for cid, size in compartment_size.items():
        eval_names.setdefault(cid, float(size))

    for r in model.getListOfReactions():
        kinetic = r.getKineticLaw()
        if kinetic is None:
            raise NotImplementedError(f"Reaction {r.getId()} has no kinetic law")
        local = {p.getId(): float(p.getValue()) for p in kinetic.getListOfParameters()}
        scope = {**eval_names, **local}
        rate_law = _eval_rate_law(libsbml.formulaToL3String(kinetic.getMath()), scope)
        reactants = [
            Reactant(species_objs[sr.getSpecies()], sr.getStoichiometry())
            for sr in r.getListOfReactants()
        ]
        products = [
            Reactant(species_objs[sr.getSpecies()], sr.getStoichiometry())
            for sr in r.getListOfProducts()
        ]
        namespace[f"_rxn_{r.getId()}"] = RateLaw(
            reactants=reactants, products=products, rate_law=rate_law
        )

    model_cls = type(name, (Compartment,), namespace)
    return model_cls, species_names, parameter_names


def model_from_antimony(antimony: str, *, model_name: str | None = None):
    """Compile Antimony -> SBML -> simbio model.

    Returns ``(Model, species_names, parameter_names)``.
    """
    sbml, name = antimony_to_sbml(antimony, model_name=model_name)
    return model_from_sbml(sbml, name=name)
