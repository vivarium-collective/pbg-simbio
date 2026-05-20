"""Regression tests for SBML/Antimony loader handling of rules & boundary species.

These mirror features common in real BioModels that the from-scratch
SBML->simbio bridge must handle: assignment rules (derived quantities used in
rate laws), boundary species (referenced but not consumed), and rate laws that
would reduce to 0/0 at build time if locals were inlined as float literals.
"""

import numpy as np
import pytest

from pbg_simbio.antimony_loader import model_from_antimony


def _solve(model):
    from simbio import Simulator

    return Simulator(model).solve(save_at=np.linspace(0, 5, 6))


def test_assignment_rule_is_not_a_state_species():
    """An assignment rule (`Total := A + B`) is inlined, not integrated."""
    ant = """model m
      species A = 10, B = 0;
      Total := A + B;
      k = 0.5;
      r: A -> B; k * A / Total;
    end"""
    model, species, _params = model_from_antimony(ant)
    assert "Total" not in species
    assert set(species) == {"A", "B"}
    res = _solve(model)
    assert "Total" not in res.data_vars


def test_boundary_species_not_consumed():
    """A boundary species ($S) drives a rate law but is not consumed."""
    ant = """model m
      species A = 10, B = 0;
      $S = 2;
      k = 0.5;
      r: A -> B; k * A * S;
    end"""
    model, species, _ = model_from_antimony(ant)
    res = _solve(model)
    # A is consumed, B produced; the closed A<->B mass is conserved
    a = float(np.asarray(res["A"])[-1])
    b = float(np.asarray(res["B"])[-1])
    assert a < 10.0 and b > 0.0
    assert abs((a + b) - 10.0) < 1e-6


def test_division_rate_law_with_assignment_denominator():
    """`X / CT` where CT is an assignment rule must not raise 0/0 at build."""
    ant = """model m
      species C2 = 1, M = 0;
      CT := C2 + M;
      k = 1.0;
      r: C2 -> M; k * M / CT;
    end"""
    # CT = C2 + M = 1 at t0 (nonzero), so M/CT is well defined — and crucially
    # building the RateLaw must not evaluate a literal 0.0/0.0.
    model, species, _ = model_from_antimony(ant)
    res = _solve(model)
    assert not np.any(np.isnan(np.asarray(res["M"])))
