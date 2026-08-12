# Baseline validation record

Date: 2026-08-12
Branch: `dev-phasefield`

## Scope

This record applies only to the dimensionless scalar-Pz numerical baseline. It does not validate physical PVDF parameters.

## Local checks completed before upload

- `tests/test_morphology.py`: passed.
  - crystal, OAF and MAF are all present in the generated lamellar morphology.
  - phase fractions sum to one.
- `tests/test_electrostatics.py`: passed.
  - a spatially uniform periodic polarization field produces zero depolarization potential/field, as required by the periodic electrostatic formulation.
- `scripts/run_baseline.py`: completed on a 64 x 64 grid.
  - one triangular -E -> +E -> -E cycle was executed;
  - final average dimensionless polarization: approximately `-1.04486`;
  - average-polarization range during the cycle: approximately `[-1.04486, 1.07225]`;
  - `state_final.npz`, `pe_history.csv`, `P_snapshots.npz` and metadata were produced locally.

## What this establishes

The current implementation is sufficient to test numerical coupling between:

1. a fixed crystal/OAF/MAF morphology;
2. phase-dependent local constitutive coefficients;
3. a heterogeneous `eps_b(r)` electrostatic solve;
4. scalar TDGL relaxation under a prescribed electric-field history.

## What remains unvalidated

- physical units and absolute field/polarization scales;
- phase-resolved PVDF Landau/constitutive parameters;
- OAF/MAF dielectric values;
- real lamellar length scales and phase fractions;
- HHTT, free-volume or folded-boundary couplings;
- vector polarization, elastic coupling and frequency-dependent dynamics.
