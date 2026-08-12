# Parameter provenance

This file separates literature/atomistic evidence from executable phase-field parameters.

## Status labels

- `DIRECT`: numerical value reported for the same material/observable and usable after unit conversion.
- `INFERRED`: obtained by fitting or inversion from reported data; the fitting procedure must be recorded.
- `ATOMISTIC`: obtained from DFT, MD or MLP calculations produced for this project.
- `PLACEHOLDER`: numerical-verification value only; it must not be interpreted physically.

## Current executable baseline

| Parameter group | Status | Source | Comment |
|---|---|---|---|
| crystal `a2,a4,a6` | PLACEHOLDER | numerical baseline | Dimensionless double-well verification only |
| OAF `a2,a4,a6` | PLACEHOLDER | numerical baseline | Weaker polar well used only to test phase contrast |
| MAF `a2,a4,a6` | PLACEHOLDER | numerical baseline | Paraelectric-like local response used only to test phase contrast |
| phase `eps_b` | PLACEHOLDER | numerical baseline | Used to verify heterogeneous Poisson solver |
| `kappa` | PLACEHOLDER | numerical baseline | Constant isotropic gradient penalty |
| TDGL mobility/time step | PLACEHOLDER | numerical baseline | Numerical relaxation scale only |

No HHTT, free-volume or folded-boundary coupling is active in the baseline.

## Next calibration targets

1. phase fractions and OAF thickness / distribution;
2. crystal, OAF and MAF dielectric response;
3. phase-resolved polarization response or a constrained fit to bulk P-E data;
4. HHTT content -> crystal/OAF/SC allocation;
5. free-volume radius/fraction -> frequency-dependent dielectric response.
