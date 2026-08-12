# Parameter provenance

This file separates literature/atomistic evidence from executable phase-field parameters.

## Status labels

- `DIRECT`: numerical value reported for the same material/observable and usable after unit conversion.
- `INFERRED`: obtained by fitting, conservation, or inversion from reported data; the inference must be recorded.
- `INFERRED_CONDITIONAL`: inferred only after fixing another non-identified parameter or geometry hypothesis.
- `SOURCE_MODEL_ASSUMPTION`: a numerical assumption used by a cited paper, not a direct measurement.
- `ATOMISTIC_REFERENCE`: atomistic/theoretical value reported in the cited literature and used only as a bound/reference.
- `ATOMISTIC`: obtained from DFT, MD or MLP calculations produced for this project.
- `PLACEHOLDER`: numerical-verification value only; it must not be interpreted physically.

## v0 numerical baseline

| Parameter group | Status | Source | Comment |
|---|---|---|---|
| crystal `a2,a4,a6` | PLACEHOLDER | numerical baseline | Dimensionless double-well verification only |
| OAF `a2,a4,a6` | PLACEHOLDER | numerical baseline | Weaker polar well used only to test phase contrast |
| MAF `a2,a4,a6` | PLACEHOLDER | numerical baseline | Paraelectric-like local response used only to test phase contrast |
| phase `eps_b` | PLACEHOLDER | numerical baseline | Background permittivity used to verify heterogeneous Poisson solver |
| `kappa` | PLACEHOLDER | numerical baseline | Constant isotropic gradient penalty |
| TDGL mobility/time step | PLACEHOLDER | numerical baseline | Numerical relaxation scale only |

No HHTT, free-volume or folded-boundary coupling is active in this baseline.

## v0.1 literature constraints

### Poled BOPVDF phase fractions and global response

| Quantity | Value | Status | Source / interpretation |
|---|---:|---|---|
| beta-crystal fraction | 0.52 | DIRECT | Huang et al., Nature Communications 12, 675 (2021), WAXD crystallinity |
| OAF fraction | ~0.28 | INFERRED | Rui et al., J. Mater. Chem. C 9, 894 (2021), three-phase hysteresis analysis |
| MAF remainder | 0.20 | INFERRED | closure `1 - f_crystal - f_OAF` |
| poled film `Ps` | 0.140 C m^-2 | DIRECT | Huang 2021 / Rui 2021 |
| beta-crystal theoretical `Ps` limit | 0.188 C m^-2 | ATOMISTIC_REFERENCE | DFT limit used by Huang 2021 |
| poled film `epsilon_r'` | 19.5 at 25 C, 10 Hz | DIRECT | Huang 2021 BDS |
| fresh film `epsilon_r'` | 11.5 at 25 C, 10 Hz | DIRECT | Huang 2021 BDS |
| poled high-field dynamic permittivity | 22.9 | DIRECT | Huang 2021, `dD/d(eps0 E)` at 300 MV/m |
| poled coercive field | 88 MV/m | DIRECT | Huang 2021 D-E loop |

With `f_crystal = 0.52`, `f_OAF = 0.28`, `f_MAF = 0.20`, `Ps,film = 0.140 C/m^2`, zero MAF spontaneous polarization, and the crystal assigned its upper bound `0.188 C/m^2`, the OAF local polarization must satisfy

`Ps,OAF >= (0.140 - 0.52*0.188)/0.28 = 0.15086 C/m^2`.

This is a lower bound, not a direct OAF measurement. Any lower crystal contribution requires a larger OAF contribution.

### Dielectric phase information

| Quantity | Value | Status | Source / interpretation |
|---|---:|---|---|
| BOPVDF amorphous-phase dielectric constant | ~21-22 at 25 C | DIRECT_REPORTED_RANGE | Yang et al., ACS Appl. Mater. Interfaces 7, 19894-19905 (2015) |
| crystal dielectric constant used in a later three-phase model | 3.0 | SOURCE_MODEL_ASSUMPTION | Rui et al., Macromolecules 55, 9705-9714 (2022), SI S2 |
| rigid OAF dielectric constant | set equal to crystal, 3.0 | SOURCE_MODEL_ASSUMPTION | same SI S2 |
| IAF dielectric constant | temperature-dependent, extrapolated from melt | SOURCE_MODEL_DERIVED | same SI S2 |
| mobile OAF dielectric constant | temperature-dependent inverse-model result | SOURCE_MODEL_DERIVED | same SI S2 |

The 2015 value `21-22` must not be silently assigned to a pure MAF region. That paper motivated a three-phase lamellar-crystal/oriented-interphase/amorphous model, whereas the later OAF/RAF literature separates amorphous subfractions more explicitly.

## Important implementation distinction

The TDGL variable `eps_b` is a **background** permittivity used in

`D = eps0 * eps_b * E + P`.

A measured low-frequency dielectric constant includes dipolar polarization that may later be represented explicitly by `P`. Therefore measured phase or film permittivity must not be copied directly into `eps_b`, or dipolar response can be counted twice.

v0.1 consequently keeps two layers separate:

1. a linear-dielectric inverse problem using effective local phase permittivities;
2. the TDGL background permittivity and local Landau response, which remain uncalibrated until local polarization-energy information is available.

## v0.1 morphology/orientation test

The first calibration constructs a laminate with the literature phase fractions and solves the periodic dielectric cell problem along both axes.

- `x`: field parallel to the idealized flat lamellae;
- `z`: field normal to the idealized flat lamellae.

The extended morphology test adds commensurate tilted laminates and periodic sinusoidal waviness.  One OAF effective permittivity is first inferred from the flat, field-parallel case and then held fixed during the geometry scan.  This is a transferability test; the inferred value remains `INFERRED`, not a direct material constant.

## v0.1 OAF relaxation scaffold

The dynamic scaffold uses

`epsilon*_OAF(omega) = epsilon_fast + Delta_epsilon/(1 + i omega tau)`.

At this stage neither `Delta_epsilon` nor `tau` is calibrated.  The 10 Hz film permittivity supplies only one real-valued constraint and cannot identify both quantities.  The code therefore fixes a series of `tau` hypotheses and infers the corresponding `Delta_epsilon`; all such outputs are labelled `INFERRED_CONDITIONAL`.

No scanned `tau` value is allowed to enter the physical TDGL configuration until frequency-resolved BDS or atomistic dynamics provides an independent constraint.

## Remaining quantities before a physical TDGL run

1. phase-resolved or atomistically derived free-energy curvature / switching barrier for crystal and OAF;
2. a defensible fast/background permittivity for each phase, separated from explicit dipolar polarization;
3. gradient coefficient or domain-wall/interphase length-scale calibration;
4. frequency dependence through OAF/RAF dipole mobility (`n`, `mu`, `g`, `tau_rot`);
5. HHTT content -> crystal/OAF/SC allocation;
6. free-volume radius/fraction -> frequency-dependent dielectric response.
