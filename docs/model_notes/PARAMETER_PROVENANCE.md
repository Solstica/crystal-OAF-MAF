# Parameter provenance

This file separates literature/atomistic evidence from executable phase-field parameters.

## Status labels

- `DIRECT`: numerical value reported for the same material/observable and usable after unit conversion.
- `DIRECT_TABULATED`: numerical value explicitly tabulated by the primary source/SI and transcribed with traceable table/figure context.
- `DIGITIZED_SOURCE`: numerical point read from a published source figure; figure/panel, sample state, unit and digitization provenance must be retained.
- `DIRECT_QUALITATIVE`: qualitative trend stated by the primary source; no numerical curve/value has been transferred into the executable model.
- `INFERRED`: obtained by fitting, conservation, or inversion from reported data; the inference must be recorded.
- `INFERRED_CONDITIONAL`: inferred only after fixing another non-identified parameter or geometry hypothesis.
- `INFERRED_SPECTRUM`: fitted from a supplied/digitized broadband spectrum; valid only for the stated fit model and spectral decomposition.
- `SOURCE_MODEL_ASSUMPTION`: a numerical assumption used by a cited paper, not a direct measurement.
- `SOURCE_MODEL_DERIVED`: a quantity calculated inside the cited paper's model from measured and assumed inputs; not a direct local measurement.
- `ATOMISTIC_REFERENCE`: atomistic/theoretical value reported in the cited literature and used only as a bound/reference.
- `ATOMISTIC`: obtained from DFT, MD or MLP calculations produced for this project.
- `NOT_YET_TRANSCRIBED`: the source/SI is known to contain the quantity, but traceable numerical values have not yet been entered.
- `UNRESOLVED_MAPPING`: the observable is associated with a related physical region/state, but its mapping to the project variable has not been validated.
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
| BOPVDF amorphous-phase dielectric constant | ~21-22 at 25 C | DIRECT | Yang et al., ACS Appl. Mater. Interfaces 7, 19894-19905 (2015); do not equate this automatically with pure MAF |
| crystal dielectric constant used in later BOPVDF model | 3.0 | SOURCE_MODEL_ASSUMPTION | Rui et al., Macromolecules 55, 9705-9714 (2022), SI S2 |
| rigid OAF dielectric constant set equal to crystal | 3.0 | SOURCE_MODEL_ASSUMPTION | same SI S2 |
| IAF dielectric constant | temperature-dependent, extrapolated from melt | SOURCE_MODEL_DERIVED | same source-model construction; numerical curve not yet transcribed |
| mobile OAF dielectric constant | temperature-dependent model inversion | SOURCE_MODEL_DERIVED | same SI S2; Figure S2 contains the resulting curve |

The 2015 value `21-22` must not be silently assigned to a pure MAF region. The BOPVDF literature uses different amorphous/interphase partitions across papers, so every transfer must preserve the source definition.

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

The extended morphology test adds commensurate tilted laminates and periodic sinusoidal waviness. One OAF effective permittivity is first inferred from the flat, field-parallel case and then held fixed during the geometry scan. This is a transferability test; the inferred value remains `INFERRED`, not a direct material constant.

## v0.1 OAF relaxation scaffold

The first dynamic scaffold uses

`epsilon*_OAF(omega) = epsilon_fast + Delta_epsilon/(1 + i omega tau)`.

At this stage neither `Delta_epsilon` nor `tau` is calibrated. The 10 Hz film permittivity supplies only one real-valued constraint and cannot identify both quantities. The code therefore fixes a series of `tau` hypotheses and infers the corresponding `Delta_epsilon`; all such outputs are labelled `INFERRED_CONDITIONAL`.

No scanned `tau` value is allowed to enter the physical TDGL configuration until frequency-resolved BDS or atomistic dynamics provides an independent constraint.

## v0.1.2 BDS dynamics bridge

Primary source: Rui et al., *Macromolecules* 55, 9705-9714 (2022), DOI `10.1021/acs.macromol.2c01110`.

| Quantity / relation | Current status | Allowed use in code |
|---|---|---|
| active-dipole concentration increases strongly from -30 to 40 C | DIRECT_QUALITATIVE | monotonic/trend constraint only; no numerical `n(T)` table yet |
| calculated dipole moment is higher in poled than unpoled BOPVDF | DIRECT_QUALITATIVE | comparison constraint only |
| Kirkwood-Frohlich `g` is higher in poled than unpoled BOPVDF | DIRECT_QUALITATIVE | comparison constraint only |
| dipole-dipole interaction increases with temperature | DIRECT_QUALITATIVE | trend constraint only |
| rotational dipole mobility increases by >4 orders of magnitude over -30 to 40 C | DIRECT_QUALITATIVE | order-of-magnitude trend constraint only |
| Debye `eps_inf, Delta_eps, tau` fitted by `fit_bds_csv.py` | INFERRED_SPECTRUM | spectrum-level constitutive fit only; retain sample state, T, frequency window and fit model |
| Kirkwood-Frohlich `g*mu^2` | INFERRED_CONDITIONAL | requires independently supplied active-dipole number density |
| time-domain auxiliary `P_rel` | numerical bridge | may be tested with physical seconds; not coupled to TDGL yet |
| TDGL seconds-per-time-unit mapping | uncalibrated | dynamic coupling disabled |

The project uses the complex-permittivity convention `epsilon* = epsilon' - i epsilon''` with positive reported loss `epsilon_loss = -Im(epsilon*)`.

The time-domain auxiliary response obeys

`tau dP_rel/dt + P_rel = eps0 Delta_epsilon E_local`.

Its exact exponential step is exact for a zero-order-hold electric field during each numerical step. The corresponding discrete harmonic transfer function is tracked separately from the continuum Debye expression so that time-discretization phase lag is not mistaken for material relaxation.

## v0.1.3 source-driven OAF subpartition

The ACS Supporting Information explicitly states that, above the glass-transition regime where devitrification occurs, the OAF contains a rigid part (ROAF) and a mobile part (MOAF). The source model writes

`epsilon_film(T) = epsilon_cr*eta_cr + epsilon_ROAF*eta_ROAF(T) + epsilon_MOAF(T)*eta_MOAF(T) + epsilon_IAF(T)*eta_IAF`.

Repository mapping:

| Project/source quantity | Status | Rule |
|---|---|---|
| total project OAF -> `ROAF + MOAF` | DIRECT_QUALITATIVE / source definition | allowed as an internal OAF subpartition |
| project MAF -> source IAF | UNRESOLVED_MAPPING | provisional only; do not claim identity |
| Figure S1 `n(T), m_d(T), g(T)` | NOT_YET_TRANSCRIBED | enter only as `DIGITIZED_SOURCE` or `DIRECT_TABULATED` |
| Figure S2 `epsilon_MOAF(T)` | NOT_YET_TRANSCRIBED | source-model-derived curve; preserve this status after digitization |
| SI `lambda(T), mu_r(T)` | NOT_YET_TRANSCRIBED | enter with full figure/panel/unit provenance |
| `epsilon_cr = 3.0` | SOURCE_MODEL_ASSUMPTION | may reproduce source algebra; not a TDGL background constant by default |
| `epsilon_ROAF = epsilon_cr` | SOURCE_MODEL_ASSUMPTION | same restriction |

v0.1.3 introduces separate MOAF and IAF linear-relaxation channels, but their physical `tau(T)` and `Delta_epsilon(T)` remain disabled until traceable source values are supplied. If the real BDS curves show broad/non-Debye relaxation, a relaxation-time distribution or HN/Prony representation will replace the single-Debye channel rather than forcing a narrow fit.

## Remaining quantities before a physical TDGL run

1. phase-resolved or atomistically derived free-energy curvature / switching barrier for crystal and OAF;
2. a defensible fast/background permittivity for each phase, separated from explicit dipolar polarization;
3. gradient coefficient or domain-wall/interphase length-scale calibration;
4. traceable frequency/temperature-resolved MOAF/IAF dynamics;
5. a validated project-MAF/source-IAF mapping;
6. a physical seconds-per-TDGL-time calibration before coupling BDS `tau` to TDGL;
7. HHTT content -> crystal/OAF/SC allocation;
8. free-volume radius/fraction -> frequency-dependent dielectric response.
