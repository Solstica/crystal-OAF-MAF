# Parameter provenance

This file separates literature/atomistic evidence from executable phase-field parameters.

## Status labels

- `DIRECT`: numerical value reported for the same material/observable and usable after unit conversion.
- `DIRECT_TABULATED`: numerical value explicitly tabulated by the primary source/SI and transcribed with traceable table/figure context.
- `DIRECT_REPORTED_TEXT`: numerical value explicitly stated in the prose/caption of the primary source. This is preferred over plot digitization when both refer to the same temperature/state.
- `DIGITIZED_SOURCE`: numerical point read from a published source figure; figure/panel, sample state, unit and digitization provenance must be retained.
- `DIRECT_QUALITATIVE`: qualitative trend stated by the primary source; no numerical curve/value has been transferred into the executable model.
- `INFERRED`: obtained by fitting, conservation, or inversion from reported data; the inference must be recorded.
- `INFERRED_CONDITIONAL`: inferred only after fixing another non-identified parameter or geometry hypothesis.
- `INFERRED_SPECTRUM`: fitted from a supplied/digitized broadband spectrum; valid only for the stated fit model and spectral decomposition.
- `SOURCE_MODEL_ASSUMPTION`: a numerical assumption used by a cited paper, not a direct measurement.
- `SOURCE_MODEL_DERIVED`: a quantity calculated inside the cited paper's model from measured and assumed inputs; not a direct local measurement.
- `SOURCE_MODEL_EXTRAPOLATION`: a quantity obtained by carrying out an extrapolation explicitly prescribed by the cited source model; it remains model-dependent and is not a direct local measurement.
- `PROJECT_REGULARIZATION_HYPOTHESIS`: a project-defined mapping introduced to reconcile source definitions or enforce physical admissibility; it must not be attributed to the source paper.
- `PROJECT_SPATIAL_HYPOTHESIS`: a project-defined spatial allocation rule used only for sensitivity analysis until an atomistic or experimental localization rule is available.
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
| IAF dielectric constant | `epsilon_IAF(T) ~= 17.4416 - 0.0432630 T_C` over project use window -30 to 40 C | SOURCE_MODEL_EXTRAPOLATION | SI S2 explicitly prescribes extrapolation from molten-PVDF main-text Figure 1B; Figure 1B was digitized in v0.1.7 |
| mobile OAF dielectric constant | temperature-dependent model inversion | SOURCE_MODEL_DERIVED | same SI S2; Figure S2 digitized in v0.1.3 |

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

The project uses the complex-permittivity convention `epsilon* = epsilon' - i epsilon''` with positive reported loss `epsilon_loss = -Im(epsilon*)`.

The time-domain auxiliary response obeys

`tau dP_rel/dt + P_rel = eps0 Delta_epsilon E_local`.

Its exact exponential step is exact for a zero-order-hold electric field during each numerical step. The corresponding discrete harmonic transfer function is tracked separately from the continuum Debye expression so that time-discretization phase lag is not mistaken for material relaxation.

## v0.1.3-v0.1.5 Rui 2022 source-driven OAF/RAF state

The SI defines the OAF subpartition `OAF = ROAF + MOAF` for its dielectric model and supplies Figure S1-S3 curves. The main article Figure 5B supplies calorimetrically determined `x_RAF(T)` and `x_MAF(T)` for melt-recrystallized BOPVDF with crystallinity `x_c = 0.59`.

Important distinction:

- Main-text Figure 5B has `x_RAF + x_MAF ~= 0.41`, because the measured crystallinity used there is `x_c = 0.59`.
- SI Section S2 subsequently approximates `eta_cr = 0.6`, `eta_IAF = 0.2`, and `eta_OAF ~= 0.2` for the dielectric inversion.
- Therefore a code check of `x_RAF + x_MAF = 0.40` against raw Figure 5B data is incorrect. The raw calorimetric state and the SI rounded dielectric state must remain separate.

Explicit main-text anchors are stored as `DIRECT_REPORTED_TEXT`:

- below -55 C: `x_RAF = 0.40`, `x_MAF = 0.008`;
- `Tg = -45.2 C`: `x_RAF = 0.331`, `x_MAF = 0.079`;
- -30 C: `x_RAF = 0.244`, `x_MAF = 0.166`;
- 40 C: `x_RAF = 0.014`, `x_MAF = 0.396`.

Other Figure 5B marker values are `DIGITIZED_SOURCE`. For the project OAF subpartition, Figure 5B is used as a temperature-dependent mobility constraint, while the exact ROAF/MOAF conversion must explicitly select either the raw calorimetric `x_c=0.59` state or the SI rounded `eta_cr=0.6` dielectric-model state.

The project does **not** equate RAF with OAF or MAF with IAF: RAF/MAF are mobility-defined fractions; OAF/IAF are structure-defined fractions.

## v0.1.6 project regularization

Because literal substitution of the raw Figure 5B fractions into the SI rounded structural fractions produces a negative implied `eta_MOAF` near -30 C, the project does not clip or silently renormalize the source data. Instead it uses normalized RAF loss and MAF gain only as a devitrification progress `q(T)` and partitions the fixed structural OAF total as

`eta_ROAF = 0.20 * (1-q)` and `eta_MOAF = 0.20 * q`.

This mapping is `PROJECT_REGULARIZATION_HYPOTHESIS` and is not attributed to Rui et al.

## v0.1.7 IAF extrapolation and four-component dielectric cell

Main-text Figure 1B molten-PVDF values from 147 to 197 C were digitized and linearly fitted because SI S2 explicitly states that `epsilon_IAF(T)` is obtained by extrapolating the melt permittivity to lower temperatures. The executable fit is approximately

`epsilon_IAF(T_C) = 17.4416 - 0.0432630*T_C`, `R^2 ~= 0.99749`.

Representative values are approximately 18.74 (-30 C), 17.44 (0 C), 16.58 (20 C), and 15.71 (40 C). These are `SOURCE_MODEL_EXTRAPOLATION`.

The v0.1.7 four-component small-signal cell uses crystal/ROAF source assumptions, source-derived `epsilon_MOAF(T)`, and source-model-extrapolated `epsilon_IAF(T)`. Spatial placement of MOAF within OAF remains unresolved by the source, so `IAF-proximal-first` and `crystal-proximal-first` are recorded only as `PROJECT_SPATIAL_HYPOTHESIS` sensitivity bounds.

None of these small-signal effective permittivities may be copied directly into TDGL `eps_b` until fast/background and explicit relaxational polarization contributions are separated.

## Remaining quantities before a physical TDGL run

1. phase-resolved or atomistically derived free-energy curvature / switching barrier for crystal and OAF;
2. a defensible fast/background permittivity for each phase, separated from explicit dipolar polarization;
3. gradient coefficient or domain-wall/interphase length-scale calibration;
4. traceable frequency/temperature-resolved MOAF/IAF relaxation times;
5. a validated spatial rule for which OAF regions devitrify first;
6. a physical seconds-per-TDGL-time calibration before coupling BDS time scales to TDGL;
7. HHTT content -> crystal/OAF/SC allocation;
8. free-volume radius/fraction -> frequency-dependent dielectric response.
