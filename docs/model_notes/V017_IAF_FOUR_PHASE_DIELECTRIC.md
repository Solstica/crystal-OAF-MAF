# v0.1.7 — IAF permittivity and four-component dielectric cell

## 1. Source statement used

Rui et al. 2022 Supporting Information Section S2 writes the small-signal dielectric accounting as

`epsilon_c(T) = epsilon_cr eta_cr + epsilon_ROAF eta_ROAF(T) + epsilon_MOAF(T) eta_MOAF(T) + epsilon_IAF(T) eta_IAF`.

The same section assumes `epsilon_cr = epsilon_ROAF = 3.0` and states that `epsilon_IAF(T)` is obtained by extrapolating the molten-PVDF permittivity (above about 147 C) to lower temperature using main-text Figure 1B.

This stage implements that source-model instruction. It does not reinterpret the extrapolated value as a directly measured local IAF dielectric constant.

## 2. Figure 1B digitization

Main-text Figure 1B reports the molten-PVDF `epsilon_c(T)` used in the source analysis. Eleven traceable points from 147 to 197 C were digitized into

`data/literature/rui2022/figure_1b_melt_eps_digitized.csv`.

A linear least-squares fit gives approximately

`epsilon_IAF(T_C) = 17.4416 - 0.0432630 T_C`

with `R^2 ~= 0.99749` over the digitized melt interval.

The resulting source-model extrapolation in the project BDS window is approximately:

| T (C) | epsilon_IAF |
|---:|---:|
| -30 | 18.74 |
| 0 | 17.44 |
| 20 | 16.58 |
| 40 | 15.71 |

These values are labelled `SOURCE_MODEL_EXTRAPOLATION`.

## 3. Four dielectric components

The v0.1.7 small-signal state uses

- crystal: `epsilon = 3.0`, `SOURCE_MODEL_ASSUMPTION`;
- ROAF: `epsilon = 3.0`, `SOURCE_MODEL_ASSUMPTION`;
- MOAF: `epsilon_MOAF(T)` from Rui 2022 SI Figure S2, `SOURCE_MODEL_DERIVED`;
- IAF: melt-Figure-1B linear extrapolation, `SOURCE_MODEL_EXTRAPOLATION`.

The ROAF/MOAF fractions are supplied by the v0.1.6 project-regularized devitrification state. This avoids forcing the raw `x_RAF/x_MAF` calorimetric fractions into the rounded SI structural fractions.

## 4. Spatial allocation

The source constrains how much OAF becomes mobile but does not identify which OAF locations devitrify first. v0.1.7 therefore introduces two explicitly labelled sensitivity hypotheses:

1. `IAF-proximal-first`: OAF closer to the IAF/MAF side is assigned a higher mobility score;
2. `crystal-proximal-first`: reverse counterfactual used to bracket spatial sensitivity.

Both are `PROJECT_SPATIAL_HYPOTHESIS`, not literature-derived constitutive laws.

## 5. Electrostatic calculation

For a four-component map, the existing periodic heterogeneous dielectric cell solver evaluates

`div[epsilon_r(r) (E0 - grad psi)] = 0`

and returns the effective small-signal permittivity and local field distribution. The v0.1.7 driver evaluates:

- an ideal laminate with field parallel to lamellae as a source-model reference;
- a 45-degree periodic wavy morphology for local-field sensitivity;
- both spatial ROAF/MOAF allocation hypotheses at `T = -30, 0, 20, 40 C`.

## 6. Strict implementation boundary

The four values above are **not** TDGL background permittivities. They contain orientational/dipolar response. Directly assigning them to `eps_b` while also evolving an explicit polarization variable would double-count part of the dielectric response.

Thus v0.1.7 is a linear small-signal/local-field calibration layer. A later dynamic stage must separate fast/background polarization from explicit relaxational polarization before coupling these data to TDGL.
