# v0.1.8 — measured-film dielectric closure audit

## 1. Why this gate is required

v0.1.7 assembled a four-component small-signal dielectric state from Rui et al. 2022:

- `epsilon_crystal = 3.0`;
- `epsilon_ROAF = 3.0`;
- `epsilon_MOAF(T)` from SI Figure S2;
- `epsilon_IAF(T)` from the SI-prescribed extrapolation of main-text Figure 1B.

Before coupling this state to a more complicated phase-field/dynamic model, it must reproduce the simplest macroscopic observable used by the source: the measured BOPVDF film `epsilon_c(T)` in main-text Figure 3A under ideal parallel-capacitor accounting.

No parameter is fitted in this audit.

## 2. Figure 3A digitization

The unpoled and poled Figure 3A marker values at `T = -30,-20,...,40 C` were digitized into

`data/literature/rui2022/figure_3a_film_eps_digitized.csv`.

Representative poled values are approximately:

| T (C) | measured epsilon_c |
|---:|---:|
| -30 | 14.412 |
| 0 | 16.238 |
| 20 | 17.120 |
| 40 | 18.226 |

The unpoled curve rises from about 11.925 at -30 C to 13.749 at 40 C.

## 3. Project-regularized closure test

Using the v0.1.6 structural partition

`eta_cr = 0.60`, `eta_IAF = 0.20`, `eta_ROAF + eta_MOAF = 0.20`,

and the v0.1.7 local dielectric terms, the ideal parallel prediction is

`epsilon_pred = eta_cr epsilon_cr + eta_ROAF epsilon_ROAF + eta_MOAF epsilon_MOAF + eta_IAF epsilon_IAF`.

The closure is not quantitative.

### Poled BOPVDF

| T (C) | measured | project prediction | relative error | project mean OAF epsilon | OAF mean required for closure |
|---:|---:|---:|---:|---:|---:|
| -30 | 14.412 | 6.148 | -57.3% | 3.00 | 44.32 |
| 0 | 16.238 | 7.757 | -52.2% | 12.34 | 54.75 |
| 20 | 17.120 | 9.320 | -45.6% | 21.02 | 60.02 |
| 40 | 18.226 | 11.224 | -38.4% | 31.41 | 66.42 |

Across all eight temperatures:

- RMSE is about `8.08` in relative-permittivity units;
- mean relative error is about `-49.6%`;
- maximum absolute relative error is about `57.3%`;
- the OAF-mean permittivity required for exact film closure exceeds the current project OAF mean by about `40.3` on average.

### Unpoled BOPVDF

The same audit gives approximately:

- RMSE `5.68`;
- mean relative error `-43.5%`;
- maximum absolute relative error `49.2%`;
- mean gap between required and current OAF-mean permittivity `28.3`.

Thus the discrepancy is present for both electrical states and is too large to attribute to marker-digitization uncertainty.

## 4. Literal SI Eq. S4 audit

The project also evaluates the algebra written in SI Eq. S4 without using it as a physical state:

`eta_ROAF = x_RAF`,

`eta_MOAF = x_MAF - 0.20`.

At -30, -20 and -10 C this gives a negative `eta_MOAF`, so the literal state is physically inadmissible and is not evaluated as a dielectric mixture.

At 0 C and above, where `eta_MOAF >= 0`, the raw Figure 5B fractions produce a total fraction of about `1.01` because the main-text calorimetry used `x_c=0.59` whereas SI S2 rounded the crystal fraction to `0.60`. Even in this admissible-temperature subset, the calculated film permittivity remains far below Figure 3A. For poled BOPVDF at 0 C, for example, the literal-SI accounting gives about `6.71` versus measured `16.24`.

Therefore the closure problem is not created only by the v0.1.6 regularization.

## 5. What this result does and does not establish

It establishes that the set of quantities currently transcribed from the 2022 main text/SI, when combined according to the equations and rounded fractions currently represented in the project, does **not** reproduce the measured Figure 3A film permittivity.

It does not by itself identify the reason. Possible items requiring source reconciliation include:

1. the exact low-temperature `epsilon_IAF(T)` construction intended through main-text Figure 1B plus ref. S1 Figure 9A;
2. the relationship between the 2022 calorimetric RAF/MAF fractions and the structural OAF/IAF fractions taken from the earlier BOPVDF work;
3. whether an additional normalization or source-specific fraction convention was used in the SI calculation but not made explicit in Eqs. S1-S4;
4. whether the plotted Figure S2 quantity should be interpreted with a more specific source convention than the current project transfer.

The discrepancy is retained as `SOURCE_MODEL_CLOSURE_UNRESOLVED`. No source value is modified to remove it.

## 6. Consequence for the project

The v0.1.7 four-component dielectric model remains useful for morphology/local-field sensitivity tests, but it is **not promoted to a quantitatively calibrated material model**.

The next source-reconciliation task should inspect Rui et al. 2021 (the SI's ref. S1), especially the definition and construction associated with its Figure 9A and the OAF/IAF fractions. Dynamic TDGL coupling remains blocked until this macroscopic closure issue is explained.
