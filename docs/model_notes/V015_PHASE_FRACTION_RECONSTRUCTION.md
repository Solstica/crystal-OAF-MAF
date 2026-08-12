# v0.1.5 — RAF/MAF to ROAF/MOAF phase-fraction reconstruction

## Source equations

Rui et al. 2022 Supporting Information Section S2 states

```text
eta_cr + eta_ROAF + eta_MOAF + eta_IAF = 1
```

and adopts the source-model values

```text
eta_cr  = 0.6
eta_IAF = 0.2
eta_OAF = eta_ROAF + eta_MOAF = 0.2
```

For temperatures above the glass-transition regime used in the source accounting,

```text
x_RAF = eta_ROAF
x_MAF = eta_MOAF + eta_IAF
```

Therefore, once Figure 5B of the main article is digitized,

```text
eta_ROAF(T) = x_RAF(T)
eta_MOAF(T) = x_MAF(T) - 0.2
```

with the closure checks

```text
x_RAF(T) + x_MAF(T) = 0.4
eta_ROAF(T) + eta_MOAF(T) = 0.2.
```

Below Tg, the SI states `x_MAF = 0`; the whole OAF is treated as rigid in the project structural partition.

## Guardrails

The following substitutions are explicitly disabled:

- `n(T) -> eta_MOAF(T)`;
- `epsilon_MOAF(T) -> eta_MOAF(T)`;
- a guessed sigmoid around Tg -> Figure 5B;
- mixing the source fixed fractions `0.6/0.2/0.2` with the separate BOPVDF `0.52/0.28/0.20` parameter set.

`n(T)` is useful as a devitrification diagnostic but is not a volume fraction. `epsilon_MOAF(T)` is a source-model-derived constitutive target and likewise does not determine the MOAF amount by itself.

## Current source gap

The uploaded Supporting Information contains the algebra above and states explicitly that the authors solve Eqs. S1-S4 using `x_RAF` and `x_MAF` values from **Figure 5B of the main article**. The SI itself does not contain Figure 5B.

The repository therefore keeps `data/literature/rui2022/figure_5b_maintext_digitized.csv` as a strict source-data gate. v0.1.5 can reconstruct the four fractions immediately after traceable Figure 5B points are supplied, but it will not fabricate those points from the SI curves.

## Next physical step after Figure 5B is available

At each temperature:

```text
x_RAF(T), x_MAF(T)
        ↓
eta_ROAF(T), eta_MOAF(T)
        ↓
crystal / ROAF / MOAF / IAF morphology state
        ↓
epsilon_MOAF(T), source kinetic descriptors
        ↓
local electrostatic / dynamic-polarization calculation
```

The final coupling to TDGL time remains disabled until a physical time-scale calibration or a validated mobility-to-relaxation relation is available.
