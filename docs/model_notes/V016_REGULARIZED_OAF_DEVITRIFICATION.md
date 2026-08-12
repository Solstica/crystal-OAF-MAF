# v0.1.6 - project-regularized OAF devitrification

## Source conflict that motivates this stage

Rui et al. main-text Figure 5B calculates mobility-defined `x_RAF(T)` and `x_MAF(T)` for melt-recrystallized PVDF using crystallinity `x_c=0.59`, hence `x_RAF+x_MAF≈0.41`.

The Supporting Information Section S2 then uses rounded structural assumptions `eta_cr=0.60`, `eta_IAF=0.20`, `eta_OAF≈0.20` and writes, above Tg,

```text
x_RAF = eta_ROAF
x_MAF = eta_MOAF + eta_IAF
```

Literal substitution is not physically admissible across the stated -30 to 40 C calculation range: at -30 C the article explicitly reports `x_MAF=0.166`, which would imply `eta_MOAF=-0.034` if `eta_IAF=0.20`. The raw Figure 5B fractions also close near 0.41 rather than the SI rounded 0.40.

The repository therefore records the inconsistency rather than clipping or silently renormalizing it.

## Project regularization

For the first non-negative temperature-driven spatial model, Figure 5B is used only to constrain the *progress* of amorphous devitrification over the same BDS window used by the source, -30 to 40 C.

Define

```text
q_RAF(T) = [x_RAF(-30)-x_RAF(T)] / [x_RAF(-30)-x_RAF(40)]
q_MAF(T) = [x_MAF(T)-x_MAF(-30)] / [x_MAF(40)-x_MAF(-30)]
q(T)     = 0.5 [q_RAF(T)+q_MAF(T)]
```

and retain the SI structural total `eta_OAF=0.20`:

```text
eta_ROAF(T) = 0.20 [1-q(T)]
eta_MOAF(T) = 0.20 q(T)
eta_IAF     = 0.20
eta_cr      = 0.60
```

This construction is labelled `PROJECT_REGULARIZATION_HYPOTHESIS`. It is not attributed to Rui et al. as their Eq. S4.

## What is source-constrained versus project-defined

Source-constrained:

- the Figure 5B RAF-loss/MAF-gain temperature trend;
- the explicit -30 C and 40 C anchors;
- `epsilon_MOAF(T)` from Figure S2;
- `lambda(T)` and `mu_r(T)` from Figure S3;
- the source paper's use of the same Figure 5B devitrification information when evaluating unpoled and poled films.

Project-defined:

- normalization of Figure 5B change to `q(T)`;
- partition of a fixed 0.20 OAF structural fraction with `q(T)`;
- any eventual spatial ordering of which OAF cells become MOAF first.

## Remaining gate before a dynamic local-field calculation

`epsilon_MOAF(T)` is an effective low-frequency/source-model dielectric target, not a TDGL background permittivity. The IAF dielectric curve used by the SI is obtained by extrapolation of molten PVDF and still needs to be transcribed/reconstructed. In addition, `mu_r(T)` is not automatically a Debye `tau(T)`.

Therefore v0.1.6 produces a temperature-dependent structural state and source kinetic descriptors, but does not yet claim a physically calibrated time-dependent TDGL run.
