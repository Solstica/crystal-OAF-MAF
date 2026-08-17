# v0.1.17 Huang-aligned lamellar-waviness sensitivity

## Question answered by this version

v0.1.16 established that a film-normal frozen polarization is tangent to an ideal flat Huang-type beta/OAF/IAF laminate whose normal lies along MD.  Such a scalar `P_FE,z` does not create a bound-charge field when it varies only across MD.

v0.1.17 asks a narrower question: if the laminate keeps the Huang orientation on average but its interfaces are locally wavy, how strongly does the local film-normal component of the interface normal activate bound charge and redistribute the self-consistent field?

The electrostatic constitutive equation is unchanged:

```text
D = eps0 * eps_b(r) * E + P_rel(r,t) + P_FE_frozen(r)
div(D) = 0
```

TDGL remains frozen and no ferroelectric switching threshold is used.

## Morphology family

The base periodic laminate is Huang-aligned:

```text
u0 = x/Lx
```

so the mean lamellar normal is along `x` (MD) and the frozen polarization is along `z` (film normal).  A periodic sinusoidal perturbation is added:

```text
u = x/Lx + A sin(2*pi*z/Lz)
```

where `A` is expressed as a fraction of one lamellar period.  The executable sweep uses

```text
A = 0, 0.01, 0.02, 0.05, 0.10, 0.15.
```

These amplitudes are `PLACEHOLDER_SENSITIVITY` values.  The current Huang/Rui source set reports stronger orientation along MD than TD and gives 2D WAXD/SAXS patterns, but the text/source tables available to this project do not provide a numerical azimuthal FWHM or a real-space lamellar-waviness distribution that can be transferred directly to `A`.

## Why waviness creates bound charge

For a discontinuity in polarization across an interface with unit normal `n`, the bound surface charge is

```text
sigma_b = Delta(P) dot n.
```

The prescribed ferroelectric source is along the film normal:

```text
P_FE = P_FE,z e_z.
```

The flat Huang-aligned limit has `n_z = 0`, so `sigma_b = 0` even when beta and OAF have different frozen-polarization magnitudes.  The wavy interface has a local slope and therefore `n_z != 0`; the same phase-polarization contrast then generates alternating bound charge along the interface.

For the phase coordinate used here,

```text
du/dz = (2*pi*A/Lz) cos(2*pi*z/Lz)
du/dx = 1/Lx
```

and the local film-normal projection of the interface normal is

```text
n_z = (du/dz) / sqrt[(du/dz)^2 + (du/dx)^2].
```

The runner reports RMS, mean-absolute and maximum `|n_z|` together with normalized local-field measures.  These geometry metrics allow later replacement of the arbitrary amplitude sweep by an experimentally calibrated orientation/waviness distribution without changing the electrostatic solver.

## Parameters held fixed

The phase thickness ratios remain those from Huang 2021 Supplementary Note 10:

```text
beta crystal = 5.78 nm
OAF total    = 3.02 nm (1.51 nm per side)
IAF          = 3.00 nm
period       = 11.8 nm
```

Their thickness ratios set geometric voxel fractions.  The WAXD beta crystallinity `0.52` remains a weight-fraction datum and is not substituted as a voxel fraction.

The Rui 2022 same-state generalized-Debye bank remains assigned to all non-crystal cells.  The beta and OAF frozen-source amplitudes remain the normalized v0.1.16 sensitivity values.  OAF frozen polarization is uniform in this version so that a change in field is attributable to morphology waviness rather than to an additional OAF-profile hypothesis.

## Interpretation boundary

This version can establish whether small local deviations from the ideal tangent-interface geometry activate a depolarization field and how the field scales with a geometric normal-projection metric.  It cannot state the actual field distribution in a measured BOPVDF specimen until at least one of the following is available:

- quantitative azimuthal orientation distribution/FWHM from WAXD or SAXS;
- real-space lamellar-normal distribution from a sufficiently resolved morphology measurement;
- another experimentally constrained descriptor that can be mapped to the local `n_z` distribution.

The absolute field magnitude also remains proportional to the underidentified phase-resolved frozen-polarization amplitudes.  The transferable output of v0.1.17 is therefore the normalized geometry dependence.

## Numerical checks

`tests/test_v0117_lamellar_waviness.py` requires:

- zero waviness reproduces the tangent-interface zero-field limit;
- finite waviness creates non-zero `z` variation in the frozen source and a finite self-consistent internal field;
- a larger waviness amplitude produces a larger RMS geometry-induced field than a small-amplitude case;
- the face-flux Gauss residual remains small.
