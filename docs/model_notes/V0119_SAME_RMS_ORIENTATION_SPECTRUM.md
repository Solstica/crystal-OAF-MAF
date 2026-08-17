# v0.1.19 Same-RMS orientation-spectrum sensitivity

## Question

v0.1.18 related the normalized local-field RMS approximately linearly to the geometry descriptor

```text
sqrt(<n_ND^2>).
```

That result was obtained inside one sinusoidal morphology family. A single second-moment descriptor can appear sufficient if every morphology differs only by the amplitude of the same spatial mode. v0.1.19 therefore holds `sqrt(<n_ND^2>)` fixed and changes the interface spectrum.

The numerical question is specific: at the same RMS lamellar-normal projection, do different spatial distributions of the same total orientational spread produce the same field RMS, p99 field and maximum field?

## Controlled quantities

The tested RMS values are

```text
0.10, 0.15, 0.20
```

all inside the v0.1.18 `PROJECT_IMAGE_DERIVED_CONSTRAINT` interval. Each spectrum is rescaled numerically until its realized RMS normal projection matches the target.

The beta/OAF/IAF thickness fractions, Rui-2022 generalized-Debye constitutive bank, frozen beta/OAF polarization amplitudes and uniform OAF source profile are unchanged. TDGL remains frozen.

## Geometry ensemble

Every interface is represented as a smooth periodic graph

```text
u = x/Lx + s * sum_k c_k sin(2*pi*k*z/Lz + phi_k).
```

The mean normal remains along MD. The perturbation depends only on the film-normal coordinate, so `du/dx` stays positive and the interface does not fold back on itself.

Five spectra are used:

- `single_k1`: v0.1.18-like one-mode control;
- `single_k2`: half the spatial wavelength at the same RMS normal spread;
- `two_mode_k1_k2`: two harmonic components;
- `smooth_multimode_A`;
- `smooth_multimode_B`.

Mode numbers, relative coefficients and phases are all `PLACEHOLDER_GEOMETRY_ENSEMBLE`. Huang 2021 constrains the preferred orientation and provides the rasterized SAXS pattern used for the RMS interval; it does not provide these real-space spectra.

## Reported diagnostics

For each morphology, v0.1.19 reports:

```text
realized sqrt(<n_ND^2>)
mean/max/95th/99th percentile |n_ND|
normalized corrector-field RMS
normalized field p90/p95/p99/max
phase-local field RMS in beta/OAF/IAF
normalized discrete divergence of the prescribed Pz source
Gauss residual
```

The field normalization remains

```text
E / (P_scale/eps0),
P_scale = max(|P_beta|, |P_OAF|).
```

The source-divergence diagnostic is expressed per normalized grid length. It is useful only for comparing project morphologies; the present grid spacing has not been mapped to a physical real-space discretization for this electrostatic audit.

## Numerical result

All 15 cases converged with Gauss residuals below `1e-7`. The target RMS normal projection was reproduced to better than `1.5e-12`, so differences below are not caused by mismatch in the controlled second moment.

At fixed `sqrt(<n_ND^2>)`, the five trial spectra produced the following spread, reported as `(max-min)/mean` across the spectrum ensemble:

| target RMS `n_ND` | field RMS spread | field p99 spread | field max spread | discrete source-divergence RMS spread |
|---:|---:|---:|---:|---:|
| 0.10 | 24.7% | 12.4% | 53.8% | 20.4% |
| 0.15 | 36.0% | 6.3% | 36.3% | 8.7% |
| 0.20 | 40.7% | 7.8% | 54.7% | 11.6% |

The single-`k=1` case reproduces the v0.1.18 geometry exactly: for target `0.10`, both versions use an internal coordinate amplitude of about `0.02267835`, generate the same phase map, and return the same normalized field RMS `0.00195630`. This cross-check removes a possible implementation difference between the v0.1.18 and v0.1.19 morphology generators.

For reference, the current re-run of v0.1.18 over its image-derived interval gives

```text
E_rms/(P_scale/eps0) = 0.00126 ... 0.00359
```

for `sqrt(<n_ND^2>) = 0.05 ... 0.20`, with a through-origin fit slope about `0.01834` and `R2_origin ~= 0.9946`. These are the current code results; older development notes or chat calculations that quoted a slope near `0.20` are superseded by the audited run retained in the CI artifact.

## Interpretation

The v0.1.19 ensemble does not support treating `sqrt(<n_ND^2>)` as a complete morphology descriptor. Field RMS changes by roughly 25-41% at fixed RMS orientation spread, and the maximum local field changes by roughly 36-55%. The spatial arrangement of interface slope therefore enters the computed field beyond its second moment in this discretized model.

The present result still requires one numerical check before it can be interpreted as a morphology effect. The phase map is pixelated on a `48 x 48` grid, and higher spatial harmonics change the discrete interface representation and the realized voxel fractions slightly. The same-RMS cases also show an 8.7-20.4% spread in the discrete polarization-source divergence. Some part of the field spread may therefore come from grid/interface aliasing rather than continuum spectral sensitivity.

The next calculation should repeat a representative same-RMS ensemble on progressively finer grids while holding the continuous harmonic functions and target RMS normal projection fixed. A converged nonzero spectrum-to-spectrum spread would justify retaining an additional morphology descriptor such as correlation length or curvature statistics. If the spread collapses with refinement, the v0.1.19 difference is primarily a pixelization artifact.

No result from v0.1.19 calibrates an experimental lamellar correlation spectrum or absolute local electric field. Those quantities still require compatible orientation/correlation data and phase-resolved remanent polarization.
