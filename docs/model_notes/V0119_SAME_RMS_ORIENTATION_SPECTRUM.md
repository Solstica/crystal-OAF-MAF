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

## Interpretation rule

Three outcomes are distinguishable.

If field RMS, p99 and maximum all remain nearly unchanged, `sqrt(<n_ND^2>)` is sufficient within the tested family.

If field RMS remains similar while p99/max differ, the RMS descriptor captures the average field scale but misses hotspot statistics. A second descriptor related to spatial correlation length, curvature or source-gradient concentration must then be retained.

If field RMS itself differs appreciably, the linear v0.1.18 relation is specific to the single-mode family and higher-order morphology information is required even for the bulk field magnitude.

No outcome from v0.1.19 calibrates an experimental lamellar correlation spectrum or absolute local electric field. Both require additional source information: raw orientation/correlation data for the morphology and phase-resolved remanent polarization for the source amplitude.
