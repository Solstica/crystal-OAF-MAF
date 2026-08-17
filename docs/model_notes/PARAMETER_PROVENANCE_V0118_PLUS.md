# Parameter provenance addendum: v0.1.18+

This addendum extends `PARAMETER_PROVENANCE.md` for geometry descriptors introduced after v0.1.17. It exists because these quantities are neither direct phase-field material constants nor quantities tabulated by the Huang/Rui sources.

## Additional status labels

- `PROJECT_IMAGE_DERIVED_CONSTRAINT`: a numerical interval estimated by this project from an identifiable published raster figure. The source figure and extraction limitation must be recorded. The interval must not be described as a directly reported measurement.
- `PLACEHOLDER_SENSITIVITY_BOUND`: an intentionally wider numerical bound used to test model sensitivity outside the source-constrained interval. It has no direct source basis.
- `PLACEHOLDER_GEOMETRY_ENSEMBLE`: a project-defined morphology spectrum/correlation pattern used to ask whether a lower-order geometry descriptor is sufficient. Its spatial spectrum must not be attributed to experiment.

## v0.1.18 lamellar-normal orientation descriptor

| Quantity | Value / range | Status | Source / interpretation |
|---|---:|---|---|
| Huang lamellar period | 11.8 nm | DIRECT_REPORTED_TEXT | Huang 2021 SI Supplementary Note 10 / Fig. 13 |
| beta layer thickness | 5.78 nm | DIRECT_REPORTED_TEXT | Huang 2021 SI Supplementary Fig. 13c |
| OAF thickness | 1.51 nm per crystal side | DIRECT_REPORTED_TEXT | Huang 2021 SI Supplementary Fig. 13c |
| IAF thickness | 3.00 nm | DIRECT_REPORTED_TEXT | Huang 2021 SI Supplementary Fig. 13c |
| preferred lamellar/crystal orientation | MD-dominant | DIRECT_QUALITATIVE | Huang 2021 main text/Methods and SI WAXD/SAXS orientation descriptions |
| raster-image first-lobe azimuthal FWHM estimate | about 9-20 deg | PROJECT_IMAGE_DERIVED_CONSTRAINT | project inspection of Huang 2021 Supplementary Fig. 13a; beamstop/cross, labels, background and pseudocolor prevent unique quantitative recovery |
| `sqrt(<n_ND^2>)` source-constrained scan | 0.05-0.20 | PROJECT_IMAGE_DERIVED_CONSTRAINT | deliberately widened around the image-level angular estimate; Huang et al. did not report these RMS values |
| `sqrt(<n_ND^2>)` upper scan | 0.25, 0.30 | PLACEHOLDER_SENSITIVITY_BOUND | sensitivity only |
| equivalent sinusoidal amplitude `A` | solved from each target RMS normal projection | INFERRED_CONDITIONAL | implementation variable for the single-mode morphology family, not a material observable |

The published source set available to the project contains no tabulated azimuthal `I(chi)` profile, numerical Herman orientation factor or directly reported azimuthal FWHM for Supplementary Fig. 13a. A later raw `I(chi)` dataset should replace the image-derived RMS interval without changing the electrostatic solver.

## v0.1.19 same-RMS morphology-spectrum test

v0.1.19 keeps `sqrt(<n_ND^2>)` fixed while varying the spatial spectrum of a smooth periodic interface. The trial spectra are `PLACEHOLDER_GEOMETRY_ENSEMBLE`. They answer a numerical identifiability question: whether two morphologies with the same RMS normal projection produce the same RMS field, peak-field statistics and phase-local field distributions.

Only the RMS normal-projection range can inherit the v0.1.18 image-derived constraint. Mode numbers, relative harmonic amplitudes, phases, random seeds and correlation lengths remain project geometry hypotheses unless a compatible SAXS azimuthal/spatial-correlation dataset becomes available.
