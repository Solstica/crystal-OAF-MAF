# Huang 2021 geometry/orientation source note

Source used by the project: Huang et al., *Nature Communications* 12, 675 (2021), together with its Supplementary Information.

## Directly reported quantities

Supplementary Note 10 / Supplementary Fig. 13 reports the poled BOPVDF three-region lamellar geometry:

- lamellar period: 11.8 nm;
- beta-PVDF layer: 5.78 nm;
- OAF: 1.51 nm on each crystal side, 3.02 nm total;
- IAF: 3.00 nm.

Supplementary Fig. 13a is a 2D edge-on SAXS pattern. The X-ray beam is along TD, MD is vertical/upward in the figure, and ND is horizontal. The authors state that the 1D SAXS curve is integrated along MD. Supplementary Fig. 1 additionally reports stronger crystal orientation along MD than TD.

The published SI does not provide a tabulated azimuthal intensity profile, Herman orientation factor, numerical azimuthal FWHM, or a real-space interface-waviness distribution that can be transferred directly to the phase-field morphology.

## Project image-derived orientation constraint for v0.1.18

The rasterized Supplementary Fig. 13a was inspected only to constrain the order of magnitude of the lamellar-normal spread. Repeated annular/ridge estimates around the first-order SAXS lobe gave an image-dependent azimuthal FWHM of roughly 9-20 degrees. Interpreting that width as a symmetric Gaussian-like spread gives a characteristic angular standard deviation of about 3.8-8.5 degrees and a small-angle film-normal normal-projection scale of order 0.07-0.15.

Because the figure contains beamstop/cross features, labels, background intensity and pseudocolor compression, v0.1.18 does not promote those numbers to a measured orientation distribution. The executable scan therefore uses a deliberately widened source-constrained interval

```text
sqrt(<n_ND^2>) = 0.05 ... 0.20
```

where `n_ND` is the film-normal component of the local lamellar/interface normal in the present 2D mapping. Values 0.25 and 0.30 are retained only as sensitivity bounds.

Provenance status:

```text
0.05-0.20 : PROJECT_IMAGE_DERIVED_CONSTRAINT
0.25-0.30 : PLACEHOLDER_SENSITIVITY_BOUND
```

This note must not be cited as if Huang et al. reported `sqrt(<n_ND^2>)` directly. A later raw `I(chi)` curve or author-provided azimuthal FWHM should replace the image-derived interval without changing the electrostatic solver.
