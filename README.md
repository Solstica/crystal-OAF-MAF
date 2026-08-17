# crystal-OAF-MAF

Phase-field development repository for semicrystalline PVDF with explicit **crystal / OAF / IAF** structural regions and separate mobility-based RAF/MAF bookkeeping where the cited source requires it.

> Historical scripts may still use `MAF` as a phase identifier. For the Huang/Rui three-region structural model, the physical region is named **IAF (isotropic amorphous fraction)**. RAF/MAF are mobility-defined fractions and are not synonyms for OAF/IAF.

## Physical scope

The executable development is deliberately stage-separated:

- fixed crystal / oriented amorphous fraction (OAF) / isotropic amorphous fraction (IAF) morphology;
- scalar polarization along the film-normal direction for the present electrostatic audits;
- phase-dependent electrostatics with spatially varying background permittivity `eps_b(r)`;
- experimentally constrained broadband amorphous relaxation represented by a positive generalized-Debye/Prony bank;
- prescribed frozen ferroelectric polarization as a separate source term before TDGL switching kinetics are activated;
- HHTT and free-volume fields reserved as explicit inputs, with couplings kept zero until calibrated from literature, experiment, MD, DFT or MLP.

The earlier inorganic relaxor code is treated only as a numerical reference. Its A/B-sublattice chemistry, PTO/STO/La endmembers, Vegard mapping and compensation-defect physics are not transferred into the PVDF model.

## Current executable stage

The active branch has reached **v0.1.19**.

- v0.1.12-v0.1.14: same-state BOPVDF broadband BDS -> Cole-Cole fit -> non-negative generalized-Debye/Prony time-domain representation.
- v0.1.15: self-consistent local-field coupling between the generalized-Debye polarization and heterogeneous Gauss electrostatics.
- v0.1.16: prescribed beta/OAF frozen ferroelectric source added to Gauss' law while TDGL remains frozen; Huang 2021 SAXS layer thicknesses are transferred as geometric ratios.
- v0.1.17: Huang-aligned lamellae are given controlled sinusoidal waviness to test how a local film-normal interface-normal component activates bound charge and redistributes the local field.
- v0.1.18: the internal waviness amplitude is replaced as the reported scan variable by `sqrt(<n_ND^2>)`, the RMS film-normal projection of the local lamellar/interface normal. A rasterized Huang Supplementary Fig. 13a estimate is recorded only as a project-level orientation constraint; larger values remain explicit sensitivity bounds.
- v0.1.19: several smooth periodic interface spectra are rescaled to the same `sqrt(<n_ND^2>)`. At fixed RMS orientation spread, the current `48 x 48` calculations show about 25-41% spread in field RMS and 36-55% spread in maximum field across spectra. These differences are not yet promoted to continuum morphology physics because interface pixelization and small voxel-fraction changes must be checked by grid refinement.

No current v0.1.16-v0.1.19 frozen-polarization amplitude is a phase-resolved measured remanent polarization. Absolute field magnitudes from these versions are sensitivity outputs until beta/OAF remanent-polarization allocation is independently constrained. The v0.1.18 `sqrt(<n_ND^2>) = 0.05-0.20` range is also not a directly reported Huang measurement; its provenance is recorded as `PROJECT_IMAGE_DERIVED_CONSTRAINT`. The real-space spectra introduced in v0.1.19 are `PLACEHOLDER_GEOMETRY_ENSEMBLE`.

## Branches

- `main`: reviewed milestones.
- `dev-phasefield`: active three-phase phase-field development.

## Repository layout

```text
crystal-OAF-MAF/
├─ docs/
│  └─ model_notes/
├─ src/
│  └─ pvdf_pf/
│     ├─ core/
│     ├─ morphology/
│     ├─ physics/
│     └─ solver/
├─ configs/
├─ scripts/
└─ tests/
```

## Development sequence

1. Verify crystal/OAF/IAF morphology and heterogeneous electrostatics numerically.
2. Calibrate same-state broadband dielectric relaxation without mixing it with TDGL time.
3. Add frozen ferroelectric sources and quantify orientation/morphology sensitivity.
4. Check grid convergence of the v0.1.19 same-RMS spectrum dependence; retain higher-order morphology descriptors only if the spread survives refinement.
5. Replace image-derived texture bounds with raw/tabulated azimuthal orientation data when available, and constrain phase-resolved remanent polarization independently.
6. Introduce field-dependent beta/OAF switching only after the static-source audit is closed.
7. Calibrate physical TDGL free-energy curvature, gradient coefficient and seconds-per-TDGL-time mapping.
8. Add HHTT and free-volume couplings with source-specific phase allocation.
9. Connect DFT / MD / MLP outputs to phase-field parameters.

**A numerical value is physical only when its provenance status and source/state are recorded in `docs/model_notes/PARAMETER_PROVENANCE.md` or the corresponding versioned provenance addendum.**
