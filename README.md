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

The active branch has reached **v0.1.18**.

- v0.1.12-v0.1.14: same-state BOPVDF broadband BDS -> Cole-Cole fit -> non-negative generalized-Debye/Prony time-domain representation.
- v0.1.15: self-consistent local-field coupling between the generalized-Debye polarization and heterogeneous Gauss electrostatics.
- v0.1.16: prescribed beta/OAF frozen ferroelectric source added to Gauss' law while TDGL remains frozen; Huang 2021 SAXS layer thicknesses are transferred as geometric ratios.
- v0.1.17: Huang-aligned lamellae are given controlled sinusoidal waviness to test how a local film-normal interface-normal component activates bound charge and redistributes the local field.
- v0.1.18: the internal waviness amplitude is replaced as the reported scan variable by `sqrt(<n_ND^2>)`, the RMS film-normal projection of the local lamellar/interface normal. A rasterized Huang Supplementary Fig. 13a estimate is recorded only as a project-level orientation constraint; larger values remain explicit sensitivity bounds.

No current v0.1.16-v0.1.18 frozen-polarization amplitude is a phase-resolved measured remanent polarization. Absolute field magnitudes from these versions are sensitivity outputs until beta/OAF remanent-polarization allocation is independently constrained. The v0.1.18 `sqrt(<n_ND^2>) = 0.05-0.20` range is also not a directly reported Huang measurement; its provenance is recorded as `PROJECT_IMAGE_DERIVED_CONSTRAINT`.

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
4. Replace image-derived texture bounds with raw/tabulated azimuthal orientation data when available, and constrain phase-resolved remanent polarization independently.
5. Introduce field-dependent beta/OAF switching only after the static-source audit is closed.
6. Calibrate physical TDGL free-energy curvature, gradient coefficient and seconds-per-TDGL-time mapping.
7. Add HHTT and free-volume couplings with source-specific phase allocation.
8. Connect DFT / MD / MLP outputs to phase-field parameters.

**A numerical value is physical only when its provenance status and source/state are recorded in `docs/model_notes/PARAMETER_PROVENANCE.md`.**
