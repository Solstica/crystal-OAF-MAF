# crystal-OAF-MAF

Phase-field development repository for semicrystalline PVDF with explicit **crystal / OAF / MAF** structure.

## Physical scope

The first executable model is intentionally narrow:

- fixed crystal / oriented amorphous fraction (OAF) / mobile amorphous fraction (MAF) morphology;
- scalar polarization field `Pz(z, x)` for the first validation stage;
- phase-dependent local Landau free energy;
- gradient energy and TDGL relaxation;
- electrostatics with spatially varying background permittivity `eps_b(r)`;
- triangular external-field schedules for P-E scans;
- HHTT and free-volume fields reserved as explicit inputs, with their couplings kept zero until calibrated from literature, experiment, MD or DFT.

The earlier inorganic relaxor code is treated only as a numerical reference. Its A/B-sublattice chemistry, PTO/STO/La endmembers, Vegard mapping and compensation-defect physics are not transferred into the PVDF model.

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

1. Verify the crystal/OAF/MAF baseline numerically.
2. Calibrate phase-local dielectric and polarization parameters.
3. Add HHTT through experimentally supported phase-allocation / local-energy couplings.
4. Add free-volume holes through dielectric / relaxation couplings.
5. Extend to vector polarization and local chain orientation.
6. Add folded-boundary / internal-strain coupling.
7. Connect DFT / MD / MLP outputs to phase-field parameters.

**Current parameter sets are dimensionless verification parameters unless a source is explicitly recorded in `docs/model_notes/PARAMETER_PROVENANCE.md`.**
