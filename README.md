# crystal-OAF-MAF

Phase-field development repository for semicrystalline PVDF/BOPVDF with an eventual **crystal / OAF / IAF** description and a DFT/MD/MLP -> continuum parameter bridge.

> Historical scripts may still use `MAF` as a phase identifier. In the Huang/Rui structural description the physical non-oriented amorphous region is treated here as IAF; RAF/MAF remain mobility-defined bookkeeping quantities and are not synonyms for OAF/IAF.

## Active branch: literature-reproduction

The current mainline is **source-first reproduction of published ferroelectric-polymer continuum models**. The previous v0.1.16-v0.1.19 frozen-source / waviness / RMS-orientation / harmonic-spectrum studies remain in history as numerical sensitivity experiments; they are not treated as reproduced material physics.

### Current literature anchors

1. **Guo et al., Nature Communications 15, 348 (2024)** — vector TDGL / polar-spiral phase-field benchmark. The one-axis Landau slice is reproduced exactly from S2/S3; the homogeneous strong/weak angular-anisotropy contrast in Fig. S16 is now reconstructed using source coefficients and a published cross-check of the contracted sixth-order convention.
2. **Ahluwalia et al., Physical Review B 78, 054110 (2008)** — now the primary multiscale-method anchor because the full article explicitly maps MD thermodynamics, domain walls, fluctuations, relaxation time and cell volume into LGD/TDGL quantities.
3. **Su et al., Nature Communications 13, 4867 (2022)** — beta-PVDF benchmark and explicit published reference for the conventional sixth-order Landau expansion and uniaxial PVDF bulk energy.
4. **Huang/Rui crystal-OAF-IAF evidence** enters only after the published polymer baselines identify which continuum coefficients can legitimately be replaced.
5. **DFT/MD/MLP** is then used to generate the atomistic observables needed to calibrate those coefficients, rather than supplying arbitrary phenomenological numbers.

Detailed evidence boundaries: `docs/model_notes/LITERATURE_REPRODUCTION_BASELINE.md`.

## Current executable checkpoints

### Guo 2024

```text
src/pvdf_pf/literature/guo2024.py
src/pvdf_pf/literature/guo2024_vector.py
scripts/reproduce_guo2024_landau_axis.py
scripts/reproduce_guo2024_s16_landau.py
tests/test_guo2024_landau_axis.py
tests/test_guo2024_vector_landau.py
docs/model_notes/GUO2024_REPRODUCTION.md
```

The branch distinguishes two claim levels:

- one-axis polynomial: exact direct transcription of the source-supported slice;
- full homogeneous S16 angular anisotropy: `SOURCE_CONSTRAINED_QUALITATIVE_S16_ANISOTROPY_REPRODUCTION`, not a pixel-exact reconstruction, because Guo does not report the plotting normalization used for the 3D energy surface.

The Guo peer-review record is also treated as part of the model provenance: the authors restrict the elastic treatment to the low-field linear regime and describe the phase-field rotational mechanism as semi-quantitative. This limitation is inherited here.

### Ahluwalia 2008

```text
src/pvdf_pf/literature/ahluwalia2008.py
scripts/reproduce_ahluwalia2008_multiscale.py
tests/test_ahluwalia2008_multiscale.py
docs/model_notes/AHLUWALIA2008_MULTISCALE_REPRODUCTION.md
```

The code preserves directly printed source values separately from quantities re-derived from rounded tables. In particular, re-applying Eqs. (2)-(3) to the rounded Table-I MD numbers does not regenerate Table II exactly; printed Table II therefore remains authoritative.

The future MLP bridge follows the paper's observable-based logic:

```text
DFT / MLP atomistic data
-> P(T), phase-energy differences, wall profiles/energies,
   elastic/dielectric response, fluctuations, relaxation times
-> calibrated LGD / gradient / electrostrictive / kinetic coefficients
-> phase-field predictions
```

## Next numerical stage

The next full solver milestone is **Guo 2024 vector TDGL under the published short-circuit boundary condition**, first with Landau + gradient + electrostatics, then with the source elastic/electrostrictive model under its stated low-field limitation. The target is the qualitative weak-anisotropy spiral texture of main Fig. 3e / Supplementary Fig. S16d before attempting field/stress-driven rotation in S17-S19.

No project-defined morphology spectrum, frozen beta/OAF source allocation, OAF decay profile or image-derived morphology constraint is permitted inside the literature-reproduction benchmark.

## Branches

- `main`: reviewed milestones.
- `dev-phasefield`: historical/experimental crystal-OAF-IAF development through v0.1.19.
- `literature-reproduction`: active source-first reproduction and multiscale-method reconstruction.

## Repository layout

```text
crystal-OAF-MAF/
├─ docs/model_notes/
├─ src/pvdf_pf/
│  ├─ literature/
│  ├─ calibration/
│  ├─ core/
│  ├─ morphology/
│  ├─ physics/
│  └─ solver/
├─ configs/
├─ scripts/
└─ tests/
```

**Rule: a quantity may enter the reproduction model only if its primary-source definition, value/state, or explicit derivation path is recorded. Missing source information is a stop condition, not an invitation to invent a parameter.**
