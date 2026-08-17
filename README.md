# crystal-OAF-MAF

Phase-field development repository for semicrystalline PVDF/BOPVDF with an eventual **crystal / OAF / IAF** description and a DFT/MD/MLP -> continuum parameter bridge.

> Historical scripts may still use `MAF` as a phase identifier. In the Huang/Rui structural description the physical non-oriented amorphous region is treated here as IAF; RAF/MAF remain mobility-defined bookkeeping quantities and are not synonyms for OAF/IAF.

## Active branch: literature-reproduction

The current mainline is **source-first reproduction of published ferroelectric-polymer continuum models**. The previous v0.1.16-v0.1.19 frozen-source / waviness / RMS-orientation / harmonic-spectrum studies remain in history as numerical sensitivity experiments; they are not treated as reproduced material physics.

### Current literature anchors

1. **Guo et al., Nature Communications 15, 348 (2024)** — vector TDGL / polar-spiral target. The one-axis Landau slice and homogeneous strong/weak angular-anisotropy contrast are reconstructed from the source. Exact full spiral reproduction is currently source-gated because the available paper/SI do not numerically specify all PDE inputs such as `kappa_ij`, `L`, FEM discretization and initialization.
2. **Ahluwalia et al., Physical Review B 78, 054110 (2008)** — primary multiscale-method and next complete solver anchor. The full paper explicitly maps MD thermodynamics, domain walls, fluctuations, relaxation time and cell volume into LGD/TDGL quantities.
3. **Su et al., Nature Communications 13, 4867 (2022)** — beta-PVDF homogeneous benchmark. Eq. (4), Eq. (5) and Supplementary Table 3 are now executable; the full domain map is source-gated because the provided source does not tabulate all gradient/kinetic inputs.
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
docs/model_notes/GUO2024_REPRODUCIBILITY_GAPS.md
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

### Su 2022

```text
src/pvdf_pf/literature/su2022.py
scripts/reproduce_su2022_pvdf_landau.py
tests/test_su2022_pvdf.py
docs/model_notes/SU2022_PVDF_REPRODUCTION.md
```

The module directly encodes the published uniaxial PVDF bulk free energy and Supplementary Table-3 coefficients. `eps_b` remains an explicit argument in the electrostatic helper because the source varies it with MXene content. The report's 300/315/330 K minima are project algebra checks at declared temperatures, not source-reported simulation outputs.

## Next numerical stage

The next **source-complete solver validation** is Ahluwalia 2008: first reproduce the homogeneous `P(T)` / intrinsic `P-E` thermodynamics from the printed LGD coefficients, then the stochastic vector TDGL calibration and switching protocol. This gives us a complete numerical backbone before inserting crystal/OAF/IAF physics.

Guo 2024 remains the topological target, but exact full spiral reproduction is not forced by guessing its unreported dielectric/kinetic/mesh/initialization inputs. Any interim Guo PDE run with project-chosen values must be labeled sensitivity analysis rather than source reproduction.

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
