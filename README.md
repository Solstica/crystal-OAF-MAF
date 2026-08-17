# crystal-OAF-MAF

Phase-field development repository for semicrystalline PVDF/BOPVDF with an eventual **crystal / OAF / IAF** description and a DFT/MD/MLP -> continuum parameter bridge.

> Historical scripts may still use `MAF` as a phase identifier. In the Huang/Rui structural description the physical non-oriented amorphous region is treated here as IAF; RAF/MAF remain mobility-defined bookkeeping quantities and are not synonyms for OAF/IAF.

## Active branch: literature-reproduction

The current mainline is no longer “add another project-defined morphology scan”. It is **source-first reproduction of published ferroelectric-polymer phase-field models**.

The previous v0.1.16-v0.1.19 frozen-source / waviness / RMS-orientation / harmonic-spectrum studies remain in repository history as numerical sensitivity experiments. They are not considered reproduced PVDF phase-field physics and are not extended on this branch.

### Reproduction hierarchy

1. **Guo et al., Nature Communications 15, 348 (2024)** — first executable benchmark because the open paper/SI gives the vector TDGL framework and strong/weak-anisotropy coefficient tables.
2. **Su et al., Nature Communications 13, 4867 (2022)** — second benchmark because the paper prints an explicit uniaxial beta-PVDF Landau form and electrostatic equation.
3. **Ahluwalia et al., Phys. Rev. B 78, 054110 (2008)** — multiscale-method anchor because MD data are used to parameterize a P(VDF-TrFE) LGD/TDGL model; exact reproduction waits for the full article/author manuscript rather than reconstructing coefficients from secondary sources.
4. Only after a published polymer phase-field baseline closes do Huang/Rui crystal-OAF-IAF data replace identified continuum inputs.
5. DFT/MD/MLP calculations are then used to generate or validate the specific continuum coefficients that remain unresolved.

The detailed selection and evidence boundaries are in `docs/model_notes/LITERATURE_REPRODUCTION_BASELINE.md`.

## Current executable checkpoint

**Guo 2024 Stage 1: source-faithful one-axis Landau benchmark.**

The source reports a vector sixth-order Landau model. Before expanding the full vector polynomial, the repository implements only the unambiguous one-axis slice

```text
f(P,T) = alpha1(T) P^2 + alpha11 P^4 + alpha111 P^6
```

using the coefficients directly transcribed from Supplementary Tables S2/S3. Cross terms vanish on this axis, so no unreported multiplicity convention is introduced.

At 25 C, the code derives the stationary minima directly from those source coefficients and records them as reproduction-derived regression anchors. These derived values are not represented as measurements or as values quoted by Guo et al.

Implementation:

```text
src/pvdf_pf/literature/guo2024.py
scripts/reproduce_guo2024_landau_axis.py
tests/test_guo2024_landau_axis.py
docs/model_notes/GUO2024_REPRODUCTION.md
```

A separate `literature-reproduction` GitHub Actions workflow runs only this source-faithful benchmark and stores its JSON report as an artifact. The original exploratory branch CI is not used as evidence that this reproduction is physically correct.

## Gate for the next stage

The next code change is **not** another morphology parameter scan. Before implementing Guo's full vector Landau surface, the exact expansion convention associated with `alpha11`, `alpha12`, `alpha111`, `alpha112`, and `alpha123` must be verified from the primary source/reference chain or source data. Only then will the short-circuit vector TDGL problem and the polar-spiral target be implemented.

No `PROJECT_IMAGE_DERIVED_CONSTRAINT`, sinusoidal waviness spectrum, frozen beta/OAF polarization allocation, or project-defined OAF decay profile is allowed inside a literature-reproduction benchmark.

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
