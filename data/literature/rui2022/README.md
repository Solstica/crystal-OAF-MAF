# Rui et al. 2022 BDS source-data staging

Primary source: G. Rui, E. Allahyarov, J. J. Thomas, P. L. Taylor, L. Zhu, *Macromolecules* **2022**, 55, 9705-9714. DOI: `10.1021/acs.macromol.2c01110`.

Supporting-information file named by ACS: `ma2c01110_si_001.pdf`.

## Source-supported structure

The SI reports:

- Figure S1: `n(T)`, `m_d(T)`, and `g(T)` for unpoled and poled BOPVDF.
- Figure S2: dielectric constant of the mobile OAF, `epsilon_MOAF(T)`.
- a source mixture model in which OAF is split into rigid OAF (ROAF) and mobile OAF (MOAF), alongside crystal and IAF.
- further SI results for `lambda(T)` and rotational mobility `mu_r(T)`.

The article abstract states that active-dipole concentration increases from -30 to 40 degC and that rotational dipole mobility increases by more than four orders of magnitude over this window.

## Digitization policy

No numerical point is entered into the repository from visual inspection unless it is accompanied by its figure/panel, sample state, unit and provenance.

Use the CSV templates in this directory. For points read from plots, set:

```text
provenance = DIGITIZED_SOURCE
```

For values explicitly tabulated in the source, set:

```text
provenance = DIRECT_TABULATED
```

Do not interpolate missing temperatures during transcription. Interpolation belongs in a later calibration step and must remain distinguishable from source data.

## Current status

The accessible ACS index exposes the SI descriptions and equations but not machine-readable numerical curves. Therefore v0.1.3 implements the exact data schema and source-model algebra now, while leaving plot-derived numerical rows empty until the SI figures are digitized from the PDF at sufficient resolution.
