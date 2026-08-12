# Rui et al. 2022 BDS source-data staging

Primary source: G. Rui, E. Allahyarov, J. J. Thomas, P. L. Taylor, L. Zhu, *Macromolecules* **2022**, 55, 9705-9714. DOI: `10.1021/acs.macromol.2c01110`.

Supporting-information file: `ma2c01110_si_001.pdf`.

## Source-supported structure

The SI reports:

- Figure S1: `n(T)`, `m_d(T)`, and `g(T)` for unpoled and poled BOPVDF.
- Figure S2: dielectric constant of the mobile OAF, `epsilon_MOAF(T)`.
- a source mixture model in which OAF is split into rigid OAF (ROAF) and mobile OAF (MOAF), alongside crystal and IAF.
- Figure S3: `lambda(T)` and rotational mobility `mu_r(T)`.

Section S2-S3 additionally states the source-model partition/assumptions `eta_cr=0.6`, constant `eta_IAF=0.2`, and `epsilon_cr=epsilon_ROAF=3.0`.

## Digitization status

The user-supplied SI was rendered at 600 dpi and Figures S1-S3 were digitized from the plotted markers.

Populated files:

- `figure_s1_digitized.csv`: 48 points.
- `figure_s2_digitized.csv`: 30 points.
- `figure_s3_digitized.csv`: 48 points.

All plot-derived values use:

```text
provenance = DIGITIZED_SOURCE
```

No guide-curve interpolation was used during transcription. See `docs/model_notes/RUI2022_SI_DIGITIZATION.md` for axis calibration, uncertainty and interpretation limits.

## Interpretation boundary

The SI explicitly notes that the Kirkwood-Frohlich treatment cannot separately solve `n(T)`, `m_d(T)`, and `g(T)` for the oriented MOAF. `epsilon_MOAF(T)` is source-model-derived, while the reported `n`, `m_d`, and `g` belong to the broader source amorphous treatment. The project therefore does not assign those quantities directly to MOAF without an additional validated mapping.

Interpolation in temperature belongs to a later calibration step and must remain distinguishable from the digitized source points.
