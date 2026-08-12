# v0.1.3 — real BDS ingestion and OAF subpartition

## Why this stage exists

v0.1.2 established a numerically correct frequency-to-time Debye bridge, but it used synthetic spectra for verification. v0.1.3 begins the source-driven stage.

The primary source is Rui et al., *Macromolecules* 2022, DOI `10.1021/acs.macromol.2c01110`. Its Supporting Information explicitly separates the oriented amorphous fraction into rigid OAF (ROAF) and mobile OAF (MOAF) and writes the film dielectric response using crystal, ROAF, MOAF and IAF contributions.

The project therefore refines its structural hierarchy as

```text
crystal / OAF / MAF
             |
             +-- OAF = ROAF + MOAF   [source-supported decomposition]

project MAF <-> source IAF            [provisional mapping; not yet identity]
```

## Source-model equation

The SI mixture identity is represented as

```text
epsilon_film(T)
 = epsilon_cr * eta_cr
 + epsilon_ROAF * eta_ROAF(T)
 + epsilon_MOAF(T) * eta_MOAF(T)
 + epsilon_IAF(T) * eta_IAF.
```

The SI sets `epsilon_cr = 3.0` and assumes `epsilon_ROAF = epsilon_cr`. These remain `SOURCE_MODEL_ASSUMPTION` values in this repository.

`epsilon_MOAF(T)` is a source-model-derived quantity. It is not a direct local dielectric measurement.

## Numerical data required from the SI

- Figure S1: `n(T)`, `m_d(T)`, `g(T)` for unpoled and poled BOPVDF.
- Figure S2: `epsilon_MOAF(T)` for unpoled and poled BOPVDF.
- Figure S3 / corresponding SI mobility plot: `lambda(T)` and `mu_r(T)`.

The public ACS index exposes the SI descriptions and model equation, but not machine-readable numerical curve data. The repository therefore contains digitization templates and refuses to promote any plot-derived value lacking explicit provenance.

## Dynamic model after digitization

The first time-domain representation will keep two relaxational channels distinct:

```text
MOAF channel: tau_MOAF(T), Delta-epsilon_MOAF(T)
IAF  channel: tau_IAF(T),  Delta-epsilon_IAF(T)
```

with

```text
tau_j dP_j/dt + P_j = epsilon0 Delta-epsilon_j E_local.
```

Total displacement will eventually be assembled as

```text
D = epsilon0 epsilon_b E + P_FE + P_MOAF + P_IAF,
```

where `P_FE` remains the nonlinear crystal/OAF ferroelectric variable and the two auxiliary terms represent dielectric relaxation.

## Guardrails

1. No temperature-dependent MOAF or IAF parameter is hard-coded from visual estimates.
2. `project MAF = source IAF` is not yet promoted from a provisional mapping.
3. A measured low-frequency film permittivity is never copied directly into TDGL `epsilon_b`.
4. BDS relaxation times in seconds remain outside the dimensionless TDGL integrator until a physical time mapping is calibrated.
5. The simple Debye channel is a first dynamical representation. If digitized spectra show broad/non-Debye peaks, the next model must use a relaxation-time distribution or HN/Prony representation rather than forcing a single Debye fit.

## v0.1.3 completion condition

This stage is ready for physical fitting only when the required Figure S1/S2/S3 curves have populated source CSVs with `DIGITIZED_SOURCE` or `DIRECT_TABULATED` provenance and units have been checked against the paper/SI.
