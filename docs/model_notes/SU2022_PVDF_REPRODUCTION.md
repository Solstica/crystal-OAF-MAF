# Su 2022 beta-PVDF reproduction log

Primary source: Yuanjie Su et al., *Nature Communications* **13**, 4867 (2022), DOI `10.1038/s41467-022-32518-3`.

## Role in this project

Su 2022 is not adopted as a model of OAF/IAF semicrystalline BOPVDF. It is used as a clean published **uniaxial PVDF phase-field algebra benchmark** and as an independent reference for the conventional sixth-order Landau polynomial convention used when expanding Guo 2024's contracted coefficients.

## Directly printed PVDF bulk free energy

The main article explicitly gives

```text
f_bulk = alpha1 Px^2 + alpha2 Py^2 + alpha3 Pz^2
       + alpha33 Pz^4 + alpha333 Pz^6.
```

The authors explain that the expression assumes uniaxial anisotropy with spontaneous polarization along `z`, the electrospinning poling direction, while `Px` and `Py` are treated as paraelectric and their higher-order terms are neglected.

Supplementary Table 3 gives

| coefficient | value |
|---|---:|
| `alpha1` | 5.647e9 |
| `alpha2` | 5.647e9 |
| `alpha3(T)` | `1.412e7 (T-315)` |
| `alpha33` | -1.842e11 |
| `alpha333` | 2.585e13 |
| `Q11` | -8.5 m^4 C^-2 |
| `Q12` | 0 |
| `Q44` | 0 |
| `s11` | 4.0e-10 m^2 N^-1 |
| `s12` | 1.11e-9 m^2 N^-1 |
| `s44` | 1.25e-9 m^2 N^-1 |

The repository transcribes these directly and does not replace them with values from Guo 2024 or Ahluwalia 2008.

## Electrostatic equation

The article explicitly uses

```text
div(eps0 * eps_b * E + P) = 0.
```

For this work `eps_b` is varied to account for the MXene-dependent dielectric background. The repository therefore keeps `eps_b` as an explicit required input; there is no hidden default dielectric constant.

## Published numerical geometry

The article states that the domain calculation uses a three-dimensional `512 nm x 512 nm x 512 nm` system, discretized by `128 x 128 x 128` grid points with `4 nm` spacing, periodic boundary conditions for polarization/electric/mechanical-displacement fields, and an applied electrospinning field of `1.2e5 V m^-1`. The calculation is reported to use the MuPRO Ferroelectric module.

These values are preserved as source metadata but are not sufficient by themselves to reproduce the domain maps.

## Algebraic checkpoint

At a deliberately declared **project test temperature** of `300 K`,

```text
alpha3 = -2.118e8 J m C^-2.
```

Solving the source polynomial analytically along `z` gives a project-derived nonnegative minimum

```text
Pz = 0.07258678047 C m^-2
f  = -2.4484662e6 J m^-3.
```

This number is **not reported by Su et al.** It exists only as a deterministic regression test of the source equation and Table-3 transcription. The reproduction report also evaluates 315 K and 330 K for algebraic continuity checks, with the same provenance label.

## Full-domain reproducibility limit

The article writes the total free energy as bulk + elastic + electric + gradient and gives the TDGL evolution equation. However, in the provided article/SI, a numerical PVDF gradient coefficient and kinetic coefficient are not tabulated. The available source therefore supports an exact Eq. (4) / Table-3 homogeneous benchmark and the form of Eq. (5), but not an exact independent recreation of the full MuPRO domain morphology.

This is recorded as

```text
SU2022_HOMOGENEOUS_PVDF_ALGEBRA = SOURCE_COMPLETE
SU2022_FULL_DOMAIN_MAP_REPRODUCTION = SOURCE_NUMERICALLY_INCOMPLETE
```

No gradient coefficient will be imported from an inorganic material or another polymer paper merely to make the full-domain simulation run.

## Implementation

```text
src/pvdf_pf/literature/su2022.py
tests/test_su2022_pvdf.py
scripts/reproduce_su2022_pvdf_landau.py
```

## Relation to the MLP bridge

Su 2022 is useful because it clarifies what a minimal uniaxial polymer free-energy model looks like, but it does not itself demonstrate atomistic calibration of these coefficients. Ahluwalia 2008 remains the stronger methodological precedent for converting atomistic observables into continuum parameters.
