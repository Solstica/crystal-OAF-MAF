# Ahluwalia 2008 multiscale reproduction log

Primary source: Rajeev Ahluwalia et al., *Physical Review B* **78**, 054110 (2008), DOI `10.1103/PhysRevB.78.054110`.

## Why this paper is the multiscale anchor

This paper does not use atomistic simulation merely as qualitative support. It explicitly transfers molecular-dynamics (MD) information into a continuum Landau-Ginzburg-Devonshire / time-dependent Ginzburg-Landau model for ideal all-trans P(VDF-TrFE) 70:30. The transfer chain is

```text
MD thermodynamics -> homogeneous LGD coefficients
MD 180-degree walls -> gradient coefficients
MD polarization fluctuations -> TDGL thermal-noise amplitudes
MD relaxation time -> TDGL physical time scale
MD cell volume -> continuum grid length scale
```

That architecture is the closest published precedent to the project's eventual DFT/MD/MLP -> phase-field bridge.

## 1. Homogeneous thermodynamics: MD -> LGD

The source uses

```text
f(Pz,T) = alpha0 (T-T0) Pz^2 / 2 - beta Pz^4 / 4 + gamma Pz^6 / 6.
```

The printed MD observables in Table I are

| quantity | source value |
|---|---:|
| `P0` | 0.111 C m^-2 |
| `Tc` | 450 K |
| `Pc` | 0.088 C m^-2 |
| `f0` | -2.12e8 J m^-3 |

The source then gives Eqs. (2)-(3) for deriving the four LGD parameters and independently prints the resulting Table II:

| coefficient | source value |
|---|---:|
| `alpha0` | 9.02e7 J m C^-2 K^-1 |
| `T0` | 252 K |
| `beta` | 9.20e12 J m^5 C^-4 |
| `gamma` | 8.88e14 J m^9 C^-6 |

### Precision audit

Applying the printed equations to the **rounded values printed in Table I** gives

| coefficient | re-derived from rounded Table I | difference from printed Table II |
|---|---:|---:|
| `alpha0` | 9.0676924e7 | +0.529% |
| `T0` | 248.2113 K | -1.503% |
| `beta` | 9.4512308e12 | +2.731% |
| `gamma` | 9.1534389e14 | +3.079% |

This is a useful reproducibility result: Table I is printed with insufficient precision to regenerate Table II exactly. Therefore **the directly printed Table II is authoritative for source reproduction**. The re-derived values are retained only as `PROJECT_REDERIVED_FROM_ROUNDED_SOURCE_VALUES` and must not overwrite Table II.

### Homogeneous P(T) checkpoint from printed Table II

The positive nonzero stationary branch of the printed Eq. (1) has

```text
P^2 = [beta + sqrt(beta^2 - 4 gamma alpha0(T-T0))] / (2 gamma).
```

Using **the directly printed Table-II coefficients**, not the re-derived rounded-Table-I fit, gives

| T | project-derived stationary `P` |
|---:|---:|
| 0 K | 0.11145019 C m^-2 |
| 300 K | 0.09932864 C m^-2 |
| 450 K | 0.08816395 C m^-2 |

The 0 K and 450 K values differ from the printed Table-I `P0=0.111` and `Pc=0.088` only by the expected rounding of the published tables.

For the first-order polynomial, exact coexistence implied by rounded Table II is

```text
Tc = 450.1332774 K
Pc = 0.08814914 C m^-2,
```

rather than exactly 450 K and 0.088 C m^-2. This ~0.13 K mismatch is another independent indication that the printed coefficients are rounded.

### Homogeneous intrinsic P-E checkpoint

For homogeneous field-controlled switching,

```text
E(P,T) = df/dP
       = alpha0(T-T0) P - beta P^3 + gamma P^5.
```

Loss of stability of the positive ferroelectric branch satisfies `dE/dP=0`. At a declared 300 K checkpoint, the printed Table-II parameters give

```text
P_spinodal = 0.07781500 C m^-2
E_switch(+ -> -) = -1.46442874e9 V m^-1.
```

The magnitude, `1.464 GV m^-1`, is consistent with the field scale of the source Fig. 3. It is a deterministic calculation from the source equation and coefficients, not a separately quoted paper number.

The homogeneous thermodynamic layer is therefore now closed sufficiently to use as the first solver regression target before stochastic spatial TDGL.

## 2. Domain-wall width: MD -> gradient coefficients

The source constructs two 180-degree domain walls and estimates

```text
xi1 = xi2 ~= 0.4 nm.
```

Using the analytical domain-wall relation in its Eq. (10), the paper reports

```text
K1 = K2 = 2.108e-8 J m^3 C^-2.
```

This is a direct atomistic-to-gradient calibration path and is particularly important for our later MLP stage: instead of guessing a continuum gradient penalty, an atomistic model can generate a domain wall, its profile/width or excess energy can be measured, and the continuum gradient coefficient can be fitted to that observable.

The source explicitly states that an analogous head-to-head/tail-to-tail wall for determining `K3` is electrostatically unstable. It therefore **arbitrarily sets `K3=K1=K2` for computational convenience**. The repository records this as a source-model assumption, not as an MD-calibrated coefficient.

## 3. MD fluctuations -> TDGL thermal noise and kinetic anisotropy

For a 300 K sudden-heating comparison, the source matches TDGL polarization fluctuations to MD and reports dimensionless noise amplitudes

```text
noise_x = 0.0711
noise_y = 1.0669
noise_z = 0.3556.
```

Its Eq. (11) defines

```text
noise_x = sqrt(m) noise_z,   Gamma_x = m Gamma_z
noise_y = sqrt(n) noise_z,   Gamma_y = n Gamma_z.
```

The printed amplitudes therefore imply the project-derived checks

```text
m ~= 0.03998  (~0.04)
n ~= 9.00169 (~9.0)
```

without introducing additional fitting.

## 4. MD relaxation -> physical TDGL time

The source observes the TDGL model equilibrating near `t*=5`, while the corresponding MD equilibration time is about 45 ps. This maps

```text
1 t* ~= 9 ps,
```

and the paper explicitly adopts about 9 ps as the smallest TDGL time scale.

The source additionally uses a `138.24 nm x 138.24 nm x 138.24 nm` TDGL cell with a smallest length of `2.16 nm`, chosen so that the continuum cell volume matches the MD simulation-cell volume. This is another explicit atomistic-to-continuum scale connection.

## 5. Electrostatics and switching boundary conditions

The source evolves all three polarization components with stochastic TDGL and enforces electrostatics through

```text
-eps0 laplacian(phi) + div(P) = 0.
```

For the switching film, the source uses periodic polarization boundary conditions in-plane, `dPi/dz=0` at film-normal boundaries, assumes perfect compensation of surface bound charge, and applies short-circuit potential boundaries before the time-dependent voltage drive.

These conditions belong to the **source benchmark only**. They are not automatically transferable to the later semicrystalline crystal/OAF/IAF model.

## 6. What this changes for the project

The full article closes the previous evidence gap. Ahluwalia 2008 is now `EXECUTABLE_MULTISCALE_ANCHOR`, not `PRIMARY_TARGET_WAITING_FULL_TEXT`.

For the future MLP bridge, the useful template is not “ML predicts a phase-field parameter directly” unless that parameter has a uniquely defined atomistic observable. The safer chain is

```text
DFT/MLP ensemble
  -> P(T), phase-energy differences, wall profiles/energies, fluctuation statistics,
     relaxation times, elastic/electrostatic observables
  -> fit the corresponding LGD / gradient / kinetic / noise coefficients
  -> phase-field model
  -> continuum observables
```

This is the central methodological lesson adopted from the paper.

## Current implementation

```text
src/pvdf_pf/literature/ahluwalia2008.py
tests/test_ahluwalia2008_multiscale.py
scripts/reproduce_ahluwalia2008_multiscale.py
```

The JSON report preserves direct printed quantities separately from values re-derived from rounded tables, records the homogeneous P(T)/intrinsic-spinodal checkpoints, and keeps `K3` explicitly marked as a source computational assumption.

## Remaining reproduction before using it as solver validation

The next stage is the spatial stochastic TDGL rather than more homogeneous fitting. It must implement the source's transverse susceptibility terms, dimensionless rescaling, gradient coefficients, electrostatic constraint, noise calibration and film switching boundary conditions before moving to Fig. 6/Fig. 7.

A successful switching reproduction would validate the numerical backbone. It still would not make the ideal all-trans 70:30 model a direct physical model of semicrystalline BOPVDF; crystal/OAF/IAF physics remains a subsequent substitution/calibration step.
