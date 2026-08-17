# Guo 2024 reproduction log

Primary source: Mengfan Guo et al., *Nature Communications* **15**, 348 (2024), DOI `10.1038/s41467-023-44395-5`.

## Stage 1 — source-faithful one-axis Landau benchmark

The published phase-field model evolves the vector polarization `P=(Px,Py,Pz)` by TDGL and contains Landau, gradient, elastic/electrostrictive and electric energies. Supplementary Tables S2 and S3 report strong- and weak-anisotropy coefficient sets.

Before reproducing the full vector polar-spiral calculation, this branch checks the part of the Landau model that is unambiguous without any reconstruction of cross-term multiplicities. Along a single polarization axis, with the other two components exactly zero, the sixth-order polynomial reduces to

```text
f(P,T) = alpha1(T) P^2 + alpha11 P^4 + alpha111 P^6.
```

All `alpha12`, `alpha112` and `alpha123` cross terms vanish identically. No value or tensor-expansion convention for those terms is assumed in Stage 1.

### Direct source inputs

The temperature is the paper's phase-field simulation temperature, 25 C.

Both S2 and S3 use

```text
alpha1(T) = 1.412e5 (T_C - 42) J m C^-2.
```

At 25 C this gives `alpha1 = -2.4004e6 J m C^-2`.

Strong-anisotropy S2 axis coefficients:

```text
alpha11  = -1.842e9  J m^5 C^-4
alpha111 =  2.585e11 J m^9 C^-6
```

Weak-anisotropy S3 axis coefficients:

```text
alpha11  = -1.842e8  J m^5 C^-4
alpha111 =  2.585e12 J m^9 C^-6
```

These values are `DIRECT_TRANSCRIPTION_FROM_GUO2024_SUPPLEMENTARY_TABLE`.

### Derived checkpoint

For `P != 0`, stationarity of the one-axis polynomial gives

```text
3 alpha111 y^2 + 2 alpha11 y + alpha1 = 0,
y = P^2.
```

Using only the source coefficients, the 25 C stable one-axis minima are:

| parameter set | `|P|min` (C m^-2) | `fmin` (J m^-3) | barrier `f(0)-fmin` (J m^-3) |
|---|---:|---:|---:|
| S2 strong anisotropy | 0.0730143463 | -25981.4041 | 25981.4041 |
| S3 weak anisotropy | 0.0240959089 | -949.833981 | 949.833981 |

These are **not values quoted by Guo et al.** They are labelled `PROJECT_REPRODUCTION_DERIVED_FROM_DIRECT_SOURCE_COEFFICIENTS` and serve only as regression anchors for the exact polynomial transcription.

## What this does and does not reproduce

Stage 1 verifies that the repository can reproduce an unambiguous mathematical slice of the published free energy without introducing a project hypothesis. It does **not** yet reproduce Fig. 3e, a spiral texture, a switching curve, a domain period, or a time scale.

The following are intentionally absent at this stage:

- expanded three-component Landau polynomial;
- gradient term implementation for this source;
- elastic equilibrium and electrostrictive eigenstrain;
- short-circuit electrostatic solve;
- kinetic coefficient/time calibration;
- finite-element nanodisk geometry;
- external electric or mechanical loading.

## Gate before Stage 2

Stage 2 may begin only after the exact expanded Landau convention associated with the reported `alpha11`, `alpha12`, `alpha111`, `alpha112`, `alpha123` coefficients is checked against the primary paper's cited phase-field formulation. The issue is not whether a conventional cubic sixth-order expression is familiar; the issue is whether its multiplicity convention is exactly the one used to generate Guo et al.'s reported phase diagram and spiral state.

Once that is verified, the next benchmark is the complete homogeneous Landau surface, followed by the published short-circuit vector TDGL problem. Until then, no guessed cross-term factors are admitted.
