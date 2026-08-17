# Guo 2024 reproduction log

Primary source: Mengfan Guo et al., *Nature Communications* **15**, 348 (2024), DOI `10.1038/s41467-023-44395-5`.

## Published phase-field problem

The paper evolves vector polarization `P=(Px,Py,Pz)` with TDGL,

```text
dPi/dt = -L delta F/delta Pi,
```

where `F` contains Landau, gradient, elastic/electrostrictive and electric energies. The simulated P(VDF-TrFE) nanodisk is 10 nm thick, radius 345 nm, at 25 C under short-circuit electrical boundary conditions and is solved by finite elements.

Supplementary Fig. S16 provides the first homogeneous-to-texture benchmark: strong-anisotropy Landau surface and polarization map in panels a,b; weak-anisotropy surface and spiral-like polarization map in panels c,d. Supplementary Tables S2/S3 contain the corresponding coefficient sets.

---

## Stage 1 — source-faithful one-axis Landau benchmark

Before reconstructing the vector polynomial, the branch first verified the part that is independent of any cross-term convention. Along one polarization axis, with the other two components exactly zero,

```text
f(P,T) = alpha1(T) P^2 + alpha11 P^4 + alpha111 P^6.
```

All cross terms vanish identically.

### Direct source inputs

Both S2 and S3 use

```text
alpha1(T) = 1.412e5 (T_C - 42) J m C^-2.
```

At 25 C, `alpha1=-2.4004e6 J m C^-2`.

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

These are `DIRECT_TRANSCRIPTION_FROM_GUO2024_SUPPLEMENTARY_TABLE`.

### Derived one-axis checkpoint

At 25 C the stable one-axis minima are

| parameter set | `|P|min` (C m^-2) | `fmin` (J m^-3) |
|---|---:|---:|
| S2 strong | 0.0730143463 | -25981.4041 |
| S3 weak | 0.0240959089 | -949.833981 |

These are **not source-quoted values**. They are `PROJECT_REPRODUCTION_DERIVED_FROM_DIRECT_SOURCE_COEFFICIENTS` regression anchors.

---

## Stage 2a — homogeneous vector Landau reconstruction for Fig. S16

The newly available primary SI closes the parameter-table gap but confirms an important detail: Guo et al. print the generic tensor form of the sixth-order Landau energy and the contracted coefficients, but do **not** print the fully expanded three-component polynomial.

To avoid inventing multiplicities, the expansion convention is cross-checked against Su et al., *Nature Communications* **13**, 4867 (2022), which explicitly prints the conventional sixth-order cubic Landau polynomial used in a closely related phase-field implementation:

```text
f = a1 sum(P_i^2)
  + a11 sum(P_i^4)
  + a12 sum_{i<j}(P_i^2 P_j^2)
  + a111 sum(P_i^6)
  + a112 sum_i[P_i^4 sum_{j!=i} P_j^2]
  + a123 Px^2 Py^2 Pz^2.
```

The repository therefore tags this expansion as

`PUBLISHED_CONVENTION_CROSSCHECK_AGAINST_SU2022_EQ3`,

not as a verbatim expanded Guo equation.

### Weak-anisotropy coefficients, direct from S3

```text
alpha11  = -1.842e8
alpha12  = -1.4736e9
alpha111 =  2.585e12
alpha112 =  9.6e12
alpha123 =  1.0857e13
```

with units as printed in S3. Elastic, electrostrictive and gradient parameters remain separate and are not involved in this homogeneous Stage 2a calculation.

### Strong-anisotropy control

S2 lists only `alpha1`, `alpha11`, and `alpha111` among the Landau coefficients. The homogeneous control reconstruction therefore includes only the source-listed Landau terms for this parameter set. This is a reproduction convention tied to the S2 table, not a claim that all omitted couplings are universally zero for P(VDF-TrFE).

### High-symmetry directional minima at 25 C

For `f(r n)=A2 r^2+A4 r^4+A6 r^6`, the analytic minimum along fixed direction `n` gives the following project-derived checkpoints:

| set | direction | `|P|min` (C m^-2) | `fmin` (J m^-3) |
|---|---|---:|---:|
| strong | [100] | 0.07301435 | -25981.40 |
| strong | [110] | 0.10325788 | -51962.81 |
| strong | [111] | 0.12646456 | -77944.21 |
| weak | [100] | 0.02409591 | -949.83 |
| weak | [110] | 0.02377743 | -953.80 |
| weak | [111] | 0.02452872 | -1029.49 |

A 4096-direction sphere sample gives, approximately,

```text
strong: radius max/min ~1.73, well-depth max/min ~3.00
weak:   radius max/min ~1.03, well-depth max/min ~1.08
```

Thus the reconstructed S2 and S3 sets produce the intended **large contrast in angular anisotropy**: S2 is strongly direction-dependent, while S3 is nearly isotropic. This is consistent with the visibly lobed strong-anisotropy surface and nearly rounded weak-anisotropy surface in Supplementary Fig. S16.

### What is and is not closed

Stage 2a closes the **homogeneous anisotropy checkpoint**. It does not claim pixel-exact reproduction of Fig. S16 because the source does not state the exact graphical radial surface definition, energy normalization or color normalization used to render that figure.

Therefore the valid statement is:

`SOURCE_CONSTRAINED_QUALITATIVE_S16_ANISOTROPY_REPRODUCTION`,

not `EXACT_FIGURE_REPRODUCTION`.

Implementation:

```text
src/pvdf_pf/literature/guo2024_vector.py
scripts/reproduce_guo2024_s16_landau.py
tests/test_guo2024_vector_landau.py
```

---

## Mechanical-model limitation from the peer-review record

The peer-review file is unusually useful for model governance. Reviewer #3 challenged the use of linear elasticity / cubic parent elastic constants for a ferroelectric polymer, especially because elastic energy contributes strongly in the reported simulations.

The authors' response does **not** claim a general nonlinear polymer-mechanics validation. They state that non-cubic symmetry is represented through electrostrictive eigenstrain and that their simulations are restricted to the linear elastic, low-field regime; simulated stresses do not exceed about 12 MPa. They explicitly characterize the phase-field calculation as effective for **semi-qualitatively** addressing the low-field rotational mechanism before destruction.

Consequently this repository will not use Guo's elastic model as evidence for quantitatively calibrated nonlinear polymer mechanics. Reproducing it is still valuable, but its claim level must remain the same as the source authors' own limitation.

---

## Next Guo checkpoint

The next Guo-specific numerical stage is the published short-circuit **vector TDGL** problem. The implementation order is deliberately incremental:

```text
homogeneous vector Landau
-> vector TDGL + published gradient coefficients + electrostatics
-> elastic/electrostrictive coupling under the source's low-field limitation
-> qualitative weak-anisotropy spiral texture (Fig. 3e / S16d)
-> field/stress evolution (S17-S19)
```

No Huang/Rui OAF/IAF substitution, project-defined morphology spectrum, frozen source allocation or MLP-derived coefficient is allowed inside this Guo reproduction benchmark.
