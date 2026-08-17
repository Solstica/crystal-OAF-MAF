# Literature reproduction baseline

This branch resets the main development logic from project-defined sensitivity scans to source-first reproduction.

## Scope reset

The v0.1.16-v0.1.19 morphology/frozen-source studies remain in repository history as numerical sensitivity experiments. They are **not** treated as literature-reproduced phase-field physics and are not extended further until a published ferroelectric-polymer phase-field baseline has been reproduced.

The next development rule is:

1. reproduce a published PVDF/P(VDF-TrFE) phase-field model using the paper's own equations, coefficients, boundary conditions and target figures;
2. document every quantity that is not explicitly reported;
3. only after the baseline closes, replace or augment a published coefficient with Huang/Rui crystal-OAF-IAF evidence or atomistic/MLP output.

## Candidate 1 — Ahluwalia et al., Phys. Rev. B 78, 054110 (2008)

**Title:** *Multiscale kinetic model for polarization switching in ferroelectric polymer thin films*  
**DOI:** 10.1103/PhysRevB.78.054110

### Why it matters

This is the closest published precedent to the long-term project goal. The authors parameterize a continuum Landau-Ginzburg-Devonshire model for ideal all-trans P(VDF-TrFE) 70:30 using molecular-dynamics data. They also set length scale, time scale and thermal-noise amplitude from MD and then solve a time-dependent Ginzburg-Landau switching model.

Therefore this paper is the **multiscale-method anchor** for the future DFT/MD/MLP -> phase-field bridge. It directly establishes that atomistic data can be used to parameterize a ferroelectric-polymer continuum switching model rather than merely serving as qualitative reference.

### Current reproducibility status

The publisher abstract and bibliographic record are verified, but the complete article text/parameter table is not present in the current local source set. We therefore do **not** transcribe coefficients from secondary citations or reconstruct them from memory.

Status: `PRIMARY_TARGET_WAITING_FULL_TEXT`.

Required before exact reproduction: full article PDF (or an author manuscript containing equations, MD-to-LGD fitting relations, scales, noise strength and numerical protocol).

## Candidate 2 — Guo et al., Nature Communications 15, 348 (2024)

**Title:** *Electrically and mechanically driven rotation of polar spirals in a relaxor ferroelectric polymer*  
**DOI:** 10.1038/s41467-023-44395-5

### Published phase-field model

The order parameter is the polarization vector

`P = (Px, Py, Pz)`.

The temporal evolution is the published TDGL equation

`dPi/dt = -L delta F/delta Pi`.

The total free energy contains the four standard contributions printed in the article:

`F = integral_V (f_Land + f_grad + f_elas + f_ele) dV`.

The article explicitly gives

- generic sixth-order Landau free energy;
- gradient energy `0.5 G_ijkl P_i,j P_k,l`;
- elastic energy `0.5 C_ijkl (eps_ij-eps0_ij)(eps_kl-eps0_kl)`;
- electrostrictive eigenstrain `eps0_ij = Q_ijkl P_k P_l`;
- electric energy `-Ei(Pi + 0.5 eps0 kappa_ij Ej)`.

The simulated P(VDF-TrFE) nanodisk is 10 nm thick with radius 345 nm, short-circuit electrical boundary condition, at 25 C. The authors used finite elements.

### Direct parameter set available in the Supplementary Information

Supplementary Table S3 gives the weak-anisotropy coefficients used to obtain the polar-spiral state:

| coefficient | published value |
|---|---:|
| `alpha1` | `1.412 (T-42) x 10^5 J m C^-2` |
| `alpha11` | `-1.842 x 10^8 J m^5 C^-4` |
| `alpha12` | `-1.4736 x 10^9 J m^5 C^-4` |
| `alpha111` | `2.585 x 10^12 J m^9 C^-6` |
| `alpha112` | `9.6 x 10^12 J m^9 C^-6` |
| `alpha123` | `1.0857 x 10^13 J m^9 C^-6` |
| `c11` | `4.88 x 10^10 J m^-3` |
| `c12` | `5.6 x 10^9 J m^-3` |
| `c44` | `2.16 x 10^10 J m^-3` |
| `Q11` | `-0.0162 m^4 C^-2` |
| `Q12` | `0.0441 m^4 C^-2` |
| `Q44` | `-0.12 m^4 C^-2` |
| `G11` | `9.96 x 10^-10 J m^3 C^-2` |
| `G12` | `0` |
| `G44/G44'` | `4.98 x 10^-10 J m^3 C^-2` |

Supplementary Table S2 separately reports the strong-anisotropy control set. We will not infer omitted cross-coefficients from that table unless the source explicitly defines them.

### Reproduction target

This paper is selected as the **first executable benchmark** because its open article and SI provide enough information to begin a source-faithful implementation now.

Reproduction order:

1. one-dimensional Landau axis slices using only directly listed coefficients;
2. weak-anisotropy Landau surface, after checking the coefficient-index expansion convention against the source/reference chain;
3. short-circuit vector TDGL without external stress;
4. reproduce the qualitative low-anisotropy polar texture shown in main Fig. 3e / Supplementary Fig. S16d;
5. only then add elastic/electrostrictive and field-driven evolution needed for Supplementary Figs. S17-S19.

Status: `EXECUTABLE_PRIMARY_BENCHMARK`.

## Candidate 3 — Su et al., Nature Communications 13, 4867 (2022)

**Title:** *High-performance piezoelectric composites via beta phase programming*  
**DOI:** 10.1038/s41467-022-32518-3

### Published PVDF phase-field model

The article explicitly uses vector TDGL and the same four free-energy classes. For PVDF it prints the uniaxial bulk energy

`f_bulk = alpha1 Px^2 + alpha2 Py^2 + alpha3 Pz^2 + alpha33 Pz^4 + alpha333 Pz^6`,

with spontaneous polarization along the electrospinning/poling `z` direction. The electrostatic equilibrium equation is

`div(eps0 eps_b E + P) = 0`.

The reported simulation cell is 512 x 512 x 512 nm^3, discretized on 128^3 nodes (4 nm spacing), with periodic boundary conditions for polarization, electric field and mechanical displacement, under `1.2 x 10^5 V/m` applied field. The paper states that PVDF constants are in Supplementary Table 3.

### Why it matters

This is the cleanest directly printed **uniaxial beta-PVDF Landau form** among the selected papers and is structurally closer to the scalar-polarization simplifications previously used in this repository. It is therefore the second executable transfer benchmark after Guo 2024.

Status: `SECONDARY_EXECUTABLE_BENCHMARK`; full SI parameter transcription still needs to be completed from the primary supplementary file.

## Selection decision

The project will **not** continue from v0.1.19 into v0.1.20 morphology inventions.

The literature-reproduction mainline is now:

`Guo 2024 exact/source-faithful benchmark -> Su 2022 beta-PVDF benchmark -> Ahluwalia 2008 MD-to-LGD bridge once full text is available -> Huang/Rui crystal-OAF-IAF substitution -> project DFT/MD/MLP parameter generation`.

The reason Guo 2024 is executed first rather than Ahluwalia 2008 is practical reproducibility, not scientific priority: Guo's complete open SI exposes the coefficients and target structures, while the full Ahluwalia parameterization is not yet in the current source set.

## Evidence boundary

Huang 2021 and Rui 2021/2022 remain the experimental source for crystal/OAF/IAF structure and amorphous dynamics. They do not by themselves provide a complete crystal/OAF/IAF TDGL free-energy parameterization. Their data must therefore enter only after a published polymer phase-field baseline has been reproduced and the exact continuum coefficient being replaced is identified.

No `PROJECT_IMAGE_DERIVED_CONSTRAINT`, sinusoidal waviness spectrum, frozen beta/OAF source allocation, or project-defined OAF decay profile is permitted in the reproduction benchmark.
