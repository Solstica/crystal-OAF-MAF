# Literature reproduction baseline

This branch resets the main development logic from project-defined sensitivity scans to source-first reproduction.

## Scope reset

The v0.1.16-v0.1.19 morphology/frozen-source studies remain in repository history as numerical sensitivity experiments. They are **not** treated as literature-reproduced phase-field physics and are not extended further until a published ferroelectric-polymer phase-field baseline has been reproduced.

The development rule is:

1. reproduce published PVDF/P(VDF-TrFE) continuum models using the papers' own equations, coefficients, boundary conditions and target observables;
2. distinguish directly printed quantities from values re-derived from rounded tables or reconstructed from a published convention;
3. preserve every source approximation as an approximation rather than silently upgrading it to a material constant;
4. only after a baseline closes, replace or augment an identified continuum input with Huang/Rui crystal-OAF-IAF evidence or atomistic/MLP output.

## Anchor A — Guo et al., Nature Communications 15, 348 (2024)

**Title:** *Electrically and mechanically driven rotation of polar spirals in a relaxor ferroelectric polymer*  
**DOI:** 10.1038/s41467-023-44395-5

### Published model

The order parameter is vector polarization `P=(Px,Py,Pz)`, evolved by TDGL. The total energy contains Landau, gradient, elastic/electrostrictive and electric terms. The source simulates a 10-nm-thick, 345-nm-radius P(VDF-TrFE) nanodisk at 25 C under short circuit using finite elements.

The SI gives strong- and weak-anisotropy coefficient sets and Supplementary Fig. S16 provides a direct energy-surface/polarization-texture target.

### Current reproduction status

- Stage 1: exact one-axis source polynomial — closed.
- Stage 2a: homogeneous three-component angular anisotropy — closed at a **source-constrained qualitative** level.
- The expanded sixth-order convention is not printed by Guo, so the repository cross-checks it against the explicit conventional polynomial printed by Su et al. 2022 rather than inventing multiplicities.
- Pixel-exact Fig. S16 surface rendering is not claimed because its graphical normalization is not reported.
- Vector TDGL/gradient/electrostatic and elastic spiral reproduction remain open.

Status: `EXECUTABLE_PRIMARY_PHASE_FIELD_BENCHMARK`.

### Important source limitation

The Guo peer-review record documents a reviewer challenge to the linear elastic/cubic-parent treatment. The authors defend it only for the low-field linear-elastic regime, with stresses no larger than about 12 MPa, and describe the phase-field result as semi-quantitative for the rotational mechanism before destruction. This scope limitation is inherited by our reproduction.

## Anchor B — Ahluwalia et al., Phys. Rev. B 78, 054110 (2008)

**Title:** *Multiscale kinetic model for polarization switching in ferroelectric polymer thin films*  
**DOI:** 10.1103/PhysRevB.78.054110

### Why it is now central

The full primary article is now available. It provides the clearest published atomistic-to-continuum transfer chain for a ferroelectric polymer:

```text
MD P(T), Pc, Tc, phase-energy difference -> homogeneous LGD coefficients
MD 180-degree domain walls -> K1, K2 gradient coefficients
MD fluctuations -> stochastic TDGL noise amplitudes / kinetic anisotropy
MD equilibration time -> physical TDGL time scale
MD cell volume -> continuum grid scale
```

The printed Table-I MD observables and Table-II LGD coefficients are now encoded separately. Re-applying the paper's Eqs. (2)-(3) to the rounded Table-I numbers gives percent-level differences from printed Table II, demonstrating that Table II must remain authoritative for exact source reproduction.

The paper directly reports `K1=K2=2.108e-8 J m^3 C^-2` from MD-estimated ~0.4 nm domain-wall widths. It explicitly states that `K3` could not be atomistically obtained from unstable head-to-head/tail-to-tail walls and was set equal to `K1=K2` for computational convenience. That distinction is preserved in code.

The source also reports matched stochastic-noise amplitudes and a mapping of approximately 9 ps per dimensionless TDGL time unit.

Status: `EXECUTABLE_MULTISCALE_ANCHOR`.

This paper is the principal methodological template for deciding **what an MLP must actually predict or reproduce before a phase-field coefficient can be calibrated**.

## Anchor C — Su et al., Nature Communications 13, 4867 (2022)

**Title:** *High-performance piezoelectric composites via beta phase programming*  
**DOI:** 10.1038/s41467-022-32518-3

### Published PVDF phase-field model

The article explicitly prints a sixth-order conventional vector Landau polynomial for the ceramic and a uniaxial PVDF bulk energy

```text
f_bulk = alpha1 Px^2 + alpha2 Py^2 + alpha3 Pz^2
       + alpha33 Pz^4 + alpha333 Pz^6,
```

with spontaneous polarization along the electrospinning/poling `z` direction. It also prints

```text
div(eps0 eps_b E + P) = 0.
```

The simulation is 512 x 512 x 512 nm^3 on 128^3 nodes (`4 nm` spacing), with periodic field boundary conditions and an applied field of `1.2e5 V/m`. The supplementary file directly supplies the PVDF coefficients.

### Why it matters

Su 2022 serves two roles:

1. a clean beta-PVDF uniaxial benchmark;
2. an explicit published check of the contracted sixth-order Landau polynomial convention used to expand the Guo S3 coefficient set.

Status: `SECONDARY_EXECUTABLE_BENCHMARK`.

## Revised reproduction order

Now that the full Ahluwalia paper is available, the sequence is no longer limited by missing source text:

```text
Guo 2024 homogeneous Landau/S16 anisotropy
-> Ahluwalia 2008 MD-to-LGD/gradient/noise/time transfer audit
-> Su 2022 beta-PVDF uniaxial benchmark
-> Guo 2024 vector TDGL + weak-anisotropy spiral benchmark
-> Huang/Rui crystal-OAF-IAF substitution at identified coefficients
-> project DFT/MD/MLP generation of unresolved atomistic observables
-> re-fit continuum coefficients and validate against experiments
```

The reason for doing the Ahluwalia transfer audit before the full Guo spiral solver is scientific: it fixes the project architecture for cross-scale parameterization, so later code does not accidentally treat phenomenological coefficients as arbitrary fitting knobs.

## What the eventual MLP bridge is allowed to do

The preferred bridge is observable-based rather than coefficient-by-coefficient black-box regression:

```text
DFT / MLP atomistic trajectories
  -> phase energies and P(T)
  -> wall profiles / excess energies
  -> elastic and dielectric response
  -> fluctuation statistics and relaxation times
  -> calibrated LGD / gradient / electrostrictive / kinetic / noise parameters
  -> phase-field observables
```

A direct MLP prediction of an LGD coefficient is acceptable only when the training target has an unambiguous source definition and the resulting continuum model is independently validated.

## Evidence boundary for crystal / OAF / IAF

Huang 2021 and Rui 2021/2022 remain the experimental sources for crystal/OAF/IAF structure and amorphous dynamics. They do not by themselves provide a complete three-phase TDGL free-energy parameterization. Their data enter only after a published polymer phase-field baseline is reproduced and the exact continuum quantity being replaced is identified.

No `PROJECT_IMAGE_DERIVED_CONSTRAINT`, sinusoidal waviness spectrum, frozen beta/OAF source allocation, or project-defined OAF decay profile is permitted in a literature-reproduction benchmark.
