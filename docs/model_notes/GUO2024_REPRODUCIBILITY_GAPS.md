# Guo 2024 full-solver reproducibility audit

This note separates what is sufficient for reproducing Supplementary Fig. S16's **homogeneous Landau anisotropy** from what is required for an exact reproduction of the published polar-spiral TDGL simulations.

## Source-complete inputs

The main paper and Supplementary Tables S2/S3 provide:

- vector order parameter `P=(Px,Py,Pz)`;
- TDGL functional form;
- homogeneous Landau tensor form;
- strong/weak Landau coefficient sets;
- gradient coefficients `G11`, `G12`, `G44/G44'`;
- elastic constants `c11`, `c12`, `c44`;
- electrostrictive coefficients `Q11`, `Q12`, `Q44`;
- nanodisk radius 345 nm and thickness 10 nm;
- temperature 25 C;
- short-circuit electrical condition;
- finite-element solution statement;
- target textures in main Fig. 3e and Supplementary Figs. S16-S19.

These are sufficient for the source-constrained homogeneous anisotropy audit already implemented.

## Inputs not numerically specified in the available primary files

The Methods define additional quantities that are required by the full PDE but are not assigned numerical values in Supplementary Tables S2/S3:

1. **Dielectric tensor `kappa_ij`.** Eq. (11) defines the electric energy using `kappa_ij`, but S2/S3 do not list it. The paper contains dielectric measurements, but does not state that a particular measured value/frequency is the `kappa_ij` used by the phase-field solver.
2. **Kinetic coefficient `L`.** Eq. (6) defines `L`, but no numerical value or nondimensionalization that removes it is reported. For a purely equilibrium texture its absolute value can rescale time, but it matters for reproducing the reported time-step-dependent dynamics and energy evolution.
3. **Finite-element discretization.** The nanodisk dimensions and FEM method are stated, but element size/order, mesh density/adaptivity and convergence criteria are not reported in the available article/SI.
4. **Initialization for the spontaneous weak-anisotropy spiral.** The published files show the converged texture, but do not specify a complete initial polarization field / random seed / initialization protocol that uniquely reproduces that state.
5. **Mechanical boundary conditions in sufficient numerical detail.** The energy functional and low-field stress scope are described, but exact displacement/traction constraints needed to reproduce each mechanical run are not fully tabulated.

## Consequence

An implementation that inserts a guessed dielectric constant, arbitrary kinetic coefficient, chosen mesh and handcrafted spiral seed may be a useful **model study**, but it is not an exact reproduction of Guo et al.

Therefore the branch adopts the following labels:

```text
Guo homogeneous Landau S16 anisotropy:
    SOURCE_CONSTRAINED_QUALITATIVE_REPRODUCTION

Guo full vector TDGL spiral with missing numerical inputs filled by us:
    NOT_ALLOWED_AS_EXACT_REPRODUCTION
```

## What can still be reproduced without inventing parameters

The following source-faithful work remains possible:

- algebraic Landau anisotropy and high-symmetry energy checks;
- dimensional/unit checks of `G`, `Q`, and `c` tables;
- a PDE implementation whose missing inputs remain explicit symbolic/configuration slots rather than hidden defaults;
- parameter-sensitivity envelopes showing which conclusions are invariant to the unreported inputs, provided these are clearly labeled as sensitivity analysis;
- validation of our numerical machinery against a more fully specified source such as Ahluwalia 2008 before returning to Guo.

## Project decision

The next exact executable benchmark should not fabricate the missing Guo inputs. The project will use Ahluwalia 2008 and the explicit Su 2022 PVDF equations to validate the continuum implementation, while Guo remains the target for weak-anisotropy topological physics. If author code, mesh files or additional numerical details become available later, the full Guo reproduction gate can be reopened.
