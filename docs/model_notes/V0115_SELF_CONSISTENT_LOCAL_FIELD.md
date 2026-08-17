# v0.1.15 - Self-consistent generalized-Debye local-field coupling

## Purpose

v0.1.14 converted the same-state combined-amorphous Cole-Cole response into a positive generalized-Debye bank. v0.1.15 places that auxiliary polarization inside Gauss' law and solves the heterogeneous local electric field self-consistently while the ferroelectric TDGL order parameter remains frozen.

The constitutive scope remains the same as v0.1.13-v0.1.14: the measured BDS response is assigned to the **combined OAF+IAF amorphous fraction**. OAF and IAF/MAF are still not given independent dynamic spectra.

## Algebraic coupling

For each positive Debye mode,

`P_j^n = a_j P_j^(n-1) + b_j eps0 Delta_eps_j E^n`

with

`a_j = exp(-dt/tau_j)` and `b_j = 1-a_j`.

Gauss' law is

`div[eps0 eps_b(r) E^n + sum_j P_j^n] = 0`.

Substitution gives

`div{eps0 [eps_b(r) + mask(r) sum_j b_j Delta_eps_j] E^n + sum_j a_j P_j^(n-1)} = 0`.

Therefore one zero-order-hold step requires one linear heterogeneous Poisson solve. No fixed-point iteration is needed. After solving `E^n`, every auxiliary polarization is recovered from its exact exponential update.

This algebra is useful because the electrostatic self-consistency and the material time integration are no longer separate approximations. The only time-discretization approximation is the already-audited zero-order hold.

## Applied-field treatment

The older depolarization helper solves only the polarization-generated corrector field. v0.1.15 additionally includes the divergence produced by a spatially heterogeneous dielectric background under a uniform applied macroscopic field,

`E(r) = E0 e_z - grad(psi)`.

This is required for crystal/amorphous laminates because a uniform electric field does not generally satisfy continuity of normal displacement across dielectric interfaces.

## Constitutive assignment

At a selected temperature and sample state:

- crystal cells: `epsilon_b = 3.0` (source-model assumption retained from the Rui analysis);
- every non-crystal cell: `epsilon_b = epsilon_infinity_amorphous(T)` from the same-state Figure-2 fit;
- every non-crystal cell: the same positive generalized-Debye bank from v0.1.14;
- crystal cells: no amorphous relaxation bank.

The morphology may still display the labels crystal/OAF/MAF, but v0.1.15 uses `OAF + MAF = 0.48` only as a combined amorphous mask. The OAF label is not assigned an OAF-specific relaxation law.

## Static-limit verification

If `dt >> max(tau_j)` and the previous auxiliary polarization is zero, all `a_j -> 0`. The self-consistent step reduces to an ordinary heterogeneous static dielectric cell with

`epsilon_amorphous,static = epsilon_infinity_amorphous + sum_j Delta_eps_j`.

This provides a direct numerical verification against the repository's independent static finite-volume dielectric solver.

Two ideal laminate orientations are checked:

1. field parallel to the layers: arithmetic (parallel-capacitor) mixing;
2. field normal to the layers: harmonic (series-capacitor) mixing.

For the parallel control geometry, the source relation used in Rui-2022 is recovered:

`epsilon_film = 0.52*3 + 0.48*epsilon_amorphous`.

For the normal geometry, the same phase properties generate a much lower effective permittivity and a larger electric field in the low-permittivity crystal. The difference is a morphology/local-field effect, not a change in fitted material parameters.

## Dynamic reference

For ideal aligned laminates the continuous generalized-Debye bank also gives an analytic frequency-domain reference:

parallel:

`epsilon_eff*(w) = f_c epsilon_c + f_a epsilon_a*(w)`

normal:

`1/epsilon_eff*(w) = f_c/epsilon_c + f_a/epsilon_a*(w)`.

v0.1.15 reports these values at representative frequencies so the subsequent spatial time-domain simulation has a known limiting case.

## Acceptance criteria

v0.1.15 requires:

- uniform media to retain the imposed macroscopic field exactly;
- the coupled zero-order-hold step to reproduce the ordinary Debye update in a uniform amorphous medium;
- the quasi-static coupled limit to match the independent static cell solver for both laminate orientations;
- the 0.52/0.28/0.20 morphology to assign exactly 0.48 of the cells to the combined-amorphous dynamic mask;
- Gauss residuals to remain numerically small.

## Interpretation boundary

v0.1.15 supports the statement that the experimentally constrained combined-amorphous relaxation can feed back on a heterogeneous local electric field in physical time.

It does not support:

- an OAF-only dynamic constitutive law;
- a distinct MAF/IAF spectrum;
- interpretation of individual Prony modes as molecular subpopulations;
- direct coupling to the dimensionless ferroelectric TDGL clock;
- a fully calibrated nonlinear P-E loop.

The next step should add a **prescribed/frozen ferroelectric polarization field** as a second explicit source in Gauss' law and quantify how crystal polarization and amorphous relaxation redistribute the same local field. Only after that electrostatic coupling is audited should a physical-time mapping for TDGL switching be attempted.
