# v0.1.14 - Cole-Cole to positive generalized-Debye time-domain bridge

## Purpose

v0.1.13 established that the same-state BOPVDF alpha relaxation is broad and is described substantially better by Cole-Cole than by one Debye mode. A direct fractional Cole-Cole constitutive law is inconvenient inside a heterogeneous local-field time integrator because it introduces memory/fractional derivatives.

v0.1.14 therefore approximates the accepted combined-amorphous Cole-Cole susceptibility with a finite positive bank of ordinary Debye modes:

`chi_am*(omega,T) ~= sum_j Delta_eps_j(T)/(1+i*omega*tau_j(T))`

with

`Delta_eps_j >= 0`

and

`sum_j Delta_eps_j = Delta_eps_am(T)`.

The resulting auxiliary variables satisfy local first-order equations:

`tau_j dP_j/dt + P_j = eps0 Delta_eps_j E_local`.

This provides a causal time-domain representation without introducing a fractional derivative into the phase-field/local-field solver.

## Scope inherited from v0.1.13

The input response is the same-state **combined OAF+IAF amorphous response** obtained from Rui-2022 Figure 2 and Figure 3A using `eta_cr=0.52` and `epsilon_cr=3`.

The Prony conversion does not make the response OAF-specific. The individual numerical modes also do not correspond to ROAF, MOAF, IAF, chain segments, or independently identified molecular processes.

## Numerical construction

The target Cole-Cole response is sampled over `1e-2 -- 1e9 Hz`, which brackets the source `1 -- 1e7 Hz` BDS window. Forty-nine Debye modes are distributed logarithmically in characteristic frequency, with a 1.5-decade support margin on both sides.

Non-negative least squares fits real and loss components simultaneously. A static-strength constraint is included and the final positive weights are normalized exactly so the total modal dielectric strength equals the v0.1.13 amorphous dielectric strength.

This makes the approximation passive at the linear-response level: no negative relaxation strengths are introduced simply to improve the fit.

## Time stepping

Each auxiliary Debye mode uses the exact exponential zero-order-hold update already used by the repository's single-mode dipolar solver. Therefore numerical material approximation and finite-time-step discretization remain separately auditable:

1. Cole-Cole -> positive Prony error;
2. continuous Prony -> discrete zero-order-hold error.

The latter is not interpreted as a material loss or relaxation mechanism.

## Acceptance criteria

For all 16 same-state Figure-2 conditions (unpoled/poled, -30 to 40 C), v0.1.14 requires:

- non-negative modal dielectric strengths;
- static dielectric-strength closure to numerical precision;
- normalized maximum complex Cole-Cole/Prony error below `5e-4` over `1e-2 -- 1e9 Hz`;
- a separately reported zero-order-hold error at the Cole-Cole characteristic frequency.

The 40 C high-frequency censoring warning from v0.1.13 remains inherited. A numerically accurate Prony representation cannot add information that is absent beyond the BDS measurement ceiling.

## Promotion boundary

v0.1.14 creates a time-domain auxiliary polarization model for the **combined amorphous response** only. It still does not justify:

- assigning the mode bank only to OAF cells;
- interpreting individual modes as structural phase fractions;
- replacing TDGL `eps_b` with the measured low-frequency permittivity;
- coupling this physical-time bank directly to the dimensionless ferroelectric TDGL clock.

The next development step is a heterogeneous electrostatic coupling test in which the generalized-Debye polarization enters Gauss/Poisson self-consistently as an explicit polarization source while the ferroelectric TDGL order parameter remains frozen or externally prescribed. This isolates local-field feedback before any attempt to calibrate the TDGL physical time scale.
