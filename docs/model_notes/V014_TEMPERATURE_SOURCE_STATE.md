# v0.1.4 — Temperature-dependent Rui 2022 source state

## Purpose

v0.1.4 converts the digitized Rui et al. 2022 Supporting Information curves into a queryable temperature-dependent state without inventing missing ROAF/MOAF fractions or a relaxation-time conversion.

For each sample state (`unpoled`, `poled`) and mobility boundary (`l1` stick, `l2` slip), the source state exposes

- active dipole density `n(T)`;
- effective dipole moment `m_d(T)`;
- Kirkwood factor `g(T)`;
- source-model-derived `epsilon_MOAF(T)`;
- interaction parameter `lambda(T)`;
- rotational mobility `mu_r(T)`.

## Interpolation

The digitized source window is -30 to 40 degC.

- `n`, `m_d`, `g`, `epsilon_MOAF`, and `lambda`: piecewise-linear interpolation.
- `mu_r`: piecewise-linear interpolation in `log10(mu_r)` because Figure S3B uses a logarithmic ordinate and the data span more than four decades.
- extrapolation outside the digitized source range is prohibited.

## Devitrification diagnostic

The source states also report

`I_dev(T) = [n(T)-n(-30 C)]/[n(40 C)-n(-30 C)]`.

This is only a normalized diagnostic of the source-reported increase in active dipoles. It must not be used as `eta_MOAF`, `x_MAF`, or any other phase volume fraction.

## Critical source boundary

The uploaded Supporting Information states that the temperature-dependent `xRAF(T)` and `xMAF(T)` used to solve the ROAF/MOAF decomposition come from Figure 5B of the main article. Those Figure 5B curves are not contained in the SI.

Therefore v0.1.4 does **not** yet create temperature-dependent ROAF/MOAF morphology maps.

Likewise, Figure S3 reports rotational mobility `mu_r`, not a Debye relaxation time `tau`. v0.1.4 does not use an ad-hoc relation such as `tau ~ 1/mu_r`.

## What can enter the next model now

The following are source-driven temperature targets:

- `epsilon_MOAF(T)` as a MOAF-specific effective dielectric target;
- `lambda(T)` and `mu_r(T)` as global BOPVDF dynamical descriptors;
- `n(T), m_d(T), g(T)` as film/amorphous-model descriptors.

`epsilon_MOAF(T)` must still not be copied directly into TDGL `eps_b`, because it already contains dipolar response that may be represented explicitly by a polarization variable.

## Next missing dataset

Digitize main-text Figure 5B (`xRAF(T)`, `xMAF(T)`) for unpoled and poled BOPVDF. Once those curves are available, the SI identities can reconstruct

- `eta_ROAF(T)`;
- `eta_MOAF(T)`;
- temperature-dependent ROAF/MOAF masks or probabilistic phase fields.

Only after a separately validated `mu_r -> tau` relation and a physical TDGL seconds mapping are available should the source mobility be coupled to time-domain TDGL dynamics.
