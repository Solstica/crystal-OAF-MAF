# v0.1.12 same-state BOPVDF BDS bridge

## Decision

v0.1.11 showed that Figure 5B devitrification weights and Figure S2 mobile-OAF permittivity cannot be transferred together into the closed BOPVDF OAF state without an irreducible inconsistency.

v0.1.12 therefore returns to the main-text BDS data measured on the same unpoled/poled BOPVDF film states.

## Source-supported observables

For each film state and temperature from -30 to 40 degC:

1. Figure 3A supplies the static/low-frequency film permittivity `epsilon_c(T)`.
2. Figure 3B supplies the alpha-relaxation peak `ln(f_d)`.
3. The article states `eta_cr=0.52` for both films and uses `epsilon_cr=3` in Eq. 20.
4. The source two-phase relation is

   `epsilon_c = eta_cr*epsilon_cr + (1-eta_cr)*epsilon_am`.

This yields an inferred `epsilon_am(T)` for the combined amorphous response.

## What is deliberately not inferred

The article explicitly explains that an ideal crystal/OAF/IAF three-phase treatment would be preferable, but the Kirkwood-Frohlich treatment could not be applied to OAF separately. The authors therefore combine OAF and IAF in the amorphous response.

Accordingly v0.1.12 does not create:

- `epsilon_OAF(T)` from Figure 3A;
- `tau_OAF(T)` from Figure 3B;
- an OAF-only `n(T)`, `m_d(T)`, or `g(T)`;
- a TDGL background-permittivity map from these measured low-field dielectric quantities.

The response object is explicitly labeled `combined_OAF_plus_IAF_amorphous_response`.

## Dynamic quantity

The alpha-relaxation time is

`tau_alpha(T) = 1 / (2*pi*f_d(T))`.

The digitized Figure 3B positive-temperature points are fitted only as a digitization check. The resulting Arrhenius activation-energy scale is required to remain close to the stronger values reported directly in the article text:

- unpoled: `57.8 +/- 1.5 kJ/mol`;
- poled: `68.3 +/- 3.6 kJ/mol`.

The text-reported values remain the preferred physical source values.

## Complex-response bridge

The present source set fixes the static combined-amorphous response and the alpha peak time scale, but not an independently validated high-frequency amorphous plateau for every temperature. Therefore absolute `epsilon_am*(omega,T)` is still underdetermined.

`same_state_bds.py` contains a one-Debye numerical bridge only for cases where `epsilon_amorphous_fast` is supplied independently:

`epsilon_am* = epsilon_am,fast + (epsilon_am,static-epsilon_am,fast)/(1+i*omega*tau_alpha)`.

This is an executable coupling form, not a claim that the measured alpha process is exactly single-Debye.

## Promotion rule

v0.1.12 is a calibration/identifiability layer. Nothing in it is promoted automatically to phase-local TDGL `eps_b`, OAF Landau coefficients, or OAF mobility.

The next admissible step is to constrain the fast-limit/relaxation-strength part of the same-state BDS response from Figure 2 or another direct spectrum source, then test whether any OAF/IAF partition is identifiable without importing the incompatible Figure 5B/S2 construction.
