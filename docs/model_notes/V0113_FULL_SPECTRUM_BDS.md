# v0.1.13 - Same-state full-spectrum BDS calibration

## Purpose

v0.1.11 showed that the Figure-5B devitrification fractions and Figure-S2 MOAF permittivity cannot be transferred together as a closed local OAF parameter set for the BOPVDF film state. v0.1.12 therefore moved to the same unpoled/poled BOPVDF state and retained only the combined amorphous response supported by the source.

v0.1.13 adds the frequency dependence by digitizing the alpha-relaxation portions of Rui-2022 Figure 2 and fitting epsilon'(f,T) and epsilon''(f,T) simultaneously.

## Source constraints

The article states that both unpoled and poled BOPVDF have crystallinity near 0.52 and uses the parallel two-phase relation

`epsilon_film = eta_cr*epsilon_cr + (1-eta_cr)*epsilon_am`

with `eta_cr=0.52` and `epsilon_cr=3`. The authors explicitly combine OAF and IAF because the Kirkwood-Frohlich treatment is not independently applicable to the oriented OAF.

Therefore v0.1.13 fits at the film level first, then maps the accepted fit to `combined OAF+IAF amorphous response`. It does not create an OAF-only spectrum.

## Relaxation models

Debye:

`eps*(w)=eps_inf+(eps_s-eps_inf)/(1+i*w*tau)`

Cole-Cole parameterization used here:

`eps*(w)=eps_inf+(eps_s-eps_inf)/(1+(i*w*tau)^beta)`

with `0<beta<=1` and `alpha_CC=1-beta`. Debye is the special case `beta=1`.

The Figure-3A static film permittivity is held fixed. The fitted parameters are `epsilon_infinity`, `tau`, and (for Cole-Cole) `beta`. A non-negative loss floor is allowed only as a nuisance term for residual background loss and is never promoted to a material parameter.

## Main numerical result

For every one of the 16 same-state spectra (-30 to 40 C, unpoled and poled), Cole-Cole is strongly preferred to single-Debye. With the repository-subsampled digitized spectra, `Delta AIC = AIC_Debye - AIC_ColeCole` remains above 50 for every state.

Representative Cole-Cole beta values are approximately:

- unpoled: 0.35 at -30 C, increasing to about 0.50 around 20-30 C;
- poled: 0.29 at -30 C, increasing to about 0.47 around 20 C.

Thus the alpha relaxation is broad; treating it as one ideal Debye relaxation is not supported by Figure 2.

The film-level fitted response is converted to the combined amorphous response by

`epsilon_am*(w,T) = [epsilon_film*(w,T)-0.52*3]/0.48`.

This mapping preserves the same tau and spectral broadening but rescales the dielectric strengths. The result remains an OAF+IAF effective response.

## Raw Figure-2 peak versus Figure-3B deconvolved peak

Figure 3B reports `f_d(T)` obtained from the authors' loss-peak/deconvolution procedure. v0.1.13 keeps that source quantity separate from the characteristic time fitted directly to the raw Figure-2 alpha window. They are compared diagnostically but are not forced to be identical, because the printed spectra contain overlapping background/ionic contributions and the paper describes a deconvolution procedure.

## Measurement ceiling

At 40 C the alpha loss peak approaches the 10^7 Hz upper measurement limit. Both unpoled and poled 40 C states are therefore marked `RIGHT_CENSORED_PEAK`. Their fitted high-frequency permittivities are trend-level extrapolations and must not be treated as strong calibration anchors.

## Promotion boundary

v0.1.13 supports a frequency- and temperature-dependent **combined amorphous** constitutive response. It still does not justify:

- epsilon_OAF*(omega,T) as an independent local material property;
- tau_OAF(T) as an independently measured OAF relaxation time;
- using measured low-frequency epsilon as TDGL `eps_b`;
- direct promotion into a physical-time TDGL without separating fast background and explicit polarization variables.

The next step is to convert the accepted Cole-Cole combined-amorphous response into a causal time-domain representation (Prony/generalized-Debye auxiliary modes) that can be coupled to heterogeneous local fields without using a nonlocal fractional derivative inside the phase-field solver.
