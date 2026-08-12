# v0.1.2 — BDS dynamics bridge

## Purpose

v0.1.2 converts the previous single-frequency dielectric identifiability test into a workflow that can accept broadband dielectric spectroscopy (BDS) data and connect a fitted relaxation time to an explicit time-domain auxiliary polarization.

The scientific boundary is strict: the 2022 Macromolecules paper analyzes rigid-amorphous-fraction (RAF) devitrification in unpoled and poled BOPVDF. RAF dynamics are useful constraints for the amorphous/interphase part of the model, but RAF is not silently relabeled as OAF.

Primary source: Rui et al., *Macromolecules* 2022, 55, 9705–9714, DOI 10.1021/acs.macromol.2c01110.

The accessible primary-source record supports the following qualitative statements:

- the concentration of active dipoles rises substantially from -30 to 40 °C;
- calculated dipole moment and Kirkwood–Fröhlich g-factor are higher in poled BOPVDF;
- dipole–dipole interaction increases strongly with temperature;
- rotational dipole mobility increases by more than four orders of magnitude over -30 to 40 °C;
- the Supporting Information contains results for n(T), m_d(T), g(T), the three-phase MOAF dielectric calculation, lambda(T), and mu_r(T).

No numerical n(T), m_d(T), g(T), lambda(T), or mu_r(T) curves are entered into the executable parameter database until those SI values are transcribed or digitized with traceable provenance.

## Spectrum-level scaffold

For one Debye process,

```text
epsilon*(omega) = epsilon_inf + Delta_epsilon / (1 + i omega tau)
```

with repository convention

```text
epsilon* = epsilon' - i epsilon''
```

so measured dielectric loss is `epsilon_loss = -Im(epsilon*) > 0`.

`src/pvdf_pf/calibration/bds.py` fits `epsilon_inf`, `Delta_epsilon`, and `tau` jointly to epsilon' and epsilon''. This fit is deliberately labelled a spectrum-level scaffold: real BOPVDF spectra may require distributed/non-Debye processes and paper-specific decomposition.

## Kirkwood–Fröhlich layer

When an independently supported active-dipole number density N is supplied, the code can evaluate

```text
((eps_s-eps_inf)(2 eps_s+eps_inf)) / (eps_s (eps_inf+2)^2)
  = N g mu^2 / (9 eps0 kB T)
```

and therefore constrain the product `g*mu^2`. The code can infer `g` if `mu` is independently supplied, or `mu` if `g` is independently supplied. It does not infer N, g, and mu simultaneously from one dielectric spectrum.

## Time-domain auxiliary polarization

`src/pvdf_pf/physics/dipolar.py` implements

```text
tau dP_rel/dt + P_rel = eps0 Delta_epsilon E_local
```

with an exact exponential step for piecewise-constant local electric field. An OAF mask can restrict this auxiliary response spatially.

This variable is not yet coupled into the main TDGL solver. Current TDGL time is dimensionless, while BDS relaxation times are in seconds. Directly inserting `tau_s` into the TDGL loop would create a false physical time scale.

## Real-data interface

`scripts/fit_bds_csv.py` accepts CSV data with required columns:

```text
frequency_Hz,epsilon_real,epsilon_loss
```

and optional grouping columns:

```text
sample_state,temperature_C
```

It writes one Debye fit per sample-state/temperature group. These outputs remain `SPECTRUM_FIT_ONLY` until the RAF/OAF/MAF decomposition is supplied.

## Exit criteria for v0.1.2

Before promoting BDS parameters into the physical phase-field configuration, all of the following are required:

1. numerical BDS spectra or digitized source curves over multiple frequencies and temperatures;
2. explicit separation of crystal, MAF/IAF, RAF/OAF and electrode/conductive contributions used by the source paper;
3. traceable n(T), m_d(T), g(T), lambda(T), and/or mu_r(T) values from the source/SI where used;
4. a documented mapping between RAF and the OAF variable used in this project, or a decision to model RAF and OAF as separate interphase states;
5. a calibrated physical time mapping for the TDGL order parameter before dynamic coupling.

Until then, the frequency-domain and time-domain modules are validated numerical bridges, not a fully calibrated OAF dynamics model.
