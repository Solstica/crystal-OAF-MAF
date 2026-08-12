# v0.1 morphology and OAF-relaxation extension

## Purpose

This stage tests two questions that remained unresolved after the first static inverse problem:

1. Does one phase-permittivity set inferred from an ideal flat laminate transfer to tilted or wavy crystal/OAF/MAF morphologies?
2. How much of the low-frequency dielectric response must be assigned to an explicit OAF dipolar relaxation rather than a fixed local permittivity?

No new phase-field Landau coefficient is calibrated in this stage.

## Periodic morphology family

The morphology generator uses an integer-winding phase coordinate

`u = (m_z z/L_z + m_x x/L_x) mod 1`

so tilted lamellae remain exactly periodic on the simulation cell.  A second generator adds a periodic sinusoidal perturbation to create controlled waviness.  These geometries are sensitivity models; they are not claimed to be measured PVDF morphologies.

The scan first infers an OAF effective permittivity from the flat-lamella, field-parallel case using the existing v0.1 assumptions.  That inverse value is then held fixed while morphology is rotated or distorted.  Loss of agreement with the film target is interpreted as non-transferability of the static inverse parameter, not as a solver failure.

## Explicit OAF Debye scaffold

The first dynamic constitutive scaffold is

`epsilon*_OAF(omega) = epsilon_fast + Delta_epsilon / (1 + i omega tau)`.

At v0.1, `Delta_epsilon` and `tau` are not physical parameters.  A single measured value of `epsilon'(10 Hz)` cannot identify both quantities.  The code therefore scans a set of `tau` hypotheses and, for each one, infers the corresponding `Delta_epsilon` required to reproduce the 10 Hz film target in the ideal parallel-laminate limit.

Every result from this operation is labelled `INFERRED_CONDITIONAL`.

## Why this matters for the later TDGL model

The low-frequency measured dielectric constant includes slow dipolar response.  The TDGL electrostatic term uses a background permittivity in

`D = epsilon0 * epsilon_b * E + P`.

Copying the measured 10 Hz phase or film dielectric constant into `epsilon_b` while also evolving `P` would double count part of the dipolar response.  The dynamic OAF module is therefore kept separate from the TDGL background until a frequency-resolved calibration is available.

## Acceptance criteria for the next step

Before promoting v0.1 to a physical TDGL parameter set, the project still needs:

- measured or atomistically constrained OAF/RAF relaxation information over frequency and temperature;
- a defensible fast/background dielectric constant separated from the explicit polarization channel;
- local polarization-energy curvature or switching barriers for crystal and OAF;
- morphology statistics sufficient to replace the ideal periodic laminate family.
