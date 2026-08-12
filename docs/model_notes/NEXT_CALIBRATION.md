# Next calibration target: crystal / OAF / MAF baseline

The next development step is to replace only the parameters that can be supported by PVDF-specific evidence. Defect couplings remain disabled during this stage.

## Required quantities

### Morphology
- crystal fraction;
- OAF/RAF fraction;
- MAF fraction;
- lamellar / interphase thickness or a defensible reduced geometry.

### Dielectric response
- crystal dielectric response;
- OAF/RAF dielectric response;
- MAF dielectric response;
- frequency and temperature at which each value is defined.

### Polarization response
- bulk P-E/D-E data for the selected PVDF material state;
- any phase-resolved polarization/susceptibility evidence;
- if phase-resolved coefficients are unavailable, document the inverse-fit assumptions explicitly.

## Calibration order

1. Fix the material state to be modeled (PVDF, BOPVDF, P(VDF-TrFE), processing state).
2. Fix morphology fractions and geometry.
3. Fit/assign dielectric contrast.
4. Fit the minimum local constitutive coefficients needed to reproduce the selected bulk response.
5. Verify morphology/mesh sensitivity.
6. Only then activate HHTT or free-volume couplings.

## Acceptance criteria for v0.1

- every non-placeholder physical parameter has a provenance label;
- simulated zero-field and small-field response is numerically stable;
- the homogeneous-limit solver reproduces its analytical/reference limit;
- changing OAF fraction produces a controlled and interpretable change in effective dielectric/polarization response;
- no defect parameter is used to compensate for an incorrect baseline.
