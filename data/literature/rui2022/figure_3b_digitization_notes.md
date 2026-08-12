# Rui 2022 Figure 3B digitization notes

Source: G. Rui et al., *Macromolecules* 2022, 55, 9705-9714, Figure 3B, DOI `10.1021/acs.macromol.2c01110`.

## Quantity

The plotted markers are `ln(f_d)` versus `1000/T`, where `f_d` is the alpha-relaxation loss-peak frequency obtained from the Figure 2 BDS spectra.

The project stores only the temperatures from -30 to 40 degC used by the same-state BOPVDF calibration sequence. `f_d` and `tau=1/(2*pi*f_d)` are derived in code from the digitized `ln(f_d)` values.

## Axis calibration

The embedded Figure 3 raster was extracted from the supplied PDF. Marker centers were digitized from the plotted red/blue source markers after calibrating the vertical axis against the labeled `ln(f_d)=0,5,10,15` ticks. The temperature identity of each point was assigned from `1000/(T+273.15)` and the plotted marker sequence.

A conservative interpretation uncertainty of roughly `+/-0.1 to 0.2` in `ln(f_d)` should be assumed. No guide-curve interpolation is stored in the CSV.

## Internal check

A straight-line fit of the digitized positive-temperature points in `ln(f_d)` versus `1000/T` recovers activation energies close to the values reported in the article text (`57.8 +/- 1.5 kJ/mol` unpoled; `68.3 +/- 3.6 kJ/mol` poled). The check is used only to validate digitization scale and point identity; the project continues to treat the article-reported activation energies as the stronger source values.

## Interpretation boundary

Figure 3B characterizes the alpha relaxation of the BOPVDF amorphous response. It does not separately resolve OAF and IAF kinetics. v0.1.12 therefore assigns these peak frequencies only to the combined amorphous response used in the source two-phase treatment, not to a standalone OAF material law.
