# Rui 2022 main-text Figure 5B extraction

Source: Rui et al., *Macromolecules* 2022, 55, 9705-9714, DOI 10.1021/acs.macromol.2c01110, main-text Figure 5B.

The source figure plots calorimetric rigid-amorphous (`x_RAF`) and mobile-amorphous (`x_MAF`) fractions for melt-recrystallized BOPVDF, assuming crystallinity `x_c = 0.59`.

## Directly reported values

The article text explicitly reports:

| T (C) | x_RAF | x_MAF | provenance |
|---:|---:|---:|---|
| below -55 | 0.400 | 0.008 | DIRECT_REPORTED_TEXT |
| -45.2 (Tg) | 0.331 | 0.079 | DIRECT_REPORTED_TEXT |
| -30 | 0.244 | 0.166 | DIRECT_REPORTED_TEXT |
| 40 | 0.014 | 0.396 | DIRECT_REPORTED_TEXT |

The sums around 0.41 are consistent with `1 - x_c = 0.41` up to rounding.

## Plot digitization

Intermediate points were read from the source-quality PDF rendering of Figure 5B using the colored source markers. The graph axes were calibrated with the labeled temperature and fraction ticks. The direct text values above override approximate graph readings at matching temperatures.

No smoothing curve was digitized as if it were raw data. The CSV stores marker-level values and retains `DIGITIZED_SOURCE` provenance.

## Raw calorimetry versus SI dielectric-model approximation

Do not enforce `x_RAF + x_MAF = 0.40` on the main-text Figure 5B data. The main figure uses measured `x_c = 0.59`, so its amorphous total is about 0.41. In Supporting Information Section S2, the dielectric inversion uses rounded structural assumptions (`eta_cr = 0.6`, `eta_IAF = 0.2`, `eta_OAF ~= 0.2`). Those are a separate source-model approximation.

v0.1.5 therefore preserves both levels:

1. **raw mobility-defined calorimetric fractions** from Figure 5B (`x_RAF`, `x_MAF`);
2. **rounded OAF/IAF dielectric-model fractions** used only when reproducing the SI algebra.

RAF/MAF are mobility-defined and OAF/IAF are structure-defined. They must not be renamed into one another.
