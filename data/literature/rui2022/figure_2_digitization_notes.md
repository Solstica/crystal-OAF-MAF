# Rui 2022 Figure 2 BDS digitization

Source: Rui et al., *Macromolecules* 2022, 55, 9705-9714, Figure 2.

The source PDF contains Figure 2 as an embedded 1458 x 1196 raster image. The four panels were calibrated directly from plot borders/ticks:

- A: unpoled epsilon', x log10(f/Hz) from -2 to 7, y linear; 5 permittivity units = 87 image pixels.
- B: unpoled epsilon'', x log10(f/Hz) from -2 to 7, y log10(epsilon'') from -1 to 1.
- C: poled epsilon', same axis calibration as A.
- D: poled epsilon'', same axis calibration as B.

Curve colors were calibrated from the printed legend. For each temperature, RGB-nearest colored pixels were sampled in a local x-neighborhood. Only the alpha-relaxation window was retained; low-frequency ion/electrode-polarization-dominated regions were intentionally excluded. Points with poor RGB matching were rejected before the stored spectrum was subsampled to 10 representative frequencies per sample/temperature for repository-level reproducibility and CI.

The retained temperatures are -30, -20, -10, 0, 10, 20, 30, and 40 C for both unpoled and poled BOPVDF. Each row keeps separate RGB digitization scores for epsilon' and epsilon'' and uses `provenance=DIGITIZED_SOURCE`.

Approximate coordinate uncertainty is set by the source line thickness/antialiasing (roughly 1-2 pixels). This corresponds to about 0.06-0.12 in epsilon' and a few percent in epsilon'' depending on vertical position. The digitized data are therefore suitable for relaxation-model discrimination and parameter trends, not for claiming higher precision than the printed figure supports.

At 40 C, the alpha peak is close to the 10^7 Hz measurement ceiling. Those fits are explicitly flagged as right-censored and are not strong anchors for epsilon_infinity.
