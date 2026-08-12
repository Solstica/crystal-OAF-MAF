# Rui et al. 2022 SI digitization record

Primary source: Rui et al., *Macromolecules* **2022**, 55, 9705-9714. Supporting Information `ma2c01110_si_001.pdf`.

## Source quantities

The SI provides:

- Figure S1A: active-dipole number density `n(T)` for unpoled and poled BOPVDF.
- Figure S1B: effective dipole moment `m_d(T)`.
- Figure S1C: Kirkwood-Frohlich factor `g(T)`.
- Figure S2: mobile-OAF dielectric constant `epsilon_MOAF(T)`.
- Figure S3A: dipole-dipole interaction parameter `lambda(T)`.
- Figure S3B: rotational dipole mobility `mu_r(T)` for stick (`l=1`) and slip (`l=2`) boundary conditions.

The source model also states `eta_cr = 0.6`, assumes constant `eta_IAF = 0.2`, and uses `epsilon_cr = epsilon_ROAF = 3.0`. These are source-model assumptions/partitions and are not treated as universal local material constants.

## Digitization procedure

The user-supplied SI PDF was rendered at 600 dpi. Marker centroids were extracted from the plotted colors after calibrating the axes from plot borders/ticks:

- Figure S1 and S3A: linear axes.
- Figure S2: linear temperature and dielectric axes.
- Figure S3B: logarithmic `mu_r` axis calibrated from the `10^3` through `10^8` major ticks.

No fitted guide curve was sampled; only plotted marker positions were transcribed.

All rows are therefore tagged:

`DIGITIZED_SOURCE`

rather than `DIRECT_TABULATED`.

## Precision / uncertainty

The stored precision is chosen to preserve the plotted coordinates, not to claim experimental precision. Estimated digitization uncertainty from marker width and rasterization is approximately:

- `n`: +/- 0.01 in the plotted `10^27 dipoles m^-3` scale;
- `m_d`: +/- 0.003 D;
- `g`: +/- 0.005;
- `epsilon_MOAF`: +/- 0.3;
- `lambda`: +/- 0.005;
- `mu_r`: approximately +/- 0.03-0.05 decade (roughly 7-12%).

These are digitization uncertainties only.

## Important source-model boundary

Section S2-S3 explicitly states that Kirkwood-Frohlich theory was not developed for the oriented polymer chains with permanent dipoles in the MOAF. The authors could obtain `epsilon_MOAF`, but could not independently solve `n(T)`, `m_d(T)`, and `g(T)` for the MOAF. Their reported `n`, `m_d`, and `g` therefore describe the broader amorphous/OAF+IAF treatment used in the source model and must not be assigned directly to MOAF without an additional mapping assumption.

## Extracted qualitative checks

The digitized data reproduce the source-level trends:

- `n(T)` increases over -30 to 40 degC;
- poled BOPVDF has larger `n` and `g` than unpoled BOPVDF over the plotted range;
- `epsilon_MOAF` is higher in poled BOPVDF than unpoled BOPVDF;
- all four `mu_r(T)` series increase by more than four orders of magnitude from -30 to 40 degC.

The three CSV files in this directory are now the traceable numerical inputs for v0.1.3.
