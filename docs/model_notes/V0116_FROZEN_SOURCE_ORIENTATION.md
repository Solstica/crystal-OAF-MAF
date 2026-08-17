# v0.1.16 frozen-source orientation audit

## Purpose

v0.1.16 adds a prescribed ferroelectric polarization source to the v0.1.15 self-consistent generalized-Debye electrostatics while keeping the TDGL order parameter frozen.  The version is an orientation/electrostatic audit.  It does not introduce ferroelectric switching kinetics or claim a phase-resolved remanent-polarization calibration.

The solved constitutive relation is

```text
D = eps0 * eps_b(r) * E + P_rel(r,t) + P_FE_frozen(r)
div(D) = 0
```

The generalized-Debye bank remains the Rui 2022 combined OAF+IAF response already calibrated in v0.1.13-v0.1.15.

## Huang 2021 SAXS geometry

Primary source: Yanfei Huang et al., *Nature Communications* 12, 675 (2021), Supplementary Notes 10-11 and Supplementary Fig. 13.

The SI reports the highly poled BOPVDF lamellar period and the source-model layer decomposition:

| quantity | source value | project use | status |
|---|---:|---|---|
| lamellar period | 11.8 nm | reference period | DIRECT_REPORTED_TEXT |
| beta-crystal thickness | 5.78 nm | geometry ratio | SOURCE_MODEL_DERIVED |
| total OAF thickness | 3.02 nm | geometry ratio | SOURCE_MODEL_DERIVED |
| OAF thickness per crystal side | 1.51 nm | profile support | SOURCE_MODEL_DERIVED |
| IAF thickness | 3.00 nm | geometry ratio | SOURCE_MODEL_DERIVED |
| beta crystallinity | 0.52 | retained as WAXD weight fraction | DIRECT_REPORTED_TEXT |
| minimum OAF fraction used in the SI construction | 0.25 | retained as source-model weight fraction | SOURCE_MODEL_DERIVED |

The executable voxel fractions are obtained from the thickness ratios, not from the WAXD weight fraction:

```text
f_beta,vol = 5.78 / 11.8 = 0.4898305085
f_OAF,vol  = 3.02 / 11.8 = 0.2559322034
f_IAF,vol  = 3.00 / 11.8 = 0.2542372881
```

This distinction matters because the SI equations use phase densities when converting the WAXD/OAF weight fractions into layer thicknesses.  Substituting 0.52 directly as voxel occupancy would mix weight and geometric fractions.

## Polarization direction relative to the lamellae

Supplementary Fig. 13 is an edge-on SAXS pattern with the X-ray beam along TD and MD vertical.  The one-dimensional SAXS profile used for the 11.8 nm period is integrated along MD.  The poled BOPVDF has its macroscopic dipole moment along the film-normal direction.  v0.1.16 therefore treats the idealized Huang lamellar normal as MD and the scalar ferroelectric polarization as film-normal `Pz`.

Within this idealized orientation, `Pz` is tangent to the flat beta/OAF/IAF interfaces.  A piecewise `Pz(x)` then satisfies

```text
div(Pz e_z) = dPz/dz = 0
```

although `Pz` changes across the interfaces in `x`.  Consequently an ideal perfectly aligned flat laminate does not acquire a depolarization field from this scalar frozen source alone.  Tilted or wavy lamellae give the interface normal a film-normal component and can generate bound charge.  This orientation statement is an `INFERRED_GEOMETRY` consequence of the reported SAXS/poling directions; it is tested explicitly rather than treated as a fitted material parameter.

## Frozen source magnitude

The available measurements constrain bulk remanent/spontaneous polarization but do not uniquely identify separate beta-crystal and OAF remanent polarizations.  v0.1.16 therefore uses

```text
P_beta = 0.010 C/m^2
mean(P_OAF) = 0.006 C/m^2
P_IAF = 0
```

only as normalized sensitivity amplitudes.  Their status is `PLACEHOLDER_SENSITIVITY`; numerical local fields may be rescaled in the present linear electrostatic model, but these amplitudes must not be reported as measured local polarizations.

A later version can replace them after a temperature/state-matched phase allocation is identified.  The high-field OAF spontaneous-polarization estimates in Huang/Rui are not silently reused as zero-field remanent polarization.

## OAF spatial profile

Two OAF profiles are compared at fixed discrete OAF mean:

1. `uniform`: constant OAF frozen polarization;
2. `interface_decay`: exponential magnitude decay from each beta/OAF interface toward the IAF side.

The decay profile is `PROJECT_SPATIAL_HYPOTHESIS`.  Huang Supplementary Note 11 supports stronger crystal-proximal orientational constraint qualitatively, but its MD box used enlarged lateral beta-cell dimensions (`a=2.1 nm`, `b=1.2 nm` versus approximately `a0=0.86 nm`, `b0=0.49 nm`) because the denser systems did not relax within the 2 ns simulation window.  The MD calculation therefore does not calibrate an absolute OAF decay length or local polarization magnitude.

## Constitutive transfer limitation

v0.1.16 combines the Huang 2021 SAXS geometry with the Rui 2022 generalized-Debye bank established in the previous versions.  This is a cross-source geometry/constitutive transfer used only to test orientation and electrostatic coupling.  It is not evidence that the two experiments provide one phase-resolved parameter set.

All non-crystal cells still share the same combined OAF+IAF relaxation bank.  Independent OAF and IAF relaxation spectra remain unresolved.

## Numerical checks

`tests/test_v0116_frozen_source.py` requires:

- zero frozen source reproduces the v0.1.15 step;
- the ideal Huang-aligned laminate with film-normal `Pz` produces no spurious bound-charge field;
- the same source in a normal-to-layer control generates a finite internal field while satisfying Gauss' law;
- the interface-decay profile preserves the requested discrete mean OAF polarization.

The face-flux Gauss residual is the conservation check.  The separately reported cell-centered constitutive-flux mismatch can remain finite at discontinuous interfaces because it compares a face-flux discretization with a cell-centered reconstruction.
