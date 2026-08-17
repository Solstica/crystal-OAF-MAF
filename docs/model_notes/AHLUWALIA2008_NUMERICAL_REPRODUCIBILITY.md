# Ahluwalia 2008 numerical reproducibility audit

Primary source: Rajeev Ahluwalia et al., *Physical Review B* **78**, 054110 (2008), DOI `10.1103/PhysRevB.78.054110`.

## Purpose

The full paper supplies enough information to reconstruct the **dimensionless TDGL coefficient set and the atomistic-to-continuum scale mapping**, but it does not specify every numerical choice needed to recreate the exact stochastic trace in Fig. 5b point by point. This note fixes that boundary before implementing the spatial solver.

## 1. Source rescaling around Eq. (11)

The source introduces a polarization scale satisfying

```text
P_scale^2 = beta / gamma,
```

and rescales polarization as `P=P_scale(u,v,w)`. With printed Table-II values,

```text
P_scale = 0.1017858554 C m^-2.
```

At 300 K the printed LGD polynomial gives the homogeneous remnant branch

```text
Pr = 0.09932864466 C m^-2,
w_eq = Pr/P_scale = 0.9758590156.
```

The source's initial MD-match state `Pz=0.114 C m^-2` therefore corresponds to

```text
w_initial = 1.1199984472.
```

## 2. Dimensionless thermodynamic/electrostatic coefficients at 300 K

Using the equations printed with Eq. (11), the directly printed Table-II coefficients and the source assumption `chi_xx=chi_yy=chi_zz`, with `chi_zz` evaluated from the inverse homogeneous Landau curvature at `Pr`, gives

```text
alpha'      = 0.04542396975
chi         = 6.08941803e-12 m F
alpha_xx'   = alpha_yy' = 1.7229057576
E'          = 1.1849187758
```

These are `PROJECT_DERIVED_FROM_SOURCE_EQUATIONS_AND_DIRECT_VALUES` checkpoints, not independently tabulated paper values.

## 3. The characteristic length is not the 2.16 nm grid spacing

The source reports the matched 300 K z-noise amplitude

```text
eps_z_tilde = 0.3556
```

and prints the rescaling relation

```text
eps_z_tilde = sqrt[2 k_B T / (beta P_scale^4 xi^3)].
```

Inverting this relation with the printed Table-II coefficients gives

```text
xi = 4.04816146e-10 m = 0.404816146 nm.
```

This is labelled

```text
PROJECT_INFERRED_FROM_SOURCE_EQ11_AND_MATCHED_NOISE.
```

It should **not** be replaced by the separately reported `2.16 nm` TDGL grid spacing. The two scales have different roles. With the inferred `xi`, the source grid/cell correspond to

```text
Delta x = 2.16 nm  ->  Delta x' = 5.33575555
L       = 138.24 nm -> L'       = 341.4883555.
```

The inferred `xi≈0.405 nm` is also numerically close to the ~0.4 nm MD domain-wall/interchain scale reported in the same paper. That proximity is an observation from the reconstruction; the source does not explicitly state that the two quantities are identical.

## 4. Dimensionless gradient and kinetic coefficients

Using `K1=K2=2.108e-8 J m^3 C^-2` and the source's computational choice `K3=K1=K2`, the inferred length produces

```text
K1' = K2' = K3' = 1.3495602114.
```

From the source matched noises

```text
eps_x = 0.0711
eps_y = 1.0669
eps_z = 0.3556,
```

and `eps_x=sqrt(m) eps_z`, `eps_y=sqrt(n) eps_z`, we recover

```text
m = Gamma_x/Gamma_z = 0.039977506
n = Gamma_y/Gamma_z = 9.001687368.
```

The source also matches `t*≈5` in TDGL to ~45 ps in MD, so

```text
1 t* ~= 9 ps.
```

## 5. What remains numerically under-specified

The paper does **not** report a numerical integration time step for Eq. (11), a random seed, or a complete stochastic time-discretization convention sufficient to regenerate the exact point sequence in Fig. 5b. Therefore the repository will not claim

```text
EXACT_POINTWISE_REPRODUCTION_OF_FIG5B.
```

A future spatial solver must instead use the label

```text
PROJECT_NUMERICAL_IMPLEMENTATION_OF_SOURCE_EQUATIONS
```

for its integrator/time-step/seed, while keeping the physics coefficients, noise amplitudes, boundary conditions and physical scale mapping tied to the source.

The valid validation targets are then statistical and physical rather than pixelwise: equilibrium mean polarization, fluctuation scale, approach to equilibrium near `t*~5`, Gauss-law residual, invariance with numerical refinement, and the reported 9 ps mapping.

## 6. Next implementation gate

The next code stage can now implement Eq. (11)-Eq. (12) without confusing the characteristic rescaling length with the grid spacing. The sequence should be

```text
dimensionless coefficient identities            CLOSED
-> deterministic periodic Eq.11/Eq.12 solver
-> numerical convergence of the deterministic solver
-> stochastic noise with declared integrator/dt/seed
-> 300 K statistical MD-match audit
-> switching-film boundary conditions and Fig.6-8 comparison
```

This gives a defensible solver validation route while preserving the source's own reproducibility limits.
