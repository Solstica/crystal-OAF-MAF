# PVDF crystal/OAF/MAF phase-field model design

## 1. Modeling objective

The model describes how semicrystalline PVDF morphology controls local polarization response. The initial target is:

```
crystal / OAF / MAF morphology
        ↓
local dielectric environment
        ↓
polarization evolution
        ↓
P-E response
```

The model does not initially predict breakdown. Defects are introduced only after the three-phase baseline is validated.

## 2. State variables

Primary order parameter:

```
P(r,t)
```

First version:

```
P = Pz
```

Phase field:

```
chi(r) = {crystal, OAF, MAF}
```

Additional future fields:

```
f_HHTT(r)       molecular connection defect field
f_FV(r)         free-volume field
f_strain(r)     internal strain field
```

## 3. Free energy

The total functional is separated into:

```
F = F_local + F_gradient + F_electrostatic
```

Local term:

```
f_local = a2(chi)P^2 + a4(chi)P^4 + a6(chi)P^6
```

Gradient term controls domain/interface width:

```
f_gradient = kappa/2 |grad P|^2
```

Electrostatic contribution will use spatial dielectric response:

```
div(epsilon(r) grad(phi)) = div(P)
```

## 4. Calibration order

Parameters are introduced in the following order:

1. Literature-supported crystal/OAF/MAF phase response.
2. Experimental dielectric and polarization validation.
3. HHTT coupling.
4. Free-volume coupling.
5. Strain/folded-boundary coupling.
6. Atomistic parameter transfer.

## 5. Parameter provenance rule

Every physical parameter must record:

- source paper;
- experimental or computational origin;
- direct value or inferred value;
- uncertainty.

Dimensionless values are only for solver verification.
