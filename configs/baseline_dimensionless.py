from pvdf_pf.physics.landau import PhasePolynomial

# Numerical verification parameters only. They are not calibrated to PVDF.
PHASE_PARAMS = {
    "crystal": PhasePolynomial(a2=-1.0, a4=1.0, a6=0.2, eps_b=3.0),
    "oaf": PhasePolynomial(a2=-0.35, a4=1.0, a6=0.2, eps_b=5.0),
    "maf": PhasePolynomial(a2=+0.50, a4=1.0, a6=0.2, eps_b=2.5),
}

GRID = {"nz": 64, "nx": 64, "dz": 1.0, "dx": 1.0}
MORPHOLOGY = {"crystal_thickness": 10, "oaf_thickness": 3, "period": 24, "offset": 0}
TDGL = {
    "mobility": 1.0,
    "dt": 0.02,
    "kappa": 0.8,
    "init_scale": 1e-3,
    "seed": 123,
    "electrostatic_every": 10,
    "poisson_tol": 1e-8,
    "poisson_maxiter": 400,
    "snapshot_every": 100,
}
FIELD = {"n_steps": 800, "amplitude": 1.2, "cycles": 1}
