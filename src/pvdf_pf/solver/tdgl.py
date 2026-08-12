import numpy as np

from pvdf_pf.core.spectral import semi_implicit_gradient_step
from pvdf_pf.physics.electrostatics import solve_depolarization_scalar
from pvdf_pf.physics.landau import local_energy_derivative


def run_tdgl_scalar(
    grid,
    maps: dict[str, np.ndarray],
    field_schedule: np.ndarray,
    *,
    mobility: float = 1.0,
    dt: float = 0.02,
    kappa: float = 1.0,
    init_scale: float = 1e-3,
    seed: int = 0,
    electrostatic_every: int = 5,
    poisson_tol: float = 1e-8,
    poisson_maxiter: int = 400,
    snapshot_every: int = 100,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[int, np.ndarray]]:
    """Run the first scalar-Pz TDGL baseline on a fixed three-phase morphology."""
    if electrostatic_every < 1:
        raise ValueError("electrostatic_every must be >= 1")

    rng = np.random.default_rng(seed)
    P = init_scale * rng.normal(size=grid.shape)
    hist = {"step": [], "Eext": [], "Pavg": [], "Prms": []}
    snapshots: dict[int, np.ndarray] = {}
    Ez_dep = np.zeros_like(P)

    for step, Eext in enumerate(np.asarray(field_schedule, dtype=float)):
        if step % electrostatic_every == 0:
            _, Ez_dep, _ = solve_depolarization_scalar(
                P,
                maps["eps_b"],
                grid,
                eps0=1.0,
                tol=poisson_tol,
                maxiter=poisson_maxiter,
            )

        E_total_z = Eext + Ez_dep
        dF_nongrad = local_energy_derivative(P, maps, E_total_z)
        P = semi_implicit_gradient_step(P, dF_nongrad, grid, mobility, dt, kappa)

        hist["step"].append(step)
        hist["Eext"].append(float(Eext))
        hist["Pavg"].append(float(P.mean()))
        hist["Prms"].append(float(np.sqrt(np.mean(P**2))))

        if snapshot_every and (step % snapshot_every == 0 or step == len(field_schedule) - 1):
            snapshots[int(step)] = P.copy()

    history = {key: np.asarray(value) for key, value in hist.items()}
    return P, history, snapshots
