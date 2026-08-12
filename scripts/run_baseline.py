from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from configs.baseline_dimensionless import FIELD, GRID, MORPHOLOGY, PHASE_PARAMS, TDGL
from pvdf_pf.core.field_schedule import triangular_field
from pvdf_pf.core.grid import Grid2D
from pvdf_pf.io.output import save_run
from pvdf_pf.morphology.three_phase import lamellar_three_phase, phase_fractions
from pvdf_pf.physics.landau import build_phase_parameter_maps
from pvdf_pf.solver.tdgl import run_tdgl_scalar


def main() -> None:
    grid = Grid2D(**GRID)
    phase_map = lamellar_three_phase(grid, **MORPHOLOGY)
    maps = build_phase_parameter_maps(phase_map, PHASE_PARAMS)
    field = triangular_field(**FIELD)

    P_final, history, snapshots = run_tdgl_scalar(grid, maps, field, **TDGL)

    outdir = ROOT / "outputs" / "baseline"
    save_run(
        outdir,
        phase_map,
        maps,
        P_final,
        history,
        snapshots,
        metadata={
            "status": "dimensionless numerical smoke test",
            "physical_calibration": False,
            "phase_ids": {"0": "crystal", "1": "OAF", "2": "MAF"},
            "phase_fractions": phase_fractions(phase_map),
            "defect_couplings_active": False,
        },
    )

    print(f"saved to: {outdir}")
    print(
        "Pavg(final)={:.6g}; Pavg range=({:.6g}, {:.6g})".format(
            history["Pavg"][-1], history["Pavg"].min(), history["Pavg"].max()
        )
    )


if __name__ == "__main__":
    main()
