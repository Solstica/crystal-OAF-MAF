from pathlib import Path
import json
import numpy as np


def save_run(
    outdir,
    phase_map: np.ndarray,
    maps: dict[str, np.ndarray],
    P_final: np.ndarray,
    history: dict[str, np.ndarray],
    snapshots: dict[int, np.ndarray],
    metadata: dict,
) -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(out / "state_final.npz", phase_map=phase_map, P=P_final, **maps)
    np.savetxt(
        out / "pe_history.csv",
        np.column_stack([history["step"], history["Eext"], history["Pavg"], history["Prms"]]),
        delimiter=",",
        header="step,Eext,Pavg,Prms",
        comments="",
    )
    if snapshots:
        np.savez_compressed(out / "P_snapshots.npz", **{f"step_{k:06d}": v for k, v in snapshots.items()})
    with open(out / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)
