from __future__ import annotations

import traceback

from run_v0115_local_field_coupling import main


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        message = f"{type(exc).__name__}: {exc}".replace("\n", " ")
        print(
            "::error file=scripts/run_v0115_local_field_coupling.py,"
            f"title=v0.1.15 smoke failure::{message}"
        )
        traceback.print_exc()
        raise
