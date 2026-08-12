from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Grid2D:
    """Periodic x-z grid used by the first scalar-Pz model."""

    nz: int
    nx: int
    dz: float = 1.0
    dx: float = 1.0

    @property
    def shape(self) -> tuple[int, int]:
        return (self.nz, self.nx)

    def kgrid(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        kz = 2.0 * np.pi * np.fft.fftfreq(self.nz, d=self.dz)
        kx = 2.0 * np.pi * np.fft.rfftfreq(self.nx, d=self.dx)
        kz2, kx2 = np.meshgrid(kz, kx, indexing="ij")
        return kz2, kx2, kz2**2 + kx2**2
