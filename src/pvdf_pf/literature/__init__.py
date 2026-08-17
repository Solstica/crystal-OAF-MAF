"""Source-faithful literature reproduction modules.

Modules in this package must distinguish quantities printed by a primary source
from values derived inside this repository. Project-defined morphology or
constitutive hypotheses do not belong here.
"""

from .guo2024 import (
    GUO2024_STRONG_ANISOTROPY_AXIS,
    GUO2024_WEAK_ANISOTROPY_AXIS,
    AxisLandauParameters,
    axis_equilibrium,
    landau_axis_density,
    landau_axis_derivative,
)

__all__ = [
    "AxisLandauParameters",
    "GUO2024_STRONG_ANISOTROPY_AXIS",
    "GUO2024_WEAK_ANISOTROPY_AXIS",
    "landau_axis_density",
    "landau_axis_derivative",
    "axis_equilibrium",
]
