# geometry/drum.py
from dataclasses import dataclass
from math import pi, acos, sin
from common.units import ureg, Q_

@dataclass
class DrumGeom:
    diameter: Q_       # length
    length: Q_         # length
    nozzle_data: dict
    water_level: Q_    # length, measured from bottom (0..diameter)

    def submerged_area(self) -> Q_:
        """
        Wetted internal lateral area of the drum (shell side).
        """
        R = (self.diameter / 2).to("m")
        L = self.length.to("m")
        h = self.water_level.to("m")

        if h.magnitude <= 0:
            return Q_(0.0, "m^2")
        if h >= 2*R:
            return pi * self.diameter.to("m") * L

        theta = 2.0 * acos((R - h) / R)  # rad
        return R * theta * L  # m^2

    def liquid_volume(self) -> Q_:
        """
        Liquid volume in the horizontal drum.
        """
        R = (self.diameter / 2).to("m")
        h = self.water_level.to("m")
        L = self.length.to("m")

        if h.magnitude <= 0:
            return Q_(0.0, "m^3")
        if h >= 2*R:
            return pi * R**2 * L

        theta = 2.0 * acos((R - h) / R)
        A = R**2 * (theta - sin(theta)) / 2.0
        return A * L
