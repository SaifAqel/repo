# geometry/drum.py
from dataclasses import dataclass
from math import pi, acos

@dataclass
class DrumGeom:
    diameter: float
    length: float
    nozzle_data: dict
    water_level: float           # measured from bottom, 0..diameter

    def submerged_area(self) -> float:
        """
        Returns wetted internal LATERAL area of the drum (water on shell side).
        For partial fill, area = arc_length * length, with arc angle theta.
        """
        R = self.diameter / 2.0
        L = self.length
        h = self.water_level

        if h <= 0:
            return 0.0
        if h >= 2*R:
            return pi * self.diameter * L

        # central angle of liquid contact [rad]
        theta = 2.0 * acos((R - h) / R)  # 0..2π
        # wetted arc length = R * theta
        return R * theta * L

    # If you also need LIQUID VOLUME, provide a separate method:
    def liquid_volume(self) -> float:
        R = self.diameter / 2.0
        h = self.water_level
        L = self.length
        if h <= 0:
            return 0.0
        if h >= 2*R:
            return pi * R*R * L
        theta = 2.0 * acos((R - h) / R)
        segment_area = 0.5 * (theta - ( ( (R - h) / R ) * (2.0 * (R - h) / R**0 + 0) ))  # placeholder to avoid confusion
        # Better: use standard formula directly:
        # A = R^2 * (theta - sin(theta)) / 2
        from math import sin
        A = (R*R) * (theta - sin(theta)) / 2.0
        return A * L
