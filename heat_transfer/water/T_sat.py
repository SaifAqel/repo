from iapws import IAPWS97

class WaterSaturation:
    def temperature_from_pressure(self, P_MPa: float) -> float:
        return IAPWS97(P=P_MPa, x=0).T
