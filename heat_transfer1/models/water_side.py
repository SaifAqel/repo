# models/water_side.py
from dataclasses import dataclass
from common.units import ureg, Q_

@dataclass
class DrumBoilingSide:
    model
    drum_geom

    def h(self, G: float, x: float, qpp: float, P: float, Dh: float):
        G   = Q_(G,   ureg.kg / (ureg.meter**2 * ureg.second))
        x   = Q_(x,   ureg.dimensionless)
        qpp = Q_(qpp, ureg.watt / ureg.meter**2)
        P   = Q_(P,   ureg.pascal)
        Dh  = Q_(Dh,  ureg.meter)
        h_val = self.model.h(G.magnitude, x.magnitude, qpp.magnitude, P.magnitude, Dh.magnitude)
        return Q_(h_val, ureg.watt / (ureg.meter**2 * ureg.kelvin))

    def Tsat(self, thermo, P: float):
        P = Q_(P, ureg.pascal)
        Tsat_val = thermo.water_saturation(P.magnitude)["Tsat"]
        return Q_(Tsat_val, ureg.kelvin)
