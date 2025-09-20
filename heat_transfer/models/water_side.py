# models/water_side.py
from dataclasses import dataclass
@dataclass
class DrumBoilingSide:
    model
    drum_geom
    def h(self,G: float, x: float, qpp: float, P: float, Dh: float) -> float:
        return self.model.h(G,x,qpp,P,Dh)
    def Tsat(self, thermo, P: float) -> float:
        return thermo.water_saturation(P)["Tsat"]
