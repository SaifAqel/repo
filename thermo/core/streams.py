from dataclasses import dataclass
from typing import Optional
from .composition import Composition, mix_molar_mass
from common.units import ureg, Q_

@dataclass(frozen=True)
class GasStream:
    T: Q_               # temperature [K]
    P: Q_               # pressure [Pa]
    composition: Composition       # basis set externally
    flow_mass: Optional[Q_] = None   # [kg/s]
    flow_mol: Optional[Q_] = None    # [mol/s]

    def as_mole_fraction(self, M) -> dict:
        return self.composition.to_mole(M).fractions

    def as_mass_fraction(self, M) -> dict:
        return self.composition.to_mass(M).fractions

    def molar_flow(self, M) -> Q_:
        if self.flow_mol is not None:
            return self.flow_mol.to(ureg.mole/ureg.second)
        x = self.as_mole_fraction(M)
        return self.mass_flow(M) / mix_molar_mass(x, M)

    def mass_flow(self, M) -> Q_:
        if self.flow_mass is not None:
            return self.flow_mass.to(ureg.kg/ureg.second)
        x = self.as_mole_fraction(M)
        return self.flow_mol * mix_molar_mass(x, M)
