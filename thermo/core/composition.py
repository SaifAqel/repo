from dataclasses import dataclass
from typing import Dict

@dataclass(frozen=True)
class Composition:
    fractions: Dict[str, float]
    basis: str  # "mass" | "mole"

    def to_mole(self, M: Dict[str,float]) -> "Composition":
        if self.basis == "mole": return self
        n = {sp: self.fractions[sp]/M[sp] for sp in self.fractions}
        tot = sum(n.values()); x = {sp: v/tot for sp,v in n.items()}
        return Composition(x, "mole")

    def to_mass(self, M: Dict[str,float]) -> "Composition":
        if self.basis == "mass": return self
        m = {sp: self.fractions[sp]*M[sp] for sp in self.fractions}
        tot = sum(m.values()); w = {sp: v/tot for sp,v in m.items()}
        return Composition(w, "mass")

def mix_molar_mass(x_mol: Dict[str,float], M: Dict[str,float]) -> float:
    return sum(x_mol[sp]*M[sp] for sp in x_mol)
