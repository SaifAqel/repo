import CoolProp.CoolProp as CP

class CoolPropThermoProvider:
    def __init__(self, fluid_map: dict):
        self._map = fluid_map
    def __getitem__(self, species: str) -> str:
        return self._map[species]
