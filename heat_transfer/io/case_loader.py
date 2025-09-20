# io/case_loader.py
import yaml
from typing import Dict, Callable
from heat_transfer.thermo.provider import CoolPropThermoProvider
from heat_transfer.geometry.tube import SingleTubeGeom
from heat_transfer.geometry.tube_bank import TubeBankGeom
from heat_transfer.geometry.reversal import ReversalChamberGeom
from heat_transfer.models.axial_marcher import AxialMarcher
from heat_transfer.models.stage_runner import StageRunner
from heat_transfer.transport.htc_gas import GnielinskiModel, DittusBoelterModel, ChurchillBernsteinModel
from heat_transfer.transport.emissivity import LecknerModel, SmithModel, WSGGM_SimpleModel
from heat_transfer.transport.htc_water_boil import ChenBoilingModel, ThomModel, GungorWintertonModel
from heat_transfer.flow.friction import ColebrookWhite, Churchill, Haaland
from heat_transfer.core.types import GasState
def _pick(model_name: str, table: Dict[str,Callable]):
    if model_name not in table:
        raise ValueError(model_name)
    return table[model_name]()
def load_case(settings_path: str, case_path: str):
    with open(settings_path,"r") as f:
        settings = yaml.safe_load(f)
    with open(case_path,"r") as f:
        case = yaml.safe_load(f)
    thermo = CoolPropThermoProvider(settings["property_backend"]["coolprop"]["backend"])
    htc_gas_map = {"gnielinski": GnielinskiModel, "dittus_boelter": DittusBoelterModel, "churchill_bernstein": ChurchillBernsteinModel}
    eps_map = {"leckner": LecknerModel, "smith": SmithModel, "wsggm_simple": WSGGM_SimpleModel}
    boil_map = {"chen": ChenBoilingModel, "thom": ThomModel, "gungor_winterton": GungorWintertonModel}
    fr_map = {"colebrook_white": ColebrookWhite, "churchill": Churchill, "haaland": Haaland}
    gi = case["inlet_gas"]
    props = thermo.gas_props(gi["T_in"],gi["P_in"],gi["composition"])
    geom_builders = {
        "single_tube": lambda g: SingleTubeGeom(**g),
        "tube_bank": lambda g: TubeBankGeom(**g),
        "reversal_chamber": lambda g: ReversalChamberGeom(**g)
    }
    def gas_inlet():
        A = geom_builders[case["stages"][0]["type"]](case["stages"][0]["geometry"]).flow_area()
        rho = props["rho"]
        u = gi["m_dot"]/(rho*A)
        return GasState(gi["T_in"],gi["P_in"],gi["m_dot"],gi["composition"],rho,props["cp"],props["mu"],props["k"],props["Pr"],u)
    stages = []
    for s in case["stages"]:
        geom = geom_builders[s["type"]](s["geometry"])
        N = s.get("discretization",{}).get("N_cells") or settings["global_numerics"]["N_cells_default"]
        marcher = AxialMarcher(geom,thermo,_pick(s["correlations"]["h_gas"],htc_gas_map),_pick(s["correlations"]["epsilon_gas"],eps_map),_pick(s["correlations"]["h_water_boiling"],boil_map),_pick(s["correlations"]["friction"],fr_map),geom.length(),N,s["flow"],s["fouling"]["R_fg"],s["fouling"]["R_fw"],s["losses"])
        Tsat = thermo.water_saturation(case["drums"][0]["pressure"])["Tsat"]
        stages.append({"runner": StageRunner(geom,marcher), "Tsat": Tsat, "G_water": 1.0})
    return settings, case, {"gas_inlet": gas_inlet, "stages": stages}
