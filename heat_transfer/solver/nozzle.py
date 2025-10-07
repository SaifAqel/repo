# heat_transfer/solver/nozzle.py

import numpy as np
from typing import Tuple
from common.units import Q_
from heat_transfer.calc_ops.stream_with_calc import GasStream
from heat_transfer.fluid_props.GasProps import GasProps
from heat_transfer.config.models import Nozzle

class Nozzle:
    def __init__(self, geom: Nozzle, gas_stream: GasStream, gas_props: GasProps, K: Q_):
        """
        Initialize Nozzle with configuration objects and extract parameters.
        
        :param nozzle_config: Nozzle config object containing diameter and length (length not used yet).
        :param gas_stream: GasStream object for mass flow rate, composition, etc.
        :param gas_props: GasProps object for density calculation.
        :param K: Loss coefficient (dimensionless).
        """
        self.geom = geom
        self.gas_stream = gas_stream
        self.gas_props = gas_props
        self.K = K
        self.D = self.geom.diameter
        self.A = np.pi * (self.D / 2)**2  # Cross-sectional area
        self.m_dot = self.gas_stream.mass_flow_rate

    def apply(self, T_in: Q_, p_in: Q_) -> Tuple[Q_, Q_]:
        """
        Apply nozzle pressure drop.
        
        :param T_in: Inlet temperature.
        :param p_in: Inlet pressure.
        :return: (T_out, p_out) - Temperature remains the same, pressure drops.
        """
        rho = self.gas_props.density(T_in, p_in, self.gas_stream.composition)
        v = self.m_dot / (rho * self.A)
        dp = self.K * 0.5 * rho * v * v
        return T_in, p_in - dp