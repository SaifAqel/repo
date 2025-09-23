# transport/htc_water_boil.py
from abc import ABC, abstractmethod
from common.units import Q_, ureg

class IHTCWaterBoilingModel(ABC):
    @abstractmethod
    def h(self, G: Q_, x: Q_, qpp: Q_, P: Q_, Dh: Q_) -> Q_: ...

class ChenBoilingModel(IHTCWaterBoilingModel):
    def h(self, G, x, qpp, P, Dh):
        G_m  = G.to("kg/m^2/s").magnitude
        x_m  = max(x.to("dimensionless").magnitude, 0.0)
        q_m  = abs(qpp.to("W/m^2").magnitude)
        Dh_m = Dh.to("m").magnitude
        # crude liquid-only Dittus–Boelter + Chen-style boost
        mu_l = 1e-3     # Pa*s
        k_l  = 0.6      # W/m/K
        Pr_l = 2.0
        Re_lo = max(G_m*Dh_m/mu_l, 1.0)
        h_lo = 0.023*(Re_lo**0.8)*(Pr_l**0.4)*k_l/Dh_m
        Xtt = ((1-x_m)/max(x_m,1e-6))**0.9 if x_m>0 else 1e6
        F = 1 + 2400/max(x_m**1.16, 1e-6) if x_m>0 else 1.0
        h = h_lo*(1 + 3000*Xtt**-1.16) + F*(q_m**0.79)*1e-6
        return Q_(h, "W/m^2/K")

class ThomModel(IHTCWaterBoilingModel):
    """
    Thom correlation for boiling water. Units:
    q'': W/m^2, P: Pa. Internally uses q [MW/m^2], P [MPa].
    ΔT_sat = 22.5 * q^0.5 * exp(-P/8.7)  => h = q'' / ΔT_sat
    Ref: wikidoc "Heat transfer coefficient" (Thom correlation).
    """
    def h(self, G: Q_, x: Q_, qpp: Q_, P: Q_, Dh: Q_) -> Q_:
        q_MW_m2 = qpp.to("MW/m^2").magnitude
        P_MPa   = P.to("MPa").magnitude
        dT_K = 22.5 * (q_MW_m2 ** 0.5) * np.exp(-P_MPa / 8.7)
        h = (qpp.to("W/m^2").magnitude) / max(dT_K, 1e-9)
        return Q_(h, "W/m^2/K")


class GungorWintertonModel(IHTCWaterBoilingModel):
    """
    Gungor–Winterton (1987) composite: h = S_GW*h_nb + E_GW*h_lo
    h_nb via Cooper pool boiling (W/m^2/K): 55 * p_r^0.12 * (-log10 p_r)^-0.55 * M^-0.5 * q''^0.67
    h_lo via Dittus–Boelter with liquid-only props.
    E_GW = 1 + 24000*Bo^1.16 + 1.37*Xtt^-0.86
    S_GW = [1 + 1.15e-6 * E2 * Re_l^1.17]^-1
    Low-Froude corrections (horizontal): E2 = Fr_l**(0.1 - 2*Fr_l), multiply S_GW by Fr_l**0.5 when Fr_l<0.05.
    Uses water props at ~100°C; adequate for a first-principles implementation given current method signature.
    Refs: Purdue notes and texts reproducing GW equations. 
    """
    # constants for water near 100°C
    _mu_l = 1.0e-3        # Pa*s
    _mu_v = 1.25e-5       # Pa*s
    _k_l  = 0.6           # W/m/K
    _Pr_l = 2.0           # –
    _rho_l = 958.0        # kg/m^3
    _rho_v = 0.6          # kg/m^3
    _h_fg = 2.257e6       # J/kg
    _M    = 18.01528      # kg/kmol
    _Pc   = 22.064e6      # Pa
    _g    = 9.80665       # m/s^2

    def _h_lo(self, G: float, Dh: float) -> float:
        Re_lo = max(G*Dh/self._mu_l, 1.0)
        Nu = 0.023*(Re_lo**0.8)*(self._Pr_l**0.4)
        return Nu*self._k_l/max(Dh, 1e-12)

    def _h_nb_cooper(self, qpp: float, P: float) -> float:
        pr = max(P/self._Pc, 1e-6)
        term = 55.0 * (pr**0.12) * (10.0**(-0.55*np.log10(pr))) * (self._M**-0.5)
        return term * (qpp**0.67)

    def h(self, G: Q_, x: Q_, qpp: Q_, P: Q_, Dh: Q_) -> Q_:
        # scalars
        Gm   = G.to("kg/m^2/s").magnitude
        xm   = min(max(x.to("dimensionless").magnitude, 0.0), 0.999999)
        q    = abs(qpp.to("W/m^2").magnitude)
        Ppa  = P.to("Pa").magnitude
        Dhm  = Dh.to("m").magnitude

        # dimensionless groups
        Bo  = q / max(Gm*self._h_fg, 1e-12)
        Xtt = ((1.0 - xm)/max(xm, 1e-6))**0.9 * (self._rho_l/self._rho_v)**0.5 * (self._mu_l/self._mu_v)**0.1
        Re_l = max(Gm*(1.0 - xm)*Dhm/self._mu_l, 1.0)
        Fr_l = (Gm**2)/(max(self._rho_l,1e-9)**2 * self._g * max(Dhm,1e-12))

        # GW factors
        E_GW = 1.0 + 24000.0*(Bo**1.16) + 1.37*(max(Xtt, 1e-12)**-0.86)
        E2   = Fr_l**(0.1 - 2.0*Fr_l) if Fr_l < 0.05 else 1.0
        S_GW = 1.0/(1.0 + 1.15e-6 * E2 * (Re_l**1.17))
        if Fr_l < 0.05:
            S_GW *= Fr_l**0.5

        # components
        h_lo = self._h_lo(Gm, Dhm)
        h_nb = self._h_nb_cooper(q, Ppa)

        h = S_GW*h_nb + E_GW*h_lo
        return Q_(h, "W/m^2/K")
