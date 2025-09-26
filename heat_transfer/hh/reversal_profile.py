from heat_transfer.heat_rate import WallResistance, FoulingResistance, TotalResistance
from heat_transfer.fluid_props.GasProps import HotFlueGas
from heat_transfer.hh.afad import HeatLossCalculator, GasTemperature
from heat_transfer.convection.flow import GasVelocityCalculator, GasNusselt, GasReynoldsNumber, PrandtlCalculator
from heat_transfer.convection.convection import GasConvectiveCoefficient

class ReversalChamberRunner:
    def __init__(self,
                 gas: HotFlueGas,
                 heat_loss_calc: HeatLossCalculator,
                 gas_temp: GasTemperature,
                 r_water: WallResistance,           # instance (water-side or external)
                 r_wall: WallResistance,            # instance (wall)
                 r_foul_inner: FoulingResistance,   # instance
                 r_foul_outer: FoulingResistance,   # instance
                 h_rad_func,                        # fn(T,P,X)-> h_rad [W/m^2-K]
                 cp_func,                           # fn(T,P,X)-> cp [J/kg-K]
                 mass_func,                         # fn(T,P,rho,dx)-> mass [kg]
                 h_func, inv_h_func,                # h(T)->J/kg, inv_h(h)->K
                 K_func,                            # fn(x,T,P,Re)-> minor-loss K
                 mdot, A_flow, A_surface, D_h, X,
                 nusselt_n):
        self.gas = gas
        self.hl = heat_loss_calc
        self.gt = gas_temp
        self.r_water = r_water
        self.r_wall = r_wall
        self.r_fi = r_foul_inner
        self.r_fo = r_foul_outer
        self.h_rad = h_rad_func
        self.cp = cp_func
        self.mass = mass_func
        self.h = h_func
        self.inv_h = inv_h_func
        self.K = K_func
        self.mdot = mdot
        self.A_flow = A_flow
        self.A_surface = A_surface
        self.D_h = D_h
        self.X = X
        self.nusselt_n = nusselt_n

    def run(self, L, dx, T0, P0):
        x, T, P = [0.0], [T0], [P0]
        v_list, Re_list, Pr_list, Nu_list = [], [], [], []
        h_conv_list, K_list, dp_list = [], [], []
        qflux_list, qrate_list = [], []

        n = int(L / dx)
        Ti, Pi = T0, P0

        for i in range(n):
            xi = (i + 1) * dx

            k  = self.gas.thermal_conductivity(Ti, Pi, self.X)
            mu = self.gas.viscosity(Ti, Pi, self.X)
            rho= self.gas.density(Ti, Pi, self.X)
            cp = self.cp(Ti, Pi, self.X)

            v = GasVelocityCalculator(self.mdot, rho, self.A_flow).velocity()

            Pr = PrandtlCalculator(mu, cp, k).compute()
            Re = GasReynoldsNumber(rho, v, self.D_h, mu).calculate()
            Nu = GasNusselt(Re, Pr, self.nusselt_n).calculate()
            h_conv = GasConvectiveCoefficient(k, self.D_h, Nu).h()
            h_rad = self.h_rad(Ti, Pi, self.X)

            r_water_pa = self.r_water.resistance_per_area()
            r_wall_pa  = self.r_wall.resistance_per_area()
            r_gas_pa   = 1.0 / (h_conv + h_rad)
            r_fi_pa    = self.r_fi.resistance_per_area()
            r_fo_pa    = self.r_fo.resistance_per_area()
            R_tot = TotalResistance(r_water_pa, r_wall_pa, r_gas_pa, r_fi_pa, r_fo_pa).resistance_per_area()

            q_flux = self.deltaT(xi, Ti, Pi) / R_tot
            q_rate = q_flux * self.A_surface

            heat_loss = self.hl.heat_loss(q_rate, v, dx)
            m = self.mass(Ti, Pi, rho, dx)
            Ti = self.gt.new_temperature(Ti, heat_loss, m, self.h, self.inv_h)

            K = self.K(xi, Ti, Pi, Re)
            dp = K * 0.5 * rho * v * v
            Pi = Pi - dp

            x.append(xi); T.append(Ti); P.append(Pi)
            v_list.append(v); Re_list.append(Re); Pr_list.append(Pr); Nu_list.append(Nu)
            h_conv_list.append(h_conv); K_list.append(K); dp_list.append(dp)
            qflux_list.append(q_flux); qrate_list.append(q_rate)

        return {
            "x": x, "T": T, "P": P,
            "v": v_list, "Re": Re_list, "Pr": Pr_list, "Nu": Nu_list,
            "h_conv": h_conv_list, "K": K_list, "dp": dp_list,
            "q_flux": qflux_list, "q_rate": qrate_list
        }
