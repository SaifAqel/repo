from common.units import ureg, Q_

def sensible_heat(mass_flow, cp_mass, inlet_T, T_ref):
    # Expect: mass_flow [kg/s], cp_mass [kJ/(kg*K)] or [J/(kg*K)], temps [K]
    m_dot = (mass_flow).to('kg/s')
    cp = (cp_mass).to('kJ/(kg*K)')       # converts J/(kg*K) → kJ/(kg*K) if needed
    dT = ((inlet_T) - (T_ref)).to('K')
    return (m_dot * cp * dT).to('kW')      # kJ/s → kW

def total_input_heat(fuel_sens, air_sens, power_LHV):
    return ((power_LHV).to('kW') +
            (fuel_sens).to('kW') +
            (air_sens).to('kW')).to('kW')
