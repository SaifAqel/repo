def run_boiler(cfg, props_g, props_w, cantera_yaml):
    # build initial states
    gas = GasState(T=float(cfg.gas_side.inlet_temperature),
                   P=float(cfg.gas_side.inlet_pressure),
                   m_dot=float(cfg.gas_side.mass_flow_rate),
                   X={k: float(v.m) for k,v in cfg.gas_side.composition.items()})
    # pick water state representation you use
    water = WaterState(P=float(cfg.water_side.pressure),
                       m_dot=float(cfg.water_side.mass_flow_rate),
                       h=... )  # from inlet T or quality

    logs = {}

    # PASS 1
    geom1 = PassGeomCalc(cfg.stages.pass1)   # your class
    mc1 = MarchConfig(
        L=float(cfg.stages.pass1.geometry.inner_length),
        eps_wall=float(cfg.stages.pass1.surfaces.inner.emissivity),
        eps_gas=0.2,
        roughness=float(cfg.stages.pass1.surfaces.inner.roughness),
        fouling_i=(float(cfg.stages.pass1.surfaces.inner.fouling_thickness),
                   float(cfg.stages.pass1.surfaces.inner.fouling_conductivity)),
        fouling_o=(float(cfg.stages.pass1.surfaces.outer.fouling_thickness),
                   float(cfg.stages.pass1.surfaces.outer.fouling_conductivity)),
    )
    gas, water, prof1 = continuous_pass(gas, water, geom1, props_g, props_w, mc1)
    logs["pass1"] = prof1

    # REVERSAL 1 (algebraic mix/turn + local ΔP, optional heat loss)
    gas.P -= local_reversal_loss(cfg.stages.reversal1, gas, geom1)
    gas.T -= reversal_heat_loss(cfg.stages.reversal1, gas)  # if modeled
    logs["reversal1"] = {"P": gas.P, "T": gas.T}

    # PASS 2
    geom2 = PassGeomCalc(cfg.stages.pass2)
    mc2 = mc1.__class__(  # same fields, new values
        L=float(cfg.stages.pass2.geometry.inner_length),
        eps_wall=float(cfg.stages.pass2.surfaces.inner.emissivity),
        eps_gas=0.2,
        roughness=float(cfg.stages.pass2.surfaces.inner.roughness),
        fouling_i=(float(cfg.stages.pass2.surfaces.inner.fouling_thickness),
                   float(cfg.stages.pass2.surfaces.inner.fouling_conductivity)),
        fouling_o=(float(cfg.stages.pass2.surfaces.outer.fouling_thickness),
                   float(cfg.stages.pass2.surfaces.outer.fouling_conductivity)),
    )
    gas, water, prof2 = continuous_pass(gas, water, geom2, props_g, props_w, mc2)
    logs["pass2"] = prof2

    # REVERSAL 2
    gas.P -= local_reversal_loss(cfg.stages.reversal2, gas, geom2)
    gas.T -= reversal_heat_loss(cfg.stages.reversal2, gas)
    logs["reversal2"] = {"P": gas.P, "T": gas.T}

    # PASS 3
    geom3 = PassGeomCalc(cfg.stages.pass3)
    mc3 = mc2.__class__(
        L=float(cfg.stages.pass3.geometry.inner_length),
        eps_wall=float(cfg.stages.pass3.surfaces.inner.emissivity),
        eps_gas=0.2,
        roughness=float(cfg.stages.pass3.surfaces.inner.roughness),
        fouling_i=(float(cfg.stages.pass3.surfaces.inner.fouling_thickness),
                   float(cfg.stages.pass3.surfaces.inner.fouling_conductivity)),
        fouling_o=(float(cfg.stages.pass3.surfaces.outer.fouling_thickness),
                   float(cfg.stages.pass3.surfaces.outer.fouling_conductivity)),
    )
    gas, water, prof3 = continuous_pass(gas, water, geom3, props_g, props_w, mc3)
    logs["pass3"] = prof3

    return gas, water, logs
