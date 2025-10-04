# heat_transfer/config/__init__.py
from dataclasses import replace
from .schemas import Config, Stages
from heat_transfer.calc_ops.stage_with_calc import DrumWithCalc, PassWithCalc

def with_calc(cfg: Config) -> Config:
    st = cfg.stages
    return replace(
        cfg,
        stages=Stages(
            drum=DrumWithCalc(**st.drum.__dict__),
            pass1=PassWithCalc(**st.pass1.__dict__),
            reversal1=st.reversal1,
            pass2=PassWithCalc(**st.pass2.__dict__),
            reversal2=st.reversal2,
            pass3=PassWithCalc(**st.pass3.__dict__),
        ),
    )
