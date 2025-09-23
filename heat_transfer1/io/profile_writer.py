# io/profile_writer.py
import pandas as pd
from common.units import ureg, Q_

def _flatten_quantities(obj: dict) -> dict:
    out = {}
    for k, v in obj.items():
        if isinstance(v, Q_):
            out[k] = v.magnitude
            out[f"{k}_unit"] = str(v.units)
        elif isinstance(v, dict):
            # handle nested dicts like composition
            for subk, subv in v.items():
                if isinstance(subv, Q_):
                    out[f"{k}_{subk}"] = subv.magnitude
                    out[f"{k}_{subk}_unit"] = str(subv.units)
                else:
                    out[f"{k}_{subk}"] = subv
        else:
            out[k] = v
    return out

def export_profiles(stage_result, path: str):
    rows = [_flatten_quantities(c.__dict__) for c in stage_result.cells]
    df = pd.DataFrame(rows)
    if path.endswith(".parquet"):
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, index=False)
