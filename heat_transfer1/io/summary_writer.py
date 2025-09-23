# io/summary_writer.py
import json
import numpy as np
from common.units import ureg, Q_

def make_serializable(obj):
    if isinstance(obj, dict):
        return {k: make_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [make_serializable(v) for v in obj]
    if isinstance(obj, (complex, np.complexfloating)):
        return [obj.real, obj.imag]
    if isinstance(obj, Q_):
        return {"value": obj.magnitude, "unit": str(obj.units)}
    if isinstance(obj, (np.generic,)):  # numpy scalars
        return obj.item()
    return obj

def export_summary(summary: dict, path: str):
    serializable_summary = make_serializable(summary)
    with open(path, "w") as f:
        json.dump(serializable_summary, f, indent=2)
