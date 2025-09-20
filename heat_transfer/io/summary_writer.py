# io/summary_writer.py
import json
import numpy as np


def make_serializable(obj):
    if isinstance(obj, dict):
        return {k: make_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [make_serializable(v) for v in obj]
    if isinstance(obj, (complex, np.complexfloating)):
        return [obj.real, obj.imag]
    return obj

def export_summary(summary: dict, path: str):
    serializable_summary = make_serializable(summary)
    with open(path, "w") as f:
        json.dump(serializable_summary, f, indent=2)
