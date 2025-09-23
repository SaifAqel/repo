# common/numcheck.py
import numpy as np
def assert_real_finite_q(q, name, atol_im=1e-14):
    m = np.asarray(q.to_base_units().magnitude)
    if np.iscomplexobj(m) and np.max(np.abs(np.imag(m))) > atol_im:
        raise ValueError(f"{name} complex; max|Im|={np.max(np.abs(np.imag(m)))}")
    if not np.all(np.isfinite(np.real(m))):
        raise ValueError(f"{name} non-finite")
    return q  # no mutation
