import numpy as np
import math
g = 9.80665  # m/s^2

def hydrostatic_profile(d, z, rho: float):
    z = np.asarray(z, dtype=float)
    return rho * g * (z[0] - z)

def _to_seg_array(x, nseg):
    if isinstance(x, (list, tuple, np.ndarray)):
        arr = np.asarray(x, dtype=float)
        # assume len==nseg
        return arr
    # scalar
    return np.full(nseg, float(x))

def _piecewise_friction(d, seg_id, grads_per_seg):
    d = np.asarray(d, dtype=float)
    seg_id = np.asarray(seg_id, dtype=int)
    dd = np.diff(d)
    grad_interval = grads_per_seg[seg_id[1:]]  # gradient associated with each interval
    Pf = np.concatenate(([0.0], np.cumsum(grad_interval * dd)))
    return Pf

