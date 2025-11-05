import numpy as np
import math
g = 9.80665 


def hydrostatic_profile(d, z, rho: float):
    z = np.asarray(z, dtype=float)
    return rho * g * (z[0] - z)


def _to_seg_array(x, nseg):
    if isinstance(x, (list, tuple, np.ndarray)):
        arr = np.asarray(x, dtype=float)
        return arr
    return np.full(nseg, float(x))


def _piecewise_friction(d, seg_id, grads_per_seg):
    d = np.asarray(d, dtype=float)
    seg_id = np.asarray(seg_id, dtype=int)
    dd = np.diff(d)
    grad_interval = grads_per_seg[seg_id[1:]]  # gradient associated with each interval
    Pf = np.concatenate(([0.0], np.cumsum(grad_interval * dd)))
    return Pf


def total_pressure_profile(d, z, rho: float, Pf):
    d = np.asarray(d, dtype=float)
    z = np.asarray(z, dtype=float)
    Pf = np.asarray(Pf, dtype=float)
    Ps = hydrostatic_profile(d, z, rho)
    Ptot = Ps + Pf
    return Ps, Pf, Ptot


def total_pressure_from_gradients(d, z, rho: float, seg_id, grads_per_seg):
    Pf = _piecewise_friction(d, np.asarray(seg_id, dtype=int), np.asarray(grads_per_seg, dtype=float))
    Ps = hydrostatic_profile(d, z, rho)
    Ptot = Ps + Pf
    return Ps, Pf, Ptot

