import numpy as np
import math
g = 9.80665 


def hydrostatic_profile(d, z, rho: float):
    z = np.asarray(z, dtype=float)
    return rho * g * (z[0] - z)


def _to_seg_array(x, nseg: int):
    if isinstance(x, (list, tuple, np.ndarray)):
        return np.asarray(x, dtype=float)
    return np.full(nseg, float(x), dtype=float)


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


import numpy as np
import math

def _to_seg_array(x, nseg: int):
    """Broadcast a scalar to length nseg, or pass through an array-like."""
    if isinstance(x, (list, tuple, np.ndarray)):
        return np.asarray(x, dtype=float)
    return np.full(nseg, float(x), dtype=float)

def build_gradients_per_segment(
    model: str,
    nseg: int,
    Db, Dp, Q,            # scalar or array-like (len==nseg)
    mu_p=None, tau_y=None,# Bingham params (scalar or array-like)
    k=None, n=None        # Power-law params (scalar or array-like)
):
 
    Db_s = _to_seg_array(Db,  nseg)
    Dp_s = _to_seg_array(Dp,  nseg)
    Q_s  = _to_seg_array(Q,   nseg)

    grads = np.empty(nseg, dtype=float)
    m = (model or "").lower()

    if m == "bingham":
        mu_p_s  = _to_seg_array(mu_p,  nseg)
        tau_y_s = _to_seg_array(tau_y, nseg)
        for i in range(nseg):
            Aa = math.pi * 0.25 * (Db_s[i]**2 - Dp_s[i]**2)
            v  = Q_s[i] / Aa
            Dh = Db_s[i] - Dp_s[i]
            grads[i] = (48.0 * mu_p_s[i] * v) / (Dh*Dh) + (6.0 * tau_y_s[i]) / Dh

    elif m == "powerlaw":
        k_s = _to_seg_array(k, nseg)
        n_s = _to_seg_array(n, nseg)
        for i in range(nseg):
            Aa = math.pi * 0.25 * (Db_s[i]**2 - Dp_s[i]**2)
            v  = Q_s[i] / Aa
            Dh = Db_s[i] - Dp_s[i]
            grads[i] = 4.0 * k_s[i] * (8.0 + 4.0/n_s[i])**(n_s[i]) * (v**(n_s[i])) / (Dh**(n_s[i] + 1.0))

    else:
        # leave as NaN to make misuse obvious upstream (wrapper will validate)
        grads.fill(np.nan)

    return grads
