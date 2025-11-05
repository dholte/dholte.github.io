import numpy as np
import math
g = 9.80665  # m/s^2


# Hydrostatic (vectorized)
def hydrostatic_profile(d, z, rho: float):
    z = np.asarray(z, dtype=float)
    return rho * g * (z[0] - z)

# ---------- helpers ----------
def _to_seg_array(x, nseg):
    """Return per-segment array of length nseg (broadcast scalar or pass-through list/array)."""
    if isinstance(x, (list, tuple, np.ndarray)):
        arr = np.asarray(x, dtype=float)
        # assume len==nseg
        return arr
    # scalar
    return np.full(nseg, float(x))

def _piecewise_friction(d, seg_id, grads_per_seg):
    """
    Build Pf(d) given per-segment gradients.
    Pf[0]=0; for each interval i: Pf[i] = Pf[i-1] + grad(interval_i)*Δd_i,
    where grad(interval_i) is taken from the segment of the interval (use seg_id at index i).
    """
    d = np.asarray(d, dtype=float)
    seg_id = np.asarray(seg_id, dtype=int)
    dd = np.diff(d)
    grad_interval = grads_per_seg[seg_id[1:]]  # gradient associated with each interval
    Pf = np.concatenate(([0.0], np.cumsum(grad_interval * dd)))
    return Pf

