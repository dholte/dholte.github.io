import numpy as np
import math

# for sensitivity analysis on v
def bingham_gradient(mu_p: float, tau_y: float, v: float, Db: float, Dp: float) -> float:
    Dh = Db - Dp
    return (48.0 * mu_p * v) / (Dh * Dh) + (6.0 * tau_y) / Dh

# most HDD designs begin from known Q
def bingham_gradient_from_Q(mu_p: float, tau_y: float, Q: float, Db: float, Dp: float) -> float:
    Aa = math.pi * 0.25 * (Db*Db - Dp*Dp)
    v = Q / Aa
    Dh = Db - Dp
    return (48.0 * mu_p * v) / (Dh * Dh) + (6.0 * tau_y) / Dh

# cumulative frictional pressure along displacement
def bingham_friction_profile(d, mu_p: float, tau_y: float, Q: float, Db: float, Dp: float):
    d = np.asarray(d, dtype=float)
    if d.size < 2:
        return np.zeros_like(d)
    grad = bingham_gradient_from_Q(mu_p, tau_y, Q, Db, Dp) 
    dd = np.diff(d)
    return np.concatenate(([0.0], np.cumsum(grad * dd)))

# solves for Q based on target gradient for back-calculation or calibration
def bingham_Q_from_gradient(mu_p: float, tau_y: float, dP_dL: float, Db: float, Dp: float) -> float:
    Db2_minus_Dp2 = Db*Db - Dp*Dp
    Dh = Db - Dp
    term1 = (math.pi / (192.0 * mu_p)) * dP_dL * Db2_minus_Dp2 * (Dh * Dh)
    term2 = (math.pi / (32.0  * mu_p)) * tau_y   * Db2_minus_Dp2 * Dh
    return term1 - term2

