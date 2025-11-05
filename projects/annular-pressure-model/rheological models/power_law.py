import numpy as np
import math

def powerlaw_gradient(k: float, n: float, v: float, Db: float, Dp: float) -> float:
    Dh = Db - Dp
    return 4.0 * k * (8.0 + 4.0/n)**n * (v**n) / (Dh**(n + 1.0))

def powerlaw_gradient_from_Q(k: float, n: float, Q: float, Db: float, Dp: float) -> float:
    Aa = math.pi * 0.25 * (Db*Db - Dp*Dp)
    v = Q / Aa
    Dh = Db - Dp
    return 4.0 * k * (8.0 + 4.0/n)**n * (v**n) / (Dh**(n + 1.0))

def powerlaw_friction_profile(d, k: float, n: float, Q: float, Db: float, Dp: float):
    d = np.asarray(d, dtype=float)
    if d.size < 2:
        return np.zeros_like(d)
    grad = powerlaw_gradient_from_Q(k, n, Q, Db, Dp)   # scalar
    dd = np.diff(d)
    return np.concatenate(([0.0], np.cumsum(grad * dd)))

def powerlaw_Q_from_gradient(k: float, n: float, dP_dL: float, Db: float, Dp: float) -> float:
    Db2_minus_Dp2 = Db*Db - Dp*Dp
    Dh = Db - Dp
    # Q = π * (1/(k * dPf/dL))^(1/n) * 2^(3 + 2/n) * (4 + 2/n) * (Db^2 - Dp^2) * (Db - Dp)^(1 + 1/n)
    return (math.pi
            * (1.0 / (k * dP_dL))**(1.0/n)
            * (2.0**(3.0 + 2.0/n))
            * (4.0 + 2.0/n)
            * Db2_minus_Dp2
            * (Dh**(1.0 + 1.0/n)))
