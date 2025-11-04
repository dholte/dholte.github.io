"""
annular_pressure.py
Core functions for annular pressure model.
"""

import math      # import common mathematical functions and constants
g = 9.81         # gravitational acceleration (m/s^2)

# ----------------------------
# ANNULUS GEOMETRY
# ----------------------------

def annulus_area(Db: float, Dp: float): -> float
    """Annular cross-sectional area (m²) = π/4 * (Db² - Dp²)."""
    if Db <= 0 or Dp <= 0:
        raise ValueError("Diameter values must be positive.")
    if Db <= Dp:
        raise ValueError("Borehole diameter must exceed pipe outer diameter.")
    return math.pi * 0.25 * (Db**2 - Dp**2)

def hydraulic_diameter(Db: float, Dp: float): -> float
    """Hydraulic diameter for concentric annulus (m) = Db - Dp."""
    if Db <= 0 or Dp <= 0:
        raise ValueError("Diameter values must be positive.")
    if Db - Dp <= 0:
        raise ValueError("Hydraulic diameter must be positive.")
    return Db - Dp

# ----------------------------
# PATH BUILDER
# ----------------------------

def path_builder(segments: list, theta0_deg: float, z0: float = 0.0, ds: float = 1.0):
    if ds <= 0:
        raise ValueError("ds must be positive.")

    # initialize starting conditions
    theta = math.radians(theta0_deg)  # grade angle (radians)
    z = [z0]
    s = [0.0]
    s_accum = 0.0

    for seg in segments:
        t = seg.get("type", "").lower()      
    # converts to lowercase to avoid errors due to capitalization
       
      if t == "tangent":
            L = float(seg["length"])
            if L <= 0:
                raise ValueError("Tangent length must be positive.")
            n = max(1, int(round(L / ds)))
            ds_eff = L / n
            # March along with constant theta
            for _ in range(n):
                s_accum += ds_eff
                # dz/ds = sin(theta)  (vertical component of along-path)
                z_next = z[-1] + ds_eff * math.sin(theta)
                s.append(s_accum)
                z.append(z_next)

        elif t == "arc":
            R = float(seg["R"])
            d_deg = float(seg["delta_deg"])
            if R <= 0:
                raise ValueError("Arc radius R must be positive.")
            d_rad = math.radians(d_deg)
            L = abs(R * d_rad)  # along-path arc length
            n = max(1, int(round(L / ds)))
            ds_eff = L / n
            kappa = (1.0 / R) * (1.0 if d_rad >= 0 else -1.0)  # signed curvature

            # Along an arc: dθ/ds = kappa; dz/ds = sin(θ)
            for _ in range(n):
                # advance by ds_eff with current theta
                s_accum += ds_eff
                # theta advances linearly with s
                theta_next = theta + kappa * ds_eff
                # integrate dz over the small step using midpoint theta (better than Euler)
                theta_mid = 0.5 * (theta + theta_next)
                z_next = z[-1] + ds_eff * math.sin(theta_mid)

                s.append(s_accum)
                z.append(z_next)
                theta = theta_next  # update grade angle

        else:
            raise ValueError(f"Unknown segment type: {t}")

    return s, z

# ----------------------------
# FLOW FUNCTIONS
# ----------------------------

def mean_velocity(Q: float, Aa: float) -> float:
    """Mean fluid velocity in annulus (m/s): v = Q / Aa."""
    if Q < 0:
        raise ValueError("Flow rate cannot be negative.")
    if Aa <= 0:
        raise ValueError("Annulus area must be positive.")
    return 0.0 if Q == 0 else Q / Aa


def reynolds_number(rho: float, v: float, Dh: float, mu: float) -> float:
    """Reynolds number (–): Re = rho * v * Dh / mu."""
    if any(x <= 0 for x in (rho, Dh, mu)) or v < 0:
        raise ValueError("rho, Dh, mu must be > 0 and v >= 0.")
    return 0.0 if v == 0 else (rho * v * Dh) / mu


def ff_smooth(Re: float) -> float:
    """Darcy friction factor (–): laminar 64/Re; turbulent Blasius; linear blend in transition."""
    if Re <= 0:
        return 0.0
    if Re <= 2100.0:           # laminar
        return 64.0 / Re
    if Re >= 3000.0:           # turbulent (Blasius)
        return 0.3164 * (Re ** -0.25)
    # transitional: linear blend for smoothness
    f_lam = 64.0 / 2100.0
    f_turb = 0.3164 * (3000.0 ** -0.25)
    w = (Re - 2100.0) / (3000.0 - 2100.0)
    return (1 - w) * f_lam + w * f_turb

