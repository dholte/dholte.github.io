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

def path_builder(segments: list, theta_deg: float, z0: float = 0.0, ds: float = 1.0):
    if ds <= 0:
        raise ValueError("ds must be positive.")

    # initialize starting conditions
    theta = math.radians(theta_deg)   # convert entry angle from degrees to radians
    z = [z0]                          # list creation for elevation values
    s = [0.0]                         # list creation for station values
    s_cum = 0.0                       # cumulative distance travelled along bore

    for seg in segments:
        t = seg.get("type", "").lower()      
    # defaults to empty string if key "type" doesn't exist to prevent crash
    # converts to lowercase to avoid errors due to capitalization
       
      if t == "tangent":
            L = float(seg["length"])       
            if L <= 0:
                raise ValueError("Tangent length must be positive.")
            n = max(1, int(round(L / ds)))        # determine number of steps to divide tangent based on sampling interval
            ds_eff = L / n                        # effective step size
    
            for _ in range(n):
                s_cum += ds_eff                                # advance by step size
                z_next = z[-1] + ds_eff * math.sin(theta)      # calculate next elevation
                s.append(s_cum)
                z.append(z_next)

        elif t == "arc":
            R = float(seg["R"])                    # arc radius
            d_deg = float(seg["delta_deg"])        # deflection angle
            if R <= 0:
                raise ValueError("Arc radius R must be positive.")
                
            d_rad = math.radians(d_deg)
            L = abs(R * d_rad)                     # calculate arc length
            n = max(1, int(round(L / ds)))
            ds_eff = L / n
            kappa = (1.0 / R) * (1.0 if d_rad >= 0 else -1.0) 

            # Along an arc: dθ/ds = kappa; dz/ds = sin(θ)
            for _ in range(n):
                s_cum += ds_eff
                theta_next = theta + kappa * ds_eff               # determine theta at next step based on curvature
                theta_mid = 0.5 * (theta + theta_next)            # midpoint approximation (assumes linear theta change)
                z_next = z[-1] + ds_eff * math.sin(theta_mid)

                s.append(s_cum)
                z.append(z_next)
                theta = theta_next  

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

