"""
annular_pressure.py
Core functions for annular pressure model.
"""

import math          # import common mathematical functions and constants
import warnings      # allows non-fatal warning messages
g = 9.81             # gravitational acceleration (m/s^2)

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
    d = [0.0]                         # list creation for horizontal displacement values
    d_cum = 0.0                       # cumulative horizontal displacement

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
                d_cum += ds_eff                                # advance by step size
                z_next = z[-1] + ds_eff * math.sin(theta)      # calculate next elevation
                d.append(d_cum)
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
            dtheta = d_rad / n                     # per step change in theta

            for _ in range(n):
                d_cum += ds_eff           
                theta_mid = theta + 0.5 * dtheta                  # midpoint approximation (assumes linear theta change)
                z_next = z[-1] + ds_eff * math.sin(theta_mid)
                d.append(d_cum)
                z.append(z_next)
                theta += dtheta  

        else:
            raise ValueError(f"Unknown segment type: {t}")

    return d, z

# ----------------------------
# FLOW FUNCTIONS
# ----------------------------

def annular_velocity(Q: float, Db: float, Dp: float) -> float:
    if Q <= 0: 
        raise ValueError("Volumetric flow rate must be positive.")
    v = Q / annulus_area(Db,Dp)
    if v < 0.762:                # recommended minimum annular velocity of 150 ft/min = 0.762 m/s                                                
        warnings.warn(
            f"Annular velocity of {v: .2f} is below the recommended "
            f"minimum of 0.762 m/s for effective cuttings removal.",
            UserWarning,
        )            


def target_velocity_flow(v_target: float, Db: float, Dp: float) -> float:
    if v_target <= 0:
        raise ValueError("Target annular velocity must be positive.")
    if v_target < 0.762:
        warnings.warn(
            f"Target annular velocity of {v_target:.2f} m/s is below the "
            f"recommended minimum of 0.762 m/s for effective cuttings removal.",
            UserWarning,
        )
    return v_target * annulus_area(Db, Dp)


def dynamic_reynolds(rho: float, v: float, Dh: float, mu: float) -> float:
    if rho <= 0:
        raise ValueError("Density must be positive.")
    if Dh <= 0:
        raise ValueError("Hydraulic diameter must be positive.")
    if mu <= 0:
        raise ValueError("Dynamic viscosity must be positive.")
    if v <= 0:
        raise ValueError("Velocity must be positive.")
    return (rho * v * Dh) / mu


def kinematic_reynolds(v: float, Dh: float, nu: float) -> float:
    if Dh <= 0:
        raise ValueError("Hydraulic diameter must be positive.")
    if nu <= 0:
        raise ValueError("Kinematic viscosity must be positive.")
    if v <= 0:
        raise ValueError("Velocity must be positive.")
    return (v * Dh) / nu



def ff_smooth(Re: float) -> float:
    if Re <= 0:
        raise ValueError("Reynolds number must be positive.")
    if 0 < Re <= 2100.0:           
        return 64.0 / Re
    if Re >= 3000.0:           
        return 0.3164 * (Re ** -0.25)
    # transitional: linear blend for smoothness
    f_lam = 64.0 / 2100.0
    f_turb = 0.3164 * (3000.0 ** -0.25)
    w = (Re - 2100.0) / (3000.0 - 2100.0)
    return (1 - w) * f_lam + w * f_turb

