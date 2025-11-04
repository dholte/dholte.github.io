"""
annular_pressure.py
Core functions for annular pressure model.
"""

import math          # import common mathematical functions and constants
import warnings      # allows non-fatal warning messages
g = 9.80665          # gravitational acceleration (m/s^2)

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
# VELOCITY
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

# ----------------------------
# REYNOLDS
# ----------------------------

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


# ----------------------------
# PRESSURE FUNCTIONS
# ----------------------------

def hydrostatic_profile(d, z, rho: float):
    Ps = []
    z_ref = z[0]
    for zi in z:
        Ps.append(rho * g * (z_ref - zi))
    return Ps

def friction_gradient_newtonian(rho: float, mu: float, Q: float, Db: float, Dp: float):
    if Q < 0:
        raise ValueError("Flow rate Q cannot be negative.")
    Dh = hydraulic_diameter(Db, Dp)
    Aa = annulus_area(Db, Dp)
    v  = 0.0 if Q == 0 else Q / Aa
    Re = reynolds_number(rho, v, Dh, mu)
    if v == 0.0:
        return 0.0, v, Re, Dh  # no flow, no friction
    if Re <= 2100.0:
        grad = 32.0 * mu * v / (Dh * Dh)
    else:
        f = ff_smooth(Re)  # your function
        grad = f * rho * v * v / (2.0 * Dh)
    return grad, v, Re, Dh

def friction_profile(d, rho: float, mu: float, Q: float, Db: float, Dp: float):
    if len(d) < 2:
        return [0.0]
    grad, v, Re, Dh = friction_gradient_newtonian(rho, mu, Q, Db, Dp)
    Pf = [0.0]
    acc = 0.0
    for i in range(1, len(d)):
        ds = d[i] - d[i-1]
        if ds < 0:
            raise ValueError("Displacement must be non-decreasing.")
        acc += grad * ds
        Pf.append(acc)
    return Pf

def pressure_profile_total(d, z, rho: float, mu: float, Q: float, Db: float, Dp: float, g: float = 9.80665):
    Ps = hydrostatic_profile(d, z, rho, g)
    Pf = friction_profile(d, rho, mu, Q, Db, Dp)
    if len(Ps) != len(Pf):
        raise ValueError("Ps and Pf lengths mismatch.")
    Ptot = [ps + pf for ps, pf in zip(Ps, Pf)]
    return Ps, Pf, Ptot


# ---------------------------------
# SINGLE WRAPPER FOR ALL FUNCTIONS
# ---------------------------------

def compute_pressure_profile(segments, theta_deg, z0, ds, rho, mu, Q, Db, Dp):
    # 1) path
    d, z = path_builder(segments, theta_deg=theta_deg, z0=z0, ds=ds)
    # 2) geometry + hydraulics
    Dh = hydraulic_diameter(Db, Dp)
    Aa = annulus_area(Db, Dp)
    v  = 0.0 if Q == 0 else Q / Aa
    Re = reynolds_number(rho, v, Dh, mu)
    # 3) pressures
    Ps = hydrostatic_profile(d, z, rho)
    Pf = friction_profile(d, rho, mu, Q, Db, Dp)
    Ptot = [ps + pf for ps, pf in zip(Ps, Pf)]
    # 4) quick stats
    i_max = max(range(len(Ptot)), key=lambda i: Ptot[i])
    results = {
        "d": d, "z": z,
        "Ps": Ps, "Pf": Pf, "Ptot": Ptot,
        "Dh": Dh, "Aa": Aa, "v": v, "Re": Re,
        "max_P": Ptot[i_max], "max_P_at_d": d[i_max], "max_depth_at_d": z[i_max]
    }
    # advisories
    if 0.0 < v < 0.76:
        warnings.warn("Annular velocity < 0.76 m/s (~150 ft/min).")
    if Re > 2000:
        warnings.warn(f"Re={Re:.0f} exceeds laminar regime; friction losses may be underpredicted.")
    return results



