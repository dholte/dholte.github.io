import numpy as np
import math

def build_path(
    segments: list,
    theta0_deg: float,
    z0: float = 0.0,
    ds: float = 1.0,
    return_theta: bool = False,
    return_seg_id: bool = False,
):
   
    # initial state
    theta = math.radians(theta0_deg)   # radians, positive up
    d_vals = [0.0]
    z_vals = [z0]
    theta_vals = [theta]
    segid_vals = [0]

    d_cum = 0.0

    for idx, seg in enumerate(segments, start=1):
        t = seg["type"].lower()

        if t == "tangent":
            L = float(seg["length"])
            n = int(math.ceil(L / ds))
            ds_eff = L / n

            # generate n steps; skip the first shared joint when appending
            for _ in range(n):
                d_cum += ds_eff
                z_next = z_vals[-1] + ds_eff * math.sin(theta)  # dz/ds = sin(theta)
                d_vals.append(d_cum)
                z_vals.append(z_next)
                theta_vals.append(theta)
                segid_vals.append(idx)

        elif t == "arc":
            R = float(seg["R"])
            delta_deg = float(seg["delta_deg"])
            delta_rad = math.radians(delta_deg)
            L = abs(R * delta_rad)
            n = int(math.ceil(L / ds))
            ds_eff = L / n
            # signed curvature kappa = dθ/ds
            kappa = (delta_rad / L) if L != 0.0 else 0.0

            for _ in range(n):
                theta_next = theta + kappa * ds_eff
                theta_mid = 0.5 * (theta + theta_next)          # midpoint for z integration
                d_cum += ds_eff
                z_next = z_vals[-1] + ds_eff * math.sin(theta_mid)

                d_vals.append(d_cum)
                z_vals.append(z_next)
                theta_vals.append(theta_next)
                segid_vals.append(idx)

                theta = theta_next  # advance state

        else:
            pass

    d = np.asarray(d_vals, dtype=float)
    z = np.asarray(z_vals, dtype=float)

    if return_theta or return_seg_id:
        theta_arr = np.asarray(theta_vals, dtype=float)
        if return_seg_id:
            seg_id = np.asarray(segid_vals, dtype=int)
            return d, z, theta_arr, seg_id
        return d, z, theta_arr

    return d, z
