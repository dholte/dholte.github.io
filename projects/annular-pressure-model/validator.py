class InputValidationError(ValueError):
    pass

def validate_inputs(geom: Geometry,
                    op: Operation,
                    fluid=None,  # FluidNewtonian or FluidBingham or None
                    model: str = "newtonian"):
    errors = []

    # Geometry
    if not (geom.Db > 0): errors.append("Db must be > 0.")
    if not (geom.Dp > 0): errors.append("Dp must be > 0.")
    if geom.Db <= geom.Dp: errors.append("Db must be greater than Dp (annulus must exist).")

    # Operation
    if not (op.Q > 0): errors.append("Q (flow rate) must be > 0.")
    if not (op.v_min > 0): errors.append("v_min must be > 0.")

    # Fluid
    model = model.lower()
    if model == "newtonian":
        if not isinstance(fluid, FluidNewtonian):
            errors.append("For model='newtonian', provide FluidNewtonian(rho, mu).")
        else:
            if not (fluid.rho > 0): errors.append("rho must be > 0.")
            if not (fluid.mu > 0):  errors.append("mu must be > 0.")
    elif model == "bingham":
        if not isinstance(fluid, FluidBingham):
            errors.append("For model='bingham', provide FluidBingham(mu_p, tau_y).")
        else:
            if not (fluid.mu_p > 0):  errors.append("mu_p must be > 0.")
            if not (fluid.tau_y >= 0): errors.append("tau_y must be ≥ 0.")
    else:
        errors.append("model must be either 'newtonian' or 'bingham'.")

    if errors:
        raise InputValidationError("Input validation failed:\n- " + "\n- ".join(errors))
