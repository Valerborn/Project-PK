# src/pkengine/solvers.py
import numpy as np
from scipy.integrate import solve_ivp

from .types import Regimen, Dose
from .models.one_compartment import one_compartment_first_order
from .helpers import split_regimen_by_drug


def simulate_1c_first_order(CL_L_per_h: float, V_L: float, ka_per_h: float,
                            regimen: Regimen, t_end_h: float, dt_h: float = 1.0):
    """
    Simulate a one-compartment model with first-order absorption.

    Handles dosing events (IM/SC/PO depot entries and IV bolus) as instantaneous
    state jumps at their scheduled times and handles IV infusions as continuous
    zero-order inputs inside the ODE right-hand side.

    Returns:
      t : array of time points (hours)
      C : array of concentrations (mg/L)
    """
    # Initial conditions: nothing in gut, nothing in central
    y0 = [0.0, 0.0]

    # Time grid for outputs (start after t=0 to avoid pre-dose sampling)
    t_grid = np.arange(dt_h, t_end_h + dt_h, dt_h)

    # Build segment boundaries: dose starts and infusion ends to capture RHS changes
    boundaries: list[float] = [0.0]
    for d in regimen.doses:
        if 0.0 <= d.start_h <= t_end_h:
            boundaries.append(float(d.start_h))
        if d.route == "iv_infusion":
            end_t = d.start_h + d.duration_h
            if 0.0 <= end_t <= t_end_h:
                boundaries.append(float(end_t))
    boundaries.append(float(t_end_h))
    # Unique, sorted
    boundaries = sorted(set(boundaries))

    def rhs(t, y):
        return one_compartment_first_order(t, y, ka_per_h, CL_L_per_h, V_L, regimen.doses)

    # Apply any doses at t=0 as instantaneous jumps before integrating
    for d in regimen.doses:
        if np.isclose(d.start_h, 0.0):
            if d.route in ("im", "sc", "po"):
                y0[0] += float(d.amount_mg)
            elif d.route == "iv_bolus":
                y0[1] += float(d.amount_mg)

    t_out: list[float] = []
    Ac_out: list[float] = []

    prev = boundaries[0]
    for idx in range(1, len(boundaries)):
        curr = boundaries[idx]
        if curr <= prev:
            continue

        # Evaluation points for this segment (exclude prev to avoid duplicate except for first segment)
        if len(t_out) == 0:
            t_eval_seg = t_grid[(t_grid >= prev) & (t_grid <= curr)]
        else:
            t_eval_seg = t_grid[(t_grid > prev) & (t_grid <= curr)]

        if t_eval_seg.size == 0:
            # Integrate without sampling if no eval points fall in this segment
            sol_seg = solve_ivp(rhs, t_span=(prev, curr), y0=y0, method="RK45")
            y_end = sol_seg.y[:, -1]
        else:
            sol_seg = solve_ivp(rhs, t_span=(prev, curr), y0=y0, method="RK45", t_eval=t_eval_seg)
            t_out.extend(sol_seg.t.tolist())
            Ac_out.extend(sol_seg.y[1].tolist())
            y_end = sol_seg.y[:, -1]

        # Update state for instantaneous doses at exactly curr
        y0 = [float(y_end[0]), float(y_end[1])]
        for d in regimen.doses:
            if np.isclose(d.start_h, curr):
                if d.route in ("im", "sc", "po"):
                    y0[0] += float(d.amount_mg)
                elif d.route == "iv_bolus":
                    y0[1] += float(d.amount_mg)

        prev = curr

    # Convert to numpy arrays
    t_arr = np.asarray(t_out, dtype=float)
    Ac_arr = np.asarray(Ac_out, dtype=float)
    C = Ac_arr / V_L
    C = np.maximum(C, 0.0)
    return t_arr, C


def simulate_multi(drug_params: dict[str, dict], regimen: Regimen,
                   t_end_h: float, dt_h: float = 1.0) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """
    Simulate a combined regimen that may contain multiple drug_ids.

    Each drug_id is simulated independently using its own PK parameters.

    Parameters
    ----------
    drug_params : dict[str, dict]
        Mapping from drug_id to a dict with keys:
           - "CL": clearance in L/h (float)
           - "V" : volume of distribution in L (float)
           - "ka": absorption rate constant in 1/h (float)
    regimen : Regimen
        Regimen that may contain doses for multiple drugs (distinguished by Dose.drug_id).
    t_end_h : float
        Total simulation horizon in hours.
    dt_h : float, default 1.0
        Output sampling interval in hours.

    Returns
    -------
    dict[str, tuple[np.ndarray, np.ndarray]]
        A mapping from drug_id to (t, C) arrays, where t is time in hours
        and C is concentration in mg/L.
    """
    # 1) Group the combined regimen by drug_id
    per_drug = split_regimen_by_drug(regimen)

    results: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for drug_id, sub_reg in per_drug.items():
        params = drug_params.get(drug_id)
        if params is None:
            raise KeyError(f"Missing PK params for drug_id '{drug_id}' in drug_params.")

        # 2) Run the single-drug simulator for this set of doses
        t, C = simulate_1c_first_order(
            CL_L_per_h=float(params['CL']),
            V_L=float(params['V']),
            ka_per_h=float(params['ka']),
            regimen=sub_reg,
            t_end_h=t_end_h,
            dt_h=dt_h,
        )
        results[drug_id] = (t, C)

    return results