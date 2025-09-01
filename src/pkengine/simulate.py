# src/pkengine/simulate.py
from .types import Drug, Regimen
from .solvers import simulate_1c_first_order, simulate_multi

def run_single(drug: Drug, regimen: Regimen, t_end_h: float, dt_h: float = 1.0):
    """
    High-level wrapper for simulating a single drug with its parameters.
    """
    return simulate_1c_first_order(
        CL_L_per_h=drug.CL_L_per_h,
        V_L=drug.V_L,
        ka_per_h=drug.ka_per_h,
        regimen=regimen,
        t_end_h=t_end_h,
        dt_h=dt_h,
    )

def run_multi(drug_map: dict[str, Drug], regimen: Regimen,
              t_end_h: float, dt_h: float = 1.0):
    """
    High-level wrapper for simulating a regimen with multiple drugs.
    drug_map: drug_id -> Drug object
    """
    # Convert Drug objects into the dict format simulate_multi expects
    params = {
        did: {"CL": d.CL_L_per_h, "V": d.V_L, "ka": d.ka_per_h}
        for did, d in drug_map.items()
    }
    return simulate_multi(params, regimen, t_end_h, dt_h)