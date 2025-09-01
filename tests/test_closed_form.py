

import math
import numpy as np

from pkengine.types import Regimen, Dose, Drug
from pkengine.dosing import fixed_every_n_days, single_dose, combine_regimens
from pkengine.solvers import simulate_1c_first_order, simulate_multi
from pkengine.metrics import (
    cmax, tmax, cmin, tmin, auc_trapz,
    peak_to_trough_ratio, fluctuation_index,
    peak_to_trough_ratio_ss, fluctuation_index_ss,
)


def test_iv_bolus_exponential_decay():
    """
    For a 1-compartment model with linear elimination, an IV bolus follows:
      C(t) = C0 * exp(-k t),  where k = CL / V and C0 = dose / V.
    We simulate an IV bolus and compare to the analytical solution.
    """
    CL, V, ka = 12.0, 300.0, 0.10  # ka is irrelevant for pure IV bolus
    dose_mg = 3000.0               # C0 should be ~10 mg/L
    t_end = 48.0

    # Build a regimen with a single IV bolus at t=0
    reg = Regimen(doses=(Dose(drug_id="test", route="iv_bolus", amount_mg=dose_mg, start_h=0.0),))

    # Run solver
    t, C = simulate_1c_first_order(CL_L_per_h=CL, V_L=V, ka_per_h=ka, regimen=reg, t_end_h=t_end, dt_h=0.5)

    # Expected curve
    k = CL / V
    C0 = dose_mg / V
    C_expected = C0 * np.exp(-k * t)

    # Allow small numerical tolerances
    assert np.isclose(C[0], C0, rtol=0.02, atol=0.05)
    assert np.allclose(C, C_expected, rtol=0.05, atol=0.1)


def test_single_drug_im_q3d_smoke():
    """
    Smoke test: 250 mg IM every 3 days for 8 weeks should produce a
    non-negative concentration-time profile and reasonable metrics.
    """
    CL, V, ka = 12.0, 300.0, 0.10
    interval_h = 72.0
    weeks = 8
    t_end = weeks * 7 * 24

    reg = fixed_every_n_days(amount_mg=250, every_days=3, weeks=weeks, route="im", drug_id="testosterone")

    t, C = simulate_1c_first_order(CL_L_per_h=CL, V_L=V, ka_per_h=ka, regimen=reg, t_end_h=t_end, dt_h=1.0)

    assert len(t) == len(C) and len(t) > 10
    assert np.all(C >= -1e-9)  # allow tiny numerical noise
    assert float(cmax(C)) > 0.0
    assert 0.0 <= float(tmax(t, C)) <= t_end
    assert float(auc_trapz(t, C)) > 0.0

    # Interval-aware metrics
    ptr_last = peak_to_trough_ratio(t, C, interval_h=interval_h)
    fi_last = fluctuation_index(t, C, interval_h=interval_h)
    ptr_ss = peak_to_trough_ratio_ss(t, C, interval_h=interval_h, tol=0.05)
    fi_ss = fluctuation_index_ss(t, C, interval_h=interval_h, tol=0.05)

    assert np.isfinite(ptr_last) and ptr_last > 1.0
    assert np.isfinite(fi_last) and fi_last >= 0.0
    assert np.isfinite(ptr_ss) and ptr_ss > 1.0
    assert np.isfinite(fi_ss) and fi_ss >= 0.0


def test_multi_drug_two_regimens():
    """
    Multi-drug scenario: testosterone + nandrolone with different PK params.
    The curves should differ and both be present in the results dict.
    """
    # Build regimens
    reg_testo = fixed_every_n_days(250, every_days=3, weeks=8, route="im", drug_id="testosterone")
    reg_nandro = fixed_every_n_days(200, every_days=7, weeks=8, route="im", drug_id="nandrolone")
    combined = combine_regimens(reg_testo, reg_nandro)

    # Map of PK params per drug
    params = {
        "testosterone": {"CL": 12.0, "V": 300.0, "ka": 0.10},
        "nandrolone":   {"CL": 10.0, "V": 250.0, "ka": 0.07},
    }
    t_end = 8 * 7 * 24

    results = simulate_multi(params, combined, t_end_h=t_end, dt_h=1.0)

    assert set(results.keys()) == {"testosterone", "nandrolone"}
    t_te, C_te = results["testosterone"]
    t_na, C_na = results["nandrolone"]

    # Basic shape checks
    assert len(t_te) == len(C_te) and len(t_te) > 10
    assert len(t_na) == len(C_na) and len(t_na) > 10
    assert np.all(C_te >= -1e-9)
    assert np.all(C_na >= -1e-9)

    # They shouldn't be identical because PK params differ
    assert not np.allclose(C_te, C_na)