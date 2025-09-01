# src/pkengine/models/one_compartment.py
import numpy as np

def one_compartment_first_order(t, y, ka, CL, V, doses):
    """
    One-compartment model with first-order absorption and elimination.
    Two states:
      y[0] = drug in absorption depot (mg)
      y[1] = drug in central compartment (mg)

    Parameters:
      t      : current time (h)
      y      : current state vector [A_depot, A_central]
      ka     : absorption rate constant (1/h)
      CL     : clearance (L/h)
      V      : volume of distribution (L)
      doses  : list of Dose events (with drug_id matched before calling)
    """
    A_gut, A_c = y

    # Base ODEs
    dA_gut_dt = -ka * A_gut
    dA_c_dt   = ka * A_gut - (CL / V) * A_c

    # Continuous inputs: IV infusion as zero-order rate during [start, start+duration]
    for d in doses:
        if d.route == "iv_infusion" and d.start_h <= t <= d.start_h + d.duration_h:
            rate = d.amount_mg / d.duration_h
            dA_c_dt += rate

    return [dA_gut_dt, dA_c_dt]