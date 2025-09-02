# src/pkengine/types.py
from dataclasses import dataclass
from typing import Literal, Sequence

# We keep *all* time in HOURS internally. (Easy math, avoids unit drift.)
Route = Literal["iv_bolus", "iv_infusion", "im", "sc", "po"]

@dataclass(frozen=True)
class Dose:
    """
    A single administration of a medication.

    drug_id     : identifier for the compound this dose belongs to (e.g., "testosterone")
    route       : how the dose is given (iv_bolus, iv_infusion, im, sc, po)
    amount_mg   : dose size in milligrams
    start_h     : when the dose starts (in hours from time 0)
    duration_h  : how long it runs (0 for bolus or non-infusion routes)
                  e.g., a 2-hour IV infusion has duration_h=2.0
    """
    drug_id: str
    route: Route
    amount_mg: float
    start_h: float
    duration_h: float = 0.0  # >0 only for infusions / zero-order inputs


@dataclass(frozen=True)
class Regimen:
    """
    A collection of Dose objects that defines the full schedule.

    doses : a sequence (list/tuple) of Dose entries. Order doesn't matter;
            the simulator can sort them by start_h when it runs.
    """
    doses: Sequence[Dose]

@dataclass(frozen=True)
class Drug:
    """
    Basic PK parameters for a one-compartment, first-order model.
    """
    CL_L_per_h: float
    V_L: float
    ka_per_h: float
    drug_id: str = "default"