# src/pkengine/dosing.py
from __future__ import annotations

from typing import Sequence, Tuple
import numpy as np

from .types import Dose, Regimen, Route


def single_dose(amount_mg: float, start_h: float, route: Route = "im", duration_h: float = 0.0, *, drug_id: str = "default") -> Regimen:
    """
    Create a regimen with exactly one dose.
    Examples:
      - 250 mg IM at t=0 h
      - 500 mg IV infusion starting at t=4 h over 2 hours (duration_h=2)
    drug_id       : identifier/name of the compound this dose belongs to (default "default")
    """
    _validate_positive("amount_mg", amount_mg)
    _validate_non_negative("start_h", start_h)
    if route == "iv_infusion":
        _validate_positive("duration_h", duration_h)
    else:
        _validate_zero_or_positive("duration_h", duration_h)
        if duration_h > 0:
            raise ValueError("duration_h should be 0 unless route is 'iv_infusion'.")

    return Regimen(doses=(Dose(drug_id=drug_id, route=route, amount_mg=float(amount_mg),
                               start_h=float(start_h), duration_h=float(duration_h)),))


def fixed_every_n_days(amount_mg: float, every_days: int, weeks: int,
                       route: Route = "im", start_offset_h: float = 0.0, *, drug_id: str = "default") -> Regimen:
    """
    Make a repeated schedule like: 250 mg IM every 3 days for 8 weeks.

    amount_mg      : size of each dose, mg
    every_days     : spacing between doses, in whole days (e.g., 3)
    weeks          : total regimen length in weeks
    route          : im/sc/po (or iv_bolus if you really want repeated IV boluses)
    start_offset_h : shift the very first dose by some hours (e.g., 9:00 = 9.0)
    drug_id        : identifier/name of the compound this dose belongs to (default "default")
    """
    _validate_positive("amount_mg", amount_mg)
    _validate_positive_int("every_days", every_days)
    _validate_positive_int("weeks", weeks)
    _validate_non_negative("start_offset_h", start_offset_h)

    total_days = weeks * 7
    # Dose times: 0d, every_days, 2*every_days, ... <= total_days  (then convert to hours)
    day_starts = np.arange(0, total_days + 1, every_days, dtype=float)  # e.g., [0, 3, 6, 9, ...]
    times_h = day_starts * 24.0 + float(start_offset_h)

    doses = tuple(
        Dose(drug_id=drug_id, route=route, amount_mg=float(amount_mg), start_h=float(t), duration_h=0.0)
        for t in times_h
    )
    return Regimen(doses=doses)


def iv_infusion(amount_mg: float, start_h: float, duration_h: float, *, drug_id: str = "default") -> Regimen:
    """
    A single IV infusion (zero-order input) of a given total amount over a duration.
    Example: 1000 mg starting at t=0 h over 2 h.
    drug_id       : identifier/name of the compound this dose belongs to (default "default")
    """
    _validate_positive("amount_mg", amount_mg)
    _validate_non_negative("start_h", start_h)
    _validate_positive("duration_h", duration_h)

    return Regimen(doses=(Dose(drug_id=drug_id, route="iv_infusion",
                               amount_mg=float(amount_mg),
                               start_h=float(start_h),
                               duration_h=float(duration_h)),))


def combine_regimens(*regimens: Regimen) -> Regimen:
    """
    Merge multiple regimens into one (e.g., IM schedule + an IV loading dose).
    Doses are simply concatenated; the simulator can sort them by time.
    Supports multi-drug regimens by merging doses regardless of drug_id.
    """
    all_doses: list[Dose] = []
    for r in regimens:
        all_doses.extend(r.doses)
    # Optional: sort now for human readability
    all_doses_sorted = tuple(sorted(all_doses, key=lambda d: (d.start_h, d.route)))
    return Regimen(doses=all_doses_sorted)


def from_explicit_schedule(entries: Sequence[Tuple[float, float]], route: Route = "im", *, drug_id: str = "default") -> Regimen:
    """
    Build a regimen from manual (time_h, amount_mg) entries.
    Example: entries=[(0.0, 250), (72.0, 250), (144.0, 250)]
    drug_id       : identifier/name of the compound this dose belongs to (default "default")
    """
    doses: list[Dose] = []
    for start_h, amount_mg in entries:
        _validate_positive("amount_mg", amount_mg)
        _validate_non_negative("start_h", start_h)
        doses.append(Dose(drug_id=drug_id, route=route, amount_mg=float(amount_mg), start_h=float(start_h)))
    doses.sort(key=lambda d: d.start_h)
    return Regimen(doses=tuple(doses))


# --------------------------
# Small input validators
# --------------------------
def _validate_positive(name: str, x: float) -> None:
    if not (x > 0):
        raise ValueError(f"{name} must be > 0 (got {x}).")

def _validate_zero_or_positive(name: str, x: float) -> None:
    if not (x >= 0):
        raise ValueError(f"{name} must be >= 0 (got {x}).")

def _validate_non_negative(name: str, x: float) -> None:
    if x < 0:
        raise ValueError(f"{name} must be >= 0 (got {x}).")

def _validate_positive_int(name: str, x: int) -> None:
    if not (isinstance(x, int) and x > 0):
        raise ValueError(f"{name} must be a positive integer (got {x}).")