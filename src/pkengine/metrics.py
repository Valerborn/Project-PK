# src/pkengine/metrics.py
import numpy as np
from typing import Tuple

def cmax(C: np.ndarray) -> float:
    """Global maximum concentration (mg/L)."""
    return float(np.max(C))

def tmax(t: np.ndarray, C: np.ndarray) -> float:
    """Time of maximum concentration (h)."""
    return float(t[int(np.argmax(C))])

def cmin(C: np.ndarray) -> float:
    """Global minimum concentration (mg/L)."""
    return float(np.min(C))

def tmin(t: np.ndarray, C: np.ndarray) -> float:
    """Time of minimum concentration (h)."""
    return float(t[int(np.argmin(C))])

def cmax_tmax(t: np.ndarray, C: np.ndarray) -> Tuple[float, float]:
    """Return Cmax (mg/L) and Tmax (h)."""
    idx = np.argmax(C)
    return float(C[idx]), float(t[idx])

def auc_trapz(t: np.ndarray, C: np.ndarray) -> float:
    """Area Under the Curve (AUC) via trapezoidal rule (mg*h/L)."""
    return float(np.trapezoid(C, t))

def ctrough(t: np.ndarray, C: np.ndarray, interval_h: float) -> float:
    """
    Return trough concentration (concentration at the end of each dosing interval).
    interval_h: dosing interval in hours (e.g., 72 for q3d).
    """
    n_intervals = int(t[-1] // interval_h)
    troughs = []
    for i in range(1, n_intervals + 1):
        target = i * interval_h
        idx = np.argmin(np.abs(t - target))
        troughs.append(C[idx])
    return float(troughs[-1]) if troughs else float("nan")

def cavg(t: np.ndarray, C: np.ndarray) -> float:
    """Average concentration over simulation horizon."""
    return float(np.mean(C))

def _window_indices_for_last_interval(t: np.ndarray, interval_h: float) -> np.ndarray:
    """
    Return a boolean mask for samples in the last full dosing interval.
    If no full interval fits, fall back to all samples.
    """
    if interval_h <= 0:
        return np.ones_like(t, dtype=bool)
    last_edge = (t[-1] // interval_h) * interval_h
    start = last_edge - interval_h
    if start < t[0]:
        return np.ones_like(t, dtype=bool)
    return (t >= start) & (t <= last_edge)

def peak_to_trough_ratio(t: np.ndarray, C: np.ndarray, interval_h: float | None = None) -> float:
    """
    Peak-to-Trough Ratio (PTR) = Cmax / Cmin.
    If interval_h is provided, compute over the last full dosing interval; otherwise use the full series.
    """
    if interval_h:
        mask = _window_indices_for_last_interval(t, float(interval_h))
        Cw = C[mask]
    else:
        Cw = C
    cmax_val = float(np.max(Cw))
    cmin_val = float(np.min(Cw))
    if cmin_val <= 0:
        return float('inf')
    return cmax_val / cmin_val

def fluctuation_index(t: np.ndarray, C: np.ndarray, interval_h: float | None = None) -> float:
    """
    Fluctuation Index (FI) = (Cmax - Cmin) / Cavg.
    If interval_h is provided, compute over the last full dosing interval; otherwise use the full series.
    """
    if interval_h:
        mask = _window_indices_for_last_interval(t, float(interval_h))
        Cw = C[mask]
    else:
        Cw = C
    cmax_val = float(np.max(Cw))
    cmin_val = float(np.min(Cw))
    cavg_val = float(np.mean(Cw))
    if cavg_val == 0.0:
        return float('inf')
    return (cmax_val - cmin_val) / cavg_val

def _interp_shift(C: np.ndarray, t: np.ndarray, shift_h: float) -> np.ndarray:
    """
    Interpolate C(t - shift_h) onto the same time grid as t.
    For t < shift_h, returns extrapolated values at the start.
    """
    if shift_h == 0:
        return C.copy()
    t_shift = t - shift_h
    # Use numpy interpolation; left/right fill with edge values
    return np.interp(t_shift, t, C, left=C[0], right=C[-1])

def steady_state_window_mask(t: np.ndarray, C: np.ndarray, interval_h: float,
                             tol: float = 0.05, max_lookback: int = 10) -> np.ndarray:
    """
    Find the last full dosing interval where the profile is approximately at steady state.

    A window [T-interval_h, T] is considered steady if the relative difference
    between C(t) and C(t - interval_h) across that window is <= tol (default 5%).

    Parameters
    ----------
    t : array of time points (hours)
    C : array of concentrations (mg/L)
    interval_h : dosing interval (hours)
    tol : relative tolerance for steady-state comparison (0.05 = 5%)
    max_lookback : how many intervals to look back if the last one isn't steady

    Returns
    -------
    mask : boolean array selecting the chosen steady-state interval.
           If no interval satisfies the tolerance, returns the mask for the last interval.
    """
    if interval_h <= 0:
        # No interval concept; use entire series
        return np.ones_like(t, dtype=bool)

    # Compute pointwise relative difference vs. curve one interval earlier
    C_shift = _interp_shift(C, t, shift_h=interval_h)
    denom = np.maximum(np.abs(C), 1e-12)  # avoid divide-by-zero
    rel_diff = np.abs(C - C_shift) / denom

    # Try the last interval first, then step back
    for k in range(max_lookback):
        end_edge = (t[-1] // interval_h) * interval_h - (k * interval_h)
        start_edge = end_edge - interval_h
        if start_edge < t[0]:
            break
        mask = (t >= start_edge) & (t <= end_edge)
        if np.all(rel_diff[mask] <= tol):
            return mask

    # Fallback: last full interval
    last_edge = (t[-1] // interval_h) * interval_h
    start = last_edge - interval_h
    if start < t[0]:
        return np.ones_like(t, dtype=bool)
    return (t >= start) & (t <= last_edge)

def peak_to_trough_ratio_ss(t: np.ndarray, C: np.ndarray, interval_h: float,
                            tol: float = 0.05) -> float:
    """
    Peak-to-Trough Ratio at steady state (PTR_ss).
    Uses the steady_state_window_mask to select the interval.
    """
    mask = steady_state_window_mask(t, C, interval_h, tol=tol)
    Cw = C[mask]
    cmax_val = float(np.max(Cw))
    cmin_val = float(np.min(Cw))
    if cmin_val <= 0:
        return float('inf')
    return cmax_val / cmin_val

def fluctuation_index_ss(t: np.ndarray, C: np.ndarray, interval_h: float,
                         tol: float = 0.05) -> float:
    """
    Fluctuation Index at steady state (FI_ss) = (Cmax - Cmin) / Cavg over the steady interval.
    """
    mask = steady_state_window_mask(t, C, interval_h, tol=tol)
    Cw = C[mask]
    cmax_val = float(np.max(Cw))
    cmin_val = float(np.min(Cw))
    cavg_val = float(np.mean(Cw))
    if cavg_val == 0.0:
        return float('inf')
    return (cmax_val - cmin_val) / cavg_val