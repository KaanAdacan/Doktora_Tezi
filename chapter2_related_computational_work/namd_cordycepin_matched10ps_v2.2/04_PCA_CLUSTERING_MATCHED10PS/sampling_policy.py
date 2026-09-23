#!/usr/bin/env python3
from __future__ import annotations
import numpy as np

TARGET_SAMPLING = {
    "MAPK": {"native_interval_ns": 0.010, "stride": 1, "analysis_interval_ns": 0.010},
    "ADK":  {"native_interval_ns": 0.002, "stride": 5, "analysis_interval_ns": 0.010},
}
EXPECTED_ANALYSIS_FRAMES = 10000
EXPECTED_START_NS = 0.01
EXPECTED_END_NS = 100.0


def matched_indices(times, target):
    """Return deterministic indices on the common 10-ps grid.

    MAPK native: 10 ps -> stride 1, indices 0..9999.
    ADK native :  2 ps -> stride 5, indices 4,9,14,...,49999.

    This is deterministic temporal subsampling only. It does not interpolate,
    smooth, filter, impute, remove outliers, or modify the raw trajectory.
    """
    times=np.asarray(times,dtype=float)
    cfg=TARGET_SAMPLING[target]
    stride=int(cfg["stride"])
    first=stride-1
    idx=np.arange(first,len(times),stride,dtype=int)
    selected=times[idx]
    expected=np.arange(1,EXPECTED_ANALYSIS_FRAMES+1,dtype=float)*0.01
    if len(selected)!=EXPECTED_ANALYSIS_FRAMES:
        raise ValueError(f"{target}: expected {EXPECTED_ANALYSIS_FRAMES} frames, got {len(selected)}")
    if not np.allclose(selected,expected,rtol=0,atol=1e-8):
        raise ValueError(f"{target}: selected times are not the exact 10-ps grid")
    return idx, selected
