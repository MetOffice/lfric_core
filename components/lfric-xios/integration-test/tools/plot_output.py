#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Plotting function for LFRic-XIOS integration tests. Depends on xarray and
matplotlib packages which can be loaded as part of scitools for users in
the Met Office.
"""

from pathlib import Path
import matplotlib
import xarray as xr
matplotlib.use('Agg')
import matplotlib.pyplot as plt # pylint: disable=wrong-import-position, ungrouped-imports

def _get_ts_data(file_path, field_id):
    """
    Get time series data for a given field from a netCDF file.
    """

    ds = xr.open_dataset(file_path, engine='netcdf4', decode_timedelta=False)
    ts = ds[field_id].mean(ds[field_id].dims[1::])
    time = ds[field_id].coords['time']

    return ts, time


def plot_test_output(in_file: Path, out_file: Path, varname: str, plot_file_path: Path):
    """
    Visually compare input and output data.
    """

    input_ts, input_time = _get_ts_data(in_file, varname)
    output_ts, output_time = _get_ts_data(out_file, varname)

    plt.rcParams["font.family"] = "serif"
    _, ax = plt.subplots(figsize=([10.8, 4.8]))
    ax.scatter(output_time, output_ts, c='C0', s=50)
    ax.plot(output_time, output_ts, linestyle='--', lw=2, label="Model output data")
    ax.scatter(input_time, input_ts, c='C3', marker='s', s=100, label="Input data")

    ax.set_xlabel("Date/Time")
    ax.set_ylabel("Mean model data")

    plt.legend(frameon=False)
    plt.savefig(f"{plot_file_path}", bbox_inches="tight")
    plt.close()
