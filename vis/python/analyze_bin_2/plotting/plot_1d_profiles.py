import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from scipy.signal import find_peaks

def plot_rc(data_dict, user_params):
    """
    Plots concatenated 1D radial cut (RC) profile from AthenaK blocks along x=user_params['profile_slice'].
    Uses output from extract_athenak_slice_as_dataframe.
    """
    df_quantities = data_dict['df_quantities']
    df_extents = data_dict['df_extents']
    block_nx2, block_nx1 = data_dict['block_shape']

    data_1d_all = []
    y_coords_all = []

    for block in range(df_extents.shape[0]):
        extent = df_extents.iloc[block].values  # [x_min, x_max, y_min, y_max]
        block_data_2d = df_quantities.loc[block].unstack(level=-1).values
        x_min, x_max, y_min, y_max = extent
        ny, nx = block_data_2d.shape
        x_coords = np.linspace(x_min, x_max, nx)
        y_coords = np.linspace(y_min, y_max, ny)
        # Find column index closest to specified profile slice
        slice_x = user_params['profile_slice']
        if (x_min <= slice_x <= x_max) or (x_min <= 0.0 <= x_max):
            ix = np.argmin(np.abs(x_coords - slice_x))
            data_1d = block_data_2d[:, ix]
            data_1d_all.extend(data_1d)
            y_coords_all.extend(y_coords)

    y_values = np.array(y_coords_all)
    data_values = np.array(data_1d_all)

    max_idx = np.argmax(data_values)
    y_peak = y_values[max_idx]
    data_peak = data_values[max_idx]

    plt.figure(figsize=(10, 6))
    plt.plot(y_values, data_values, 'o-', markersize=3)
    plt.axvline(x=y_peak, color='r', linestyle='--', label=f'Peak at y={y_peak * user_params.get("axes_scale", 1.0):.3f} kpc')
    plt.scatter([y_peak], [data_peak], color='red')
    plt.xlabel('y')
    plt.ylabel(f"Data at x = {user_params['profile_slice']}")
    plt.title(f'Concatenated 1D slice along x={user_params["profile_slice"]} across all valid blocks')
    plt.legend()
    plt.grid(alpha=0.3)
    out_dir = Path(user_params['output_path']) / "rc_profile"
    os.makedirs(out_dir, exist_ok=True)
    if user_params.get('loop_bin', False):
        fname = os.path.join(out_dir, f"rc_{user_params['slice_number']}.png")
    else:
        fname = os.path.join(out_dir, 'rc.png')
    plt.savefig(fname, bbox_inches="tight")
    plt.close()
    print(f"Profile saved at {fname}")

def plot_zh(data_dict, user_params):
    """
    Plots concatenated 1D vertical cut (ZH) profile from AthenaK blocks along y=user_params['profile_slice'].
    Uses output from extract_athenak_slice_as_dataframe.
    """
    df_quantities = data_dict['df_quantities']
    df_extents = data_dict['df_extents']
    block_nx2, block_nx1 = data_dict['block_shape']

    data_1d_all = []
    x_coords_all = []

    for block in range(df_extents.shape[0]):
        extent = df_extents.iloc[block].values  # [x_min, x_max, y_min, y_max]
        block_data_2d = df_quantities.loc[block].unstack(level=-1).values
        x_min, x_max, y_min, y_max = extent
        ny, nx = block_data_2d.shape
        x_coords = np.linspace(x_min, x_max, nx)
        y_coords = np.linspace(y_min, y_max, ny)
        slice_y = user_params['profile_slice']
        # Only proceed if y=slice_y is within block
        if (y_min <= slice_y <= y_max):
            iy = np.argmin(np.abs(y_coords - slice_y))
            data_1d = block_data_2d[iy, :]
            data_1d_all.extend(data_1d)
            x_coords_all.extend(x_coords)

    x_values = np.array(x_coords_all)
    data_values = np.array(data_1d_all)
    # Mask singularity at x=0
    mask = x_values != 0
    x_filtered = x_values[mask]
    data_filtered = data_values[mask]

    peaks, _ = find_peaks(data_filtered)
    if len(peaks) > 0:
        last_peak_idx = peaks[-1]
        x_peak = x_filtered[last_peak_idx]
        data_peak = data_filtered[last_peak_idx]
        label_str = f'Last peak at x = {x_peak * user_params.get("axes_scale", 1.0):.3f} kpc'
    else:
        last_peak_idx = np.argmax(data_filtered)
        x_peak = x_filtered[last_peak_idx]
        data_peak = data_filtered[last_peak_idx]
        label_str = f'Global max at x = {x_peak * user_params.get("axes_scale", 1.0):.3f} kpc'

    plt.figure(figsize=(10,6))
    plt.plot(x_values, data_values, 'o-', markersize=3)
    plt.axvline(x=x_peak, color='r', linestyle='--', label=label_str)
    plt.scatter([x_peak], [data_peak], color='red')
    plt.xlabel('x')
    plt.ylabel(f"Data at y={user_params['profile_slice']}")
    plt.legend()
    plt.grid(alpha=0.3)
    out_dir = Path(user_params['output_path']) / "zh_profile"
    os.makedirs(out_dir, exist_ok=True)
    if user_params.get('loop_bin', False):
        fname = os.path.join(out_dir, f"zh_{user_params['slice_number']}.png")
    else:
        fname = os.path.join(out_dir, 'zh.png')
    plt.savefig(fname, bbox_inches="tight")
    plt.close()
    print(f"Profile saved at {fname}")
