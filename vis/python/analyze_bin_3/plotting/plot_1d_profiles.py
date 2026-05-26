import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from scipy.signal import find_peaks
from data_processing.profile_data import get_zh_data

def plot_rc(data_dict, user_params):
    """
    Plots concatenated 1D radial cut (RC) profile from AthenaK blocks along x=user_params['profile_slice'].
    Uses output from extract_athenak_slice_as_dataframe.
    """
    variable_name=user_params['variable']
    if variable_name.startswith("derived:"):
        variable_name=variable_name.split(":", 1)[1].strip()
        user_params.update(variable=variable_name)
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
    x_filtered, data_filtered = get_zh_data(data_dict, user_params)

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
    import os
    from pathlib import Path
    import matplotlib.pyplot as plt

    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update({
        "font.size": 16,
        "font.weight": "bold",
        "axes.labelsize": 18,
        "axes.labelweight": "bold",
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
        "axes.linewidth": 1.2,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    })

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.scatter(
        x_filtered,
        data_filtered,
        color=user_params["color"],
        s=18,
        alpha=0.85,
        edgecolors="none",
        rasterized=True,
    )

    ax.set_xlabel(user_params["xlabel"], labelpad=8)
    ax.set_ylabel(user_params["cmap_label"], labelpad=8)
    ax.set_xlim(user_params["clim"][0], user_params["clim"][1])

    ax.tick_params(direction="in", length=5, width=1.1, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)

    ax.grid(alpha=0.25, linewidth=0.8)

    fig.tight_layout()

    out_dir = Path(user_params["output_path"]) / "zh_profile"
    os.makedirs(out_dir, exist_ok=True)
    if user_params.get("loop_bin", False):
        fname = os.path.join(out_dir, f"zh_{user_params['slice_number']}.png")
    else:
        fname = os.path.join(out_dir, "zh.png")

    fig.savefig(fname, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

    print(f"Profile saved at {fname}")

def plot_x_profile(data_dict, user_params):
    """
    Plots concatenated 1D profile along x-axis from (0,0,0) to (xmax,0,0) 
    from AthenaK blocks at y=0, z=0.
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
        
        # Fixed slice at y=0, z=0 (central line along x-axis)
        slice_y = 0.0
        
        # Only proceed if y=0 is within block
        if (y_min <= slice_y <= y_max):
            iy = np.argmin(np.abs(y_coords - slice_y))
            data_1d = block_data_2d[iy, :]
            data_1d_all.extend(data_1d)
            x_coords_all.extend(x_coords)

    x_values = np.array(x_coords_all)
    data_values = np.array(data_1d_all)
    
    # Only plot from x=0 to x=xmax (keep x=0, remove negative x)
    mask = x_values >= 0
    x_filtered = x_values[mask]
    data_filtered = data_values[mask]
    
    # Sort by x coordinate to ensure monotonic increasing order
    sort_idx = np.argsort(x_filtered)
    x_sorted = x_filtered[sort_idx]
    data_sorted = data_filtered[sort_idx]

    plt.figure(figsize=(10,6))
    if user_params['norm']=="log":
        data_sorted = np.log10(data_sorted)
        if user_params['xlabel'].endswith("[log]"):
            plt.plot(np.log10(user_params['axes_scale']*x_sorted), data_sorted, 'o-',color=user_params['color'], markersize=3)
        else:
            plt.plot(user_params['axes_scale']*x_sorted, data_sorted, 'o-',color=user_params['color'], markersize=3)
    else:
        plt.plot(user_params['axes_scale']*x_sorted, data_sorted, 'o-',color=user_params['color'], markersize=3)
    plt.ylim([
        user_params['clim'][0] or np.min(data_sorted),
        user_params['clim'][1] or np.max(data_sorted)
    ])
    plt.xlabel(user_params['xlabel'])
    plt.ylabel(user_params['cmap_label'])
    
    plt.grid(alpha=0.3)
    
    if user_params['variable']=="eint":
        out_dir = Path(user_params['output_path'])/ "x_profile" / "pres"
    elif user_params['variable'].startswith("derived:"):
        variable=user_params['variable'].split(":", 1)[1].strip()
        out_dir = Path(user_params['output_path'])/ "x_profile" / variable
    else:    
        out_dir = Path(user_params['output_path'])/ "x_profile" / user_params["variable"] 
    os.makedirs(out_dir, exist_ok=True)
    if user_params.get('loop_bin', False):
        # plt.title(f"t={user_params['slice_number']}"+r'$t_{ff}$')
        fname = os.path.join(out_dir, f"t_{user_params['slice_number']}.png")
    else:
        plt.title("x profile")
        fname = os.path.join(out_dir, 'x_profile.png')
    plt.savefig(fname, bbox_inches="tight")
    plt.close()
    print(f"X-profile saved at {fname}")

def plot_spherical_volume_weighted_avg_profile(profile, user_params):
    """
    Computes and plots 3D spherically-averaged profile from extract_athenak_3D_block output.
    Uses fast per-block NumPy arrays instead of slow DataFrame reconstruction.
    
    Parameters:
    -----------
    profile : dict
        Output from data_processing.block_data.compute_radial_profile with keys:
        r_centers, profile, std, count, sum
    
    user_params : dict
        Must include center_x, center_y, center_z, num_radial_bins, etc.
    """
    
    # ===== Create figure =====
    fig, ax = plt.subplots(figsize=(12, 7))
    
    axes_scale = user_params.get('axes_scale', 1.0)
    color = user_params.get('color', 'blue')
    
    # Apply scaling and normalization
    if user_params.get('norm') == "log":
        if user_params['xlabel'].endswith("[log]"):
            r_plot = np.log10(axes_scale * profile['r_centers'] + 1e-30)
        else:
            r_plot = axes_scale * profile['r_centers']
        
        #Log scale for y-axis
        radial_profile_plot = np.log10(profile['profile'] + 1e-30)
        # For log scale: if y = log10(Y), then dy ≈ dY / (Y * ln(10))
        radial_std_plot = profile['std'] / (profile['profile'] * np.log(10))

        
    else:
        radial_profile_plot = profile['profile']
        radial_std_plot = profile['std']
        r_plot = axes_scale * profile['r_centers']
    
    
    # Plot line with shaded uncertainty region
    ax.plot(r_plot, radial_profile_plot, 
            color=color, 
            linewidth=2.5,
            marker='o',
            markersize=6,
            label='Volume-weighted mean',
            zorder=3)
    
    # Add shaded region (±1σ confidence band)
    ax.fill_between(r_plot, 
                     radial_profile_plot - radial_std_plot,
                     radial_profile_plot + radial_std_plot,
                     color=color, 
                     alpha=0.2,  # ← Reduced transparency for shaded region
                     label=r'$\pm 1 \sigma$ uncertainty',
                     zorder=2)
    # # Plot with error bars
    # ax.errorbar(r_plot, radial_profile_plot, 
    #             yerr=radial_std_plot, 
    #             fmt='o-', 
    #             color=user_params.get('color', 'blue'), 
    #             markersize=6, 
    #             capsize=5, 
    #             alpha=0.75,
    #             linewidth=2,
    #             label='Volume-weighted mean')
    
    # Set y-axis limits
    clim = user_params.get('clim', [None, None])
    y_min = clim[0] if clim[0] is not None else np.nanmin(radial_profile_plot)
    y_max = clim[1] if clim[1] is not None else np.nanmax(radial_profile_plot)
    ax.set_ylim([y_min, y_max])
    
    # Labels and grid
    ax.set_xlabel(user_params.get('xlabel', r'$r$'), fontsize=12)
    ax.set_ylabel(user_params.get('cmap_label', 'Value'), fontsize=12)
    ax.grid()
    ax.legend()
    
    # Create output directory
    variable = user_params.get('variable', 'unknown')
    if variable == "eint":
        out_dir = Path(user_params['output_path']) / "spherical_profile_3d" / "pres"
    elif variable.startswith("derived:"):
        var_name = variable.split(":", 1)[1].strip()
        out_dir = Path(user_params['output_path']) / "spherical_profile_3d" / var_name
    else:
        out_dir = Path(user_params['output_path']) / "spherical_profile_3d" / variable
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Save figure
    if user_params.get('loop_bin', False):
        fname = os.path.join(out_dir, f"t_{user_params.get('slice_number', 0)}.png")
    else:
        fname = os.path.join(out_dir, 'spherical_profile_3d.png')
    
    plt.tight_layout()
    plt.savefig(fname, bbox_inches="tight", dpi=150)
    plt.close()
    
    print(f"\n 3D Spherical profile saved at {fname}")
    
    # Return data for further analysis
    return profile

def plot_spherical_mass_weighted_avg_profile(profile, user_params):
    """
    Computes and plots 3D spherically-averaged MASS-WEIGHTED profile.
    Ideal for velocities and momentum-dependent quantities.
    
    Parameters:
    -----------
    profile : dict
        Output from data_processing.block_data.compute_radial_profile with keys:
        r_centers, profile, std, count, sum
    
    user_params : dict
        Must include:
        - center_x, center_y, center_z: sphere center
        - num_radial_bins: number of radial bins
        - density_file: (optional) filename for density if not already extracted
        - Other standard params: color, norm, xlabel, cmap_label, etc.
    """    
    # ===== Create figure =====
    fig, ax = plt.subplots(figsize=(12, 7))
    
    axes_scale = user_params.get('axes_scale', 1.0)
    # Get plot color
    color = user_params.get('color', 'blue')
    
    # Apply scaling and normalization
    if user_params.get('norm') == "log":
        if user_params['xlabel'].endswith("[log]"):
            r_plot = np.log10(axes_scale * profile['r_centers'] + 1e-30)
        else:
            r_plot = axes_scale * profile['r_centers']
            
        radial_profile_plot = np.log10(np.abs(profile['profile']) + 1e-30)
        # For log scale error transformation
        radial_std_plot = profile['std'] / (np.abs(profile['profile']) * np.log(10))
        
        
    else:
        radial_profile_plot = profile['profile']
        radial_std_plot = profile['std']
        r_plot = axes_scale * profile['r_centers']
    
    
    
    # Plot line with shaded uncertainty region
    ax.plot(r_plot, radial_profile_plot, 
            color=color, 
            linewidth=2.5,
            marker='o',
            markersize=6,
            label='Mass-weighted mean',
            zorder=3)
    
    # Add shaded region (std as confidence band)
    ax.fill_between(r_plot, 
                     radial_profile_plot - radial_std_plot,
                     radial_profile_plot + radial_std_plot,
                     color=color, 
                     alpha=0.2,  # ← Reduced transparency for shaded region
                     label=r'$\pm 1 \sigma$ uncertainty',
                     zorder=2)
    
    # Set y-axis limits
    clim = user_params.get('clim', [None, None])
    y_min = clim[0] if clim[0] is not None else np.nanmin(radial_profile_plot)
    y_max = clim[1] if clim[1] is not None else np.nanmax(radial_profile_plot)
    ax.set_ylim([y_min, y_max])
    
    # Labels and grid
    ax.set_xlabel(user_params.get('xlabel', r'$r$'), fontsize=12)
    ax.set_ylabel(user_params.get('cmap_label', 'Value'), fontsize=12)
    ax.grid()
    ax.legend()
    
    # Create output directory
    variable = user_params.get('variable', 'unknown')
    if variable.startswith("derived:"):
        var_name = variable.split(":", 1)[1].strip()
        out_dir = Path(user_params['output_path']) / "spherical_profile_3d" / var_name
    else:
        out_dir = Path(user_params['output_path']) / "spherical_profile_3d" / variable
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Save figure
    if user_params.get('loop_bin', False):
        fname = os.path.join(out_dir, f"t_{user_params.get('slice_number', 0)}.png")
    else:
        fname = os.path.join(out_dir, 'spherical_profile_3d_mass_weighted.png')
    
    plt.tight_layout()
    plt.savefig(fname, bbox_inches="tight", dpi=150)
    plt.close()
    
    print(f"✓ 3D Mass-weighted profile saved at {fname}")
    
    return profile

def plot_time_scales_profile(timescales, user_params):
     # ===== Create figure =====
    fig, ax = plt.subplots(figsize=(12, 7))
    all_profiles = []
    axes_scale = user_params.get('axes_scale', 1.0)
    for time_scale in timescales:
        time_scale_profile  = time_scale['profile']
        time_scale_label    = time_scale['label']
        time_scale_color    = time_scale['color']
        # Process radii for plotting (apply log scale if needed)
        if user_params.get('norm') == "log":
            if user_params['xlabel'].endswith("[log]"):
                r_plot      = np.log10(axes_scale * time_scale_profile['r_centers'] + 1e-30)
            else:
                r_plot      = axes_scale * time_scale_profile['r_centers']
            # Log scale for y-axis
            profile_plot    = np.log10(time_scale_profile['profile'] + 1e-30)
            std_plot        = time_scale_profile['std'] / (time_scale_profile['profile'] * np.log(10))
        else:
            r_plot          = axes_scale * time_scale_profile['r_centers']
            profile_plot    = time_scale_profile['profile']
            std_plot        = time_scale_profile['std']
        all_profiles.append(profile_plot)
        # Plot
        ax.plot(r_plot, profile_plot, 
                color       = time_scale_color, 
                linestyle   = '-',  # SOLID
                linewidth   = 2.5,
                marker      = 'o',
                markersize  = 4,
                label       = time_scale_label,
                zorder      = 3)
    
        # Add shaded region
        ax.fill_between(r_plot, 
                         profile_plot - std_plot,
                         profile_plot + std_plot,
                         color = time_scale_color, 
                         alpha = 0.1,
                         zorder = 1)    
    # Set y-axis limits
    clim = user_params.get('clim', [None, None])
    all_profiles_list = np.concatenate(all_profiles)
    y_min = clim[0] if clim[0] is not None else np.nanmin(all_profiles_list)
    y_max = clim[1] if clim[1] is not None else np.nanmax(all_profiles_list)
    ax.set_ylim([y_min, y_max])
    
    # Labels and grid
    ax.set_xlabel(user_params.get('xlabel', r'$r$'), fontsize=12)
    ax.set_ylabel(user_params.get('cmap_label', 'Value'), fontsize=12)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=11, loc='best')
    
    # Create output directory
    variable = user_params.get('variable', 'unknown')
    if variable.startswith("derived:"):
        var_name = variable.split(":", 1)[1].strip()
        out_dir = Path(user_params['output_path']) / "spherical_profile_3d" / var_name
    else:
        out_dir = Path(user_params['output_path']) / "spherical_profile_3d" / variable
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Save figure
    if user_params.get('loop_bin', False):
        fname = os.path.join(out_dir, f"t_{user_params.get('slice_number', 0)}.png")
    else:
        fname = os.path.join(out_dir, 'spherical_profile_3d_mass_weighted.png')
    
    plt.tight_layout()
    plt.savefig(fname, bbox_inches="tight", dpi=150)
    plt.close()
    
    print(f"✓ 3D Mass-weighted profile saved at {fname}")