import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from scipy.signal import find_peaks
from data_processing.block_data import extract_athenak_3D_block

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

def plot_spherical_volume_weighted_avg_profile(df3D, user_params):
    """
    Computes and plots 3D spherically-averaged profile from extract_athenak_3D_block output.
    Uses fast per-block NumPy arrays instead of slow DataFrame reconstruction.
    
    Parameters:
    -----------
    df3D : dict
        Output from extract_athenak_3D_block containing:
        - blocks: list of dicts with 'data', 'x', 'y', 'z', 'extent'
        - df_extents: DataFrame with block extents
        - num_blocks: number of blocks
        - block_shape: (nz, ny, nx)
    
    user_params : dict
        Must include center_x, center_y, center_z, num_radial_bins, etc.
    """
    blocks = df3D['blocks']
    df_extents = df3D['df_extents']
    nz, ny, nx = df3D['block_shape']
    
    # Extract sphere center
    x_center = user_params.get('center_x', 0.0)
    y_center = user_params.get('center_y', 0.0)
    z_center = user_params.get('center_z', 0.0)
    num_bins = user_params.get('num_radial_bins', 100)
    
    # print(f"Computing 3D spherical average centered at ({x_center}, {y_center}, {z_center})")
    # print(f"Block shape: nz={nz}, ny={ny}, nx={nx}")
    # print(f"Number of Blocks = {len(blocks)}")
    
    radii_all = []
    data_all = []
    volumes_all = []
    
    # ===== Iterate through each block using fast NumPy access =====
    for block_idx, block_dict in enumerate(blocks):
        try:
            # Get 3D numpy array and coordinates directly (NO MultiIndex extraction!)
            data_3d = block_dict['data']          # shape (nz, ny, nx)
            x_coords = block_dict['x']             # shape (nx,)
            y_coords = block_dict['y']             # shape (ny,)
            z_coords = block_dict['z']             # shape (nz,)
            x_min, x_max, y_min, y_max, z_min, z_max = block_dict['extent']
            
            # Create 3D coordinate grids (cell centers)
            # Note: meshgrid with indexing='ij' gives (nx, ny, nz)
            # We need to transpose to match data shape (nz, ny, nx)
            xx, yy, zz = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
            
            # Reorder axes: meshgrid gives (nx, ny, nz), we want (nz, ny, nx)
            xx = np.transpose(xx, (2, 1, 0))  # (nz, ny, nx)
            yy = np.transpose(yy, (2, 1, 0))  # (nz, ny, nx)
            zz = np.transpose(zz, (2, 1, 0))  # (nz, ny, nx)
            
            # Compute 3D radial distance from center
            rr = np.sqrt((xx - x_center)**2 + (yy - y_center)**2 + (zz - z_center)**2)
            
            # Cell volume
            dx = (x_max - x_min) / nx if nx > 1 else 1.0
            dy = (y_max - y_min) / ny if ny > 1 else 1.0
            dz = (z_max - z_min) / nz if nz > 1 else 1.0
            cell_volume = dx * dy * dz
            
            # Flatten and accumulate (very fast on NumPy arrays)
            radii_all.append(rr.ravel())
            data_all.append(data_3d.ravel())
            volumes_all.append(np.full(data_3d.size, cell_volume, dtype=np.float64))
            
            # print(f"  Block {block_idx}: extent=[{x_min:.3e}, {x_max:.3e}] x " +
            #       f"[{y_min:.3e}, {y_max:.3e}] x [{z_min:.3e}, {z_max:.3e}]")
            
        except Exception as e:
            print(f"  ✗ Error processing block {block_idx}: {e}")
            continue
    
    # Concatenate all blocks (still fast with NumPy)
    radii = np.concatenate(radii_all)
    data = np.concatenate(data_all)
    volumes = np.concatenate(volumes_all)
    
    # print(f"\nTotal cells accumulated: {len(data)}")
    # print(f"Data range: [{np.nanmin(data):.3e}, {np.nanmax(data):.3e}]")
    # print(f"Radii range: [{np.nanmin(radii):.3e}, {np.nanmax(radii):.3e}]")
    
    # Remove NaN and invalid data
    mask = np.isfinite(data) & np.isfinite(radii) & (radii >= 0)
    radii = radii[mask]
    data = data[mask]
    volumes = volumes[mask]
    
    # print(f"Valid cells after filtering: {len(data)}")
    
    # Define radial bins
    r_min = np.min(radii)
    r_max = np.max(radii)
    r_bins = np.linspace(r_min, r_max, num_bins + 1)
    r_centers = (r_bins[:-1] + r_bins[1:]) / 2
    
    # print(f"Radial binning: {num_bins} bins from {r_min:.3e} to {r_max:.3e}")
    
    # Bin and compute volume-weighted averages (vectorized)
    radial_profile = np.empty(num_bins)
    radial_std = np.empty(num_bins)
    radial_count = np.zeros(num_bins, dtype=int)
    radial_sum = np.empty(num_bins)
    
    for i in range(num_bins):
        mask_shell = (radii >= r_bins[i]) & (radii < r_bins[i+1])
        n_cells = np.sum(mask_shell)
        
        if n_cells > 0:
            data_in_shell = data[mask_shell]
            volumes_in_shell = volumes[mask_shell]
            
            # Volume-weighted mean
            radial_profile[i] = np.average(data_in_shell, weights=volumes_in_shell)
            radial_sum[i] = np.sum(data_in_shell * volumes_in_shell)
            radial_std[i] = np.std(data_in_shell)
            radial_count[i] = n_cells
        else:
            radial_profile[i] = np.nan
            radial_sum[i] = np.nan
            radial_std[i] = np.nan
            radial_count[i] = 0
    
    # Filter out empty bins
    valid_mask = radial_count > 0
    r_centers_valid = r_centers[valid_mask]
    radial_profile_valid = radial_profile[valid_mask]
    radial_std_valid = radial_std[valid_mask]
    radial_sum_valid = radial_sum[valid_mask]
    radial_count_valid = radial_count[valid_mask]
    
    # print(f"Bins with data: {np.sum(valid_mask)}")
    
    # ===== Create figure =====
    fig, ax = plt.subplots(figsize=(12, 7))
    
    axes_scale = user_params.get('axes_scale', 1.0)
    
    # Apply scaling and normalization
    if user_params.get('norm') == "log":
        radial_profile_plot = np.log10(radial_profile_valid + 1e-30)
        # For log scale: if y = log10(Y), then dy ≈ dY / (Y * ln(10))
        radial_std_plot = radial_std_valid / (radial_profile_valid * np.log(10))

        if user_params['xlabel'].endswith("[log]"):
            r_plot = np.log10(axes_scale * r_centers_valid + 1e-30)
        else:
            r_plot = axes_scale * r_centers_valid
    else:
        radial_profile_plot = radial_profile_valid
        radial_std_plot = radial_std_valid
        r_plot = axes_scale * r_centers_valid
    
    # Get plot color
    color = user_params.get('color', 'blue')
    
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
    
    print(f"\n✓ 3D Spherical profile saved at {fname}")
    
    # Return data for further analysis
    return {
        "r_centers": r_centers_valid,
        "profile": radial_profile_valid,
        "std": radial_std_valid,
        "count": radial_count_valid,
        "sum": radial_sum_valid
    }

def plot_spherical_mass_weighted_avg_profile(df3D, user_params):
    """
    Computes and plots 3D spherically-averaged MASS-WEIGHTED profile.
    Ideal for velocities and momentum-dependent quantities.
    
    Parameters:
    -----------
    df3D : dict
        Output from extract_athenak_3D_block containing:
        - blocks: list of dicts with 'data', 'x', 'y', 'z', 'extent'
        - df_extents: DataFrame with block extents
        - num_blocks: number of blocks
        - block_shape: (nz, ny, nx)
    
    user_params : dict
        Must include:
        - center_x, center_y, center_z: sphere center
        - num_radial_bins: number of radial bins
        - density_file: (optional) filename for density if not already extracted
        - Other standard params: color, norm, xlabel, cmap_label, etc.
    """
    blocks = df3D['blocks']
    nz, ny, nx = df3D['block_shape']
    
    # Extract sphere center
    x_center = user_params.get('center_x', 0.0)
    y_center = user_params.get('center_y', 0.0)
    z_center = user_params.get('center_z', 0.0)
    num_bins = user_params.get('num_radial_bins', 100)
    
    # === Extract density for mass-weighting ===
    rho_params = user_params.copy()
    rho_params['variable'] = 'dens'
    rho_df3D = extract_athenak_3D_block(rho_params)
    rho_blocks = rho_df3D['blocks']
    
    radii_all = []
    data_all = []
    masses_all = []  # ← Mass weights instead of volumes!
    
    # ===== Iterate through each block =====
    for block_idx, (block_dict, rho_block_dict) in enumerate(zip(blocks, rho_blocks)):
        try:
            # Get data and density
            data_3d = block_dict['data']          # shape (nz, ny, nx)
            rho_3d = rho_block_dict['data']       # shape (nz, ny, nx)
            
            x_coords = block_dict['x']
            y_coords = block_dict['y']
            z_coords = block_dict['z']
            x_min, x_max, y_min, y_max, z_min, z_max = block_dict['extent']
            
            # Create 3D coordinate grids
            xx, yy, zz = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
            
            # Reorder axes to match data shape (nz, ny, nx)
            xx = np.transpose(xx, (2, 1, 0))
            yy = np.transpose(yy, (2, 1, 0))
            zz = np.transpose(zz, (2, 1, 0))
            
            # Compute radial distance
            rr = np.sqrt((xx - x_center)**2 + (yy - y_center)**2 + (zz - z_center)**2)
            
            # Cell volume
            dx = (x_max - x_min) / nx if nx > 1 else 1.0
            dy = (y_max - y_min) / ny if ny > 1 else 1.0
            dz = (z_max - z_min) / nz if nz > 1 else 1.0
            cell_volume = dx * dy * dz
            
            # === KEY DIFFERENCE: Mass = rho * volume ===
            cell_mass = rho_3d * cell_volume
            
            # Flatten and accumulate
            radii_all.append(rr.ravel())
            data_all.append(data_3d.ravel())
            masses_all.append(cell_mass.ravel())  # ← Mass per cell, not constant volume!
            
        except Exception as e:
            print(f"✗ Error processing block {block_idx}: {e}")
            continue
    
    # Concatenate all blocks
    radii = np.concatenate(radii_all)
    data = np.concatenate(data_all)
    masses = np.concatenate(masses_all)
    
    # Remove NaN and invalid data
    mask = np.isfinite(data) & np.isfinite(radii) & np.isfinite(masses) & (radii >= 0) & (masses > 0)
    radii = radii[mask]
    data = data[mask]
    masses = masses[mask]
    
    # Define radial bins
    r_min = np.min(radii)
    r_max = np.max(radii)
    r_bins = np.linspace(r_min, r_max, num_bins + 1)
    r_centers = (r_bins[:-1] + r_bins[1:]) / 2
    
    # === Bin and compute MASS-weighted averages ===
    radial_profile = np.empty(num_bins)
    radial_std = np.empty(num_bins)
    radial_count = np.zeros(num_bins, dtype=int)
    radial_sum = np.empty(num_bins)
    
    for i in range(num_bins):
        mask_shell = (radii >= r_bins[i]) & (radii < r_bins[i+1])
        n_cells = np.sum(mask_shell)
        
        if n_cells > 0:
            data_in_shell = data[mask_shell]
            masses_in_shell = masses[mask_shell]
            
            # === MASS-weighted mean ===
            radial_profile[i] = np.average(data_in_shell, weights=masses_in_shell)
            radial_sum[i] = np.sum(data_in_shell * masses_in_shell)
            radial_std[i] = np.std(data_in_shell)  # std of values, not weighted
            radial_count[i] = n_cells
        else:
            radial_profile[i] = np.nan
            radial_sum[i] = np.nan
            radial_std[i] = np.nan
            radial_count[i] = 0
    
    # Filter out empty bins
    valid_mask = radial_count > 0
    r_centers_valid = r_centers[valid_mask]
    radial_profile_valid = radial_profile[valid_mask]
    radial_std_valid = radial_std[valid_mask]
    radial_sum_valid = radial_sum[valid_mask]
    radial_count_valid = radial_count[valid_mask]
    
    # ===== Create figure =====
    fig, ax = plt.subplots(figsize=(12, 7))
    
    axes_scale = user_params.get('axes_scale', 1.0)
    
    # Apply scaling and normalization
    if user_params.get('norm') == "log":
        radial_profile_plot = np.log10(np.abs(radial_profile_valid) + 1e-30)
        # For log scale error transformation
        radial_std_plot = radial_std_valid / (np.abs(radial_profile_valid) * np.log(10))
        
        if user_params['xlabel'].endswith("[log]"):
            r_plot = np.log10(axes_scale * r_centers_valid + 1e-30)
        else:
            r_plot = axes_scale * r_centers_valid
    else:
        radial_profile_plot = radial_profile_valid
        radial_std_plot = radial_std_valid
        r_plot = axes_scale * r_centers_valid
    
    # Get plot color
    color = user_params.get('color', 'blue')
    
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
    
    return {
        "r_centers": r_centers_valid,
        "profile": radial_profile_valid,
        "std": radial_std_valid,
        "count": radial_count_valid,
        "sum": radial_sum_valid
    }
