import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from data_processing.profile_data import get_zh_data

def plot_comparative_profile(profile_no_rad, profile_rad, user_params):
    """
    Computes and plots comparative 3D spherically-averaged profiles:
    - Without radiative cooling (df3D) as DOTTED lines
    - With radiative cooling (raddf3D) as SOLID lines
    
    Parameters:
    -----------
    profile_no_rad : dict
        Output from data_processing.block_data.compute_radial_profile with keys:
        r_centers, profile, std, count, sum
    profile_rad : dict
        Output from data_processing.block_data.compute_radial_profile with keys:
        r_centers, profile, std, count, sum
    user_params : dict
        Must include center_x, center_y, center_z, num_radial_bins, etc.
    """
    
    # ===== Create figure =====
    fig, ax = plt.subplots(figsize=(12, 7))
    
    axes_scale = user_params.get('axes_scale', 1.0)
    color = user_params.get('color', 'blue')
    
    # Process radii for plotting (apply log scale if needed)
    if user_params.get('norm') == "log":
        if user_params['xlabel'].endswith("[log]"):
            r_plot_no_rad = np.log10(axes_scale * profile_no_rad['r_centers'] + 1e-30)
            r_plot_rad = np.log10(axes_scale * profile_rad['r_centers'] + 1e-30)
        else:
            r_plot_no_rad = axes_scale * profile_no_rad['r_centers']
            r_plot_rad = axes_scale * profile_rad['r_centers']
        
        # Log scale for y-axis
        profile_plot_no_rad = np.log10(profile_no_rad['profile'] + 1e-30)
        profile_plot_rad = np.log10(profile_rad['profile'] + 1e-30)
        std_plot_no_rad = profile_no_rad['std'] / (profile_no_rad['profile'] * np.log(10))
        std_plot_rad = profile_rad['std'] / (profile_rad['profile'] * np.log(10))
    else:
        r_plot_no_rad = axes_scale * profile_no_rad['r_centers']
        r_plot_rad = axes_scale * profile_rad['r_centers']
        profile_plot_no_rad = profile_no_rad['profile']
        profile_plot_rad = profile_rad['profile']
        std_plot_no_rad = profile_no_rad['std']
        std_plot_rad = profile_rad['std']
    
    # Plot WITHOUT radiative cooling - DOTTED line
    ax.plot(r_plot_no_rad, profile_plot_no_rad, 
            color=color, 
            linestyle=':',  # DOTTED
            linewidth=2.5,
            marker='o',
            markersize=4,
            label='Without radiative cooling',
            zorder=3)
    
    # Add shaded region for no-rad case
    ax.fill_between(r_plot_no_rad, 
                     profile_plot_no_rad - std_plot_no_rad,
                     profile_plot_no_rad + std_plot_no_rad,
                     color=color, 
                     alpha=0.1,
                     zorder=1)
    
    # Plot WITH radiative cooling - SOLID line
    ax.plot(r_plot_rad, profile_plot_rad, 
            color=color, 
            linestyle='-',  # SOLID
            linewidth=2.5,
            marker='s',
            markersize=4,
            label='With radiative cooling',
            zorder=3)
    
    # Add shaded region for rad case
    ax.fill_between(r_plot_rad, 
                     profile_plot_rad - std_plot_rad,
                     profile_plot_rad + std_plot_rad,
                     color=color, 
                     alpha=0.2,
                     zorder=2)
    
    # Set y-axis limits
    clim = user_params.get('clim', [None, None])
    all_profiles = np.concatenate([profile_plot_no_rad, profile_plot_rad])
    y_min = clim[0] if clim[0] is not None else np.nanmin(all_profiles)
    y_max = clim[1] if clim[1] is not None else np.nanmax(all_profiles)
    ax.set_ylim([y_min, y_max])
    
    # Labels and grid
    ax.set_xlabel(user_params.get('xlabel', r'$r$'), fontsize=12)
    ax.set_ylabel(user_params.get('cmap_label', 'Value'), fontsize=12)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=11, loc='best')
    
    # Create output directory
    variable = user_params.get('variable', 'unknown')
    if variable == "eint":
        out_dir = Path(user_params['output_path']) / "comparision_profile" / "pres"
    elif variable.startswith("derived:"):
        var_name = variable.split(":", 1)[1].strip()
        out_dir = Path(user_params['output_path']) / "comparision_profile" / var_name
    else:
        out_dir = Path(user_params['output_path']) / "comparision_profile" / variable
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Save figure
    if user_params.get('loop_bin', False):
        fname = os.path.join(out_dir, f"t_{user_params.get('slice_number', 0)}_comparison.png")
    else:
        fname = os.path.join(out_dir, 'spherical_profile_3d_comparison.png')
    
    plt.tight_layout()
    plt.savefig(fname, bbox_inches="tight", dpi=150)
    plt.close()
    
    print(f"\n✓ Comparative 3D Spherical profile saved at {fname}")
    
    # Return comparison data
    return {
        "r_centers_no_rad": profile_no_rad['r_centers'],
        "profile_no_rad": profile_no_rad['profile'],
        "std_no_rad": profile_no_rad['std'],
        "r_centers_rad": profile_rad['r_centers'],
        "profile_rad": profile_rad['profile'],
        "std_rad": profile_rad['std'],
    }
    
def plot_comparative_time_scales(timescales, user_params):
     # ===== Create figure =====
    fig, ax = plt.subplots(figsize=(12, 7))
    all_profiles = []
    axes_scale = user_params.get('axes_scale', 1.0)
    for time_scale in timescales:
        time_scale_profile  = time_scale['profile']
        time_scale_label    = time_scale['label']
        time_scale_color    = time_scale['color']
        time_scale_style    = time_scale['style']
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
                linestyle   = time_scale_style,
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
    if variable == "eint":
        out_dir = Path(user_params['output_path']) / "comparision_profile" / "pres"
    elif variable.startswith("derived:"):
        var_name = variable.split(":", 1)[1].strip()
        out_dir = Path(user_params['output_path']) / "comparision_profile" / var_name
    else:
        out_dir = Path(user_params['output_path']) / "comparision_profile" / variable
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Save figure
    if user_params.get('loop_bin', False):
        fname = os.path.join(out_dir, f"t_{user_params.get('slice_number', 0)}_comparison.png")
    else:
        fname = os.path.join(out_dir, 'spherical_profile_3d_comparison.png')
    
    plt.tight_layout()
    plt.savefig(fname, bbox_inches="tight", dpi=150)
    plt.close()
    
    print(f"\n✓ Comparative 3D Spherical profile saved at {fname}")



def plot_comparative_profile_uCGM_tCGM(uCGM_df, tCGM_df, uCGM_params, tCGM_params):
    if uCGM_params['profile_variable'] != 'zh':
        return

    uCGM_x, uCGM_data = get_zh_data(uCGM_df, uCGM_params)
    tCGM_x, tCGM_data = get_zh_data(tCGM_df, tCGM_params)

    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update({
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.labelweight": "bold",
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 13,
        "axes.linewidth": 1.2,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    })

    fig, ax = plt.subplots(figsize=uCGM_params.get("figsize", (10, 8)))

    sc1 = ax.scatter(
        uCGM_x,
        uCGM_data,
        color=uCGM_params.get("color", "tab:blue"),
        s=uCGM_params.get("marker_size", 20),
        alpha=0.75,
        edgecolors="none",
        rasterized=True,
        label="uCGM"
    )

    sc2 = ax.scatter(
        tCGM_x,
        tCGM_data,
        color=tCGM_params.get("color", "tab:orange"),
        s=tCGM_params.get("marker_size", 20),
        alpha=0.75,
        edgecolors="none",
        rasterized=True,
        label="tCGM"
    )

    ax.set_xlabel(uCGM_params["xlabel"], labelpad=10)
    ax.set_ylabel(uCGM_params["cmap_label"], labelpad=10)
    ax.set_xlim(*uCGM_params["clim"])

    ax.tick_params(direction="in", length=5, width=1.1, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)

    ax.grid(alpha=0.25, linewidth=0.8)

    legend = ax.legend(
        handles=[sc1, sc2],
        labels=["uniform CGM", "exponential CGM"],
        loc="lower left",
        bbox_to_anchor=(0.02, 0.02),
        frameon=True,
        fancybox=True,
        framealpha=0.95,
        edgecolor="0.8",
        borderpad=0.6,
        handletextpad=0.6,
        scatterpoints=1,
        markerscale=1.4,
    )

    legend.get_frame().set_linewidth(0.8)

    fig.tight_layout()

    out_dir = Path(uCGM_params["output_path"]) / "compare_zh_profile"
    out_dir.mkdir(parents=True, exist_ok=True)

    if uCGM_params.get("loop_bin", False):
        fname = out_dir / f"zh_{uCGM_params['slice_number']}.png"
    else:
        fname = out_dir / "zh.png"

    fig.savefig(fname, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

    print(f"Profile saved at {fname}")