# project_root/plotting/streamlines.py

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

# Ensure parent directory (where bin_convert.py lives) is on sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..', '..'))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from data_processing.slice_data import stitch_meshblocks_to_global  # assumes stitch_data in data_processing.py

def plot_streamlines(X, Y, U, V, user_params, 
    quiver_stride=4,
    quiver_scale=50,
    headwidth=3,
    headlength=5,
    headaxislength=4,
    log_vmin=1e-3):
    """
    Plot streamlines from global 2D velocity arrays.

    Parameters:
        X, Y : 2D numpy arrays of grid coordinates
        U, V : 2D numpy arrays of velocity components on the grid
        user_params (dict): Optional plotting parameters:
            - 'cmap': colormap (default 'viridis')
            - 'density': streamline density (default 1.0)
            - 'output_path': where to save the figure
            - 'input_file': base name for saving
            - 'slice_number': for naming in loop mode
            - 'dt': time step, used in title if provided
            - 'extents': [xmin,xmax,ymin,ymax] for axis limits
    """
    cmap    = user_params['cmap']
    density = user_params['streamline_density']
    figsize = user_params['figsize']
    
    # Mask out grid points without data (NaNs) so streamplot ignores them
    U_masked = np.where(np.isnan(U), np.nan, U)
    V_masked = np.where(np.isnan(V), np.nan, V)
    
    speed = np.sqrt(U_masked**2 + V_masked**2)
    max_speed = np.nanmax(speed)
    mask = speed <= 0
    speed_mask = np.where(mask, np.nan, speed)
    linewidth = (1.5 * speed / max_speed) if max_speed > 0 else 1.0

    fig, ax = plt.subplots(figsize=figsize)
    # Set up normalization
    norm = LogNorm(vmin=log_vmin, vmax=np.nanmax(speed_mask)) if user_params['norm']=="log" else None
    strm = ax.streamplot(
        X, Y, U_masked, V_masked,
        color=speed_mask,
        linewidth=linewidth,
        cmap=cmap,
        norm=norm,
        density=density,
        arrowsize=0,           # hide built-in arrows
        integration_direction='both'
    )
    
    # 2) Overlay quiver:
    #    Sample every `quiver_stride` points
    xs = X[::quiver_stride, ::quiver_stride]
    ys = Y[::quiver_stride, ::quiver_stride]
    us = U_masked[::quiver_stride, ::quiver_stride]
    vs = V_masked[::quiver_stride, ::quiver_stride]
    
    # Compute speed at quiver points
    speed_quiver = np.sqrt(us**2 + vs**2)

    # Create a mask for very small speeds
    small_speed_mask = speed_quiver < 1e-10

    # Replace these small vectors with zero-length arrows (or remove them)
    us = np.where(small_speed_mask, np.nan, us)
    vs = np.where(small_speed_mask, np.nan, vs)

    q = ax.quiver(
        xs, ys, us, vs,
        angles='xy', scale_units='xy', scale=quiver_scale,
        width=0.002,
        headwidth=headwidth,
        headlength=headlength,
        headaxislength=headaxislength,
        color='k'
    )
    cbar = fig.colorbar(strm.lines, ax=ax, norm=norm)
    cbar.set_label('Speed')

    ax.set_xlabel('x')
    ax.set_ylabel('z')

    # Title includes time if dt and slice_number are present
    if 'slice_number' in user_params and 'dt' in user_params:
        t = user_params['slice_number'] * user_params['dt']
        ax.set_title(f"Streamlines at t = {t:.2f}")
    else:
        ax.set_title(user_params.get('title', 'Streamlines'))

    ax.set_aspect('equal')
    
    # Force full-domain limits
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())

    # Then override if user explicitly provided extents
    ext = user_params['extents']
    xmin_u, xmax_u, ymin_u, ymax_u = ext
    if xmin_u is not None and xmax_u is not None:
        ax.set_xlim(xmin_u, xmax_u)
    if ymin_u is not None and ymax_u is not None:
        ax.set_ylim(ymin_u, ymax_u)

    # Save figure if output_path given
    out_path = user_params.get('output_path')
    if out_path:
        save_dir = Path(out_path) / 'streamlines'
        save_dir.mkdir(parents=True, exist_ok=True)

        if user_params.get('loop_bin', False):
            fname = f"streamlines_{user_params['slice_number']}.png"
        else:
            fname = f"streamlines_{user_params.get('input_file', 'output')}.png"

        fig.savefig(save_dir / fname, dpi=300, bbox_inches='tight')
        print(f"Saved streamlines plot to {save_dir / fname}")
    plt.close()


def plot_streamlines_from_dataframes(df_u, df_v, user_params):
    """
    Stitch u and v block data into global arrays, then plot streamlines.

    Parameters:
        df_u : DataFrame with 'extents' and 'data2d' for u-velocity
        df_v : DataFrame with 'extents' and 'data2d' for v-velocity
        user_params (dict): passed to plot_streamlines()
    """
    # Stitch data
    X, Y, U = stitch_meshblocks_to_global(df_u)
    _, _, V = stitch_meshblocks_to_global(df_v)

    # Plot
    plot_streamlines(X, Y, U, V, user_params)
