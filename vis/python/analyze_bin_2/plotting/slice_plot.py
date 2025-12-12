import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
import matplotlib.colors as plt_col
from pathlib import Path
import os
import numpy as np

def set_normalization(user_params, vmin, vmax):
    if user_params['norm'] is None:
        return plt_col.Normalize(vmin, vmax)
    elif user_params['norm'] == "log":
        return plt_col.LogNorm(vmin, vmax)
    
import numpy as np
import matplotlib.pyplot as plt
    
def plot_athenak_combined(data_dict, user_params):
    """
    Combine all blocks' 2D slices spatially and plot as a single image with colorbar and save.
    
    Parameters:
    - data_dict: dict output from extract_athenak_slice_as_dataframe, with keys:
        'df_quantities', 'df_extents', 'block_shape'
    - user_params: dict for plotting options including cmap, figsize, xlabel, ylabel,
                   axes_scale, extents, output_path, variable, input_file, slice_number, dt
    - norm: matplotlib normalization instance controlling colormap scaling

    Saves the plot as PNG to specified output path and closes figure.
    """
    df_quantities = data_dict['df_quantities']
    df_extents = data_dict['df_extents']
    block_nx2, block_nx1 = data_dict['block_shape']
    
    global_min = df_quantities.apply(np.min).min()
    global_max = df_quantities.apply(np.max).max()
    
    vmin = user_params['clim'][0] if user_params['clim'][0] is not None else global_min
    vmax = user_params['clim'][1] if user_params['clim'][1] is not None else global_max
    norm = set_normalization(user_params, vmin, vmax)

    scaled_extents = [tuple(user_params['axes_scale'] * x for x in extent) for extent in df_extents.values]

    # Initialize figure and axis
    fig, ax = plt.subplots(figsize=user_params.get('figsize', (10, 10)))
    
    for block_num, extent in enumerate(scaled_extents):
        block_data = df_quantities.loc[block_num].unstack(level=-1).values
        
        ax.imshow(
            block_data,
            cmap=user_params['cmap'],
            norm=norm,
            interpolation='none',
            origin='lower',
            extent=extent,
            aspect='auto'
        )

    all_extents = np.vstack(scaled_extents)
    user_extents = user_params.get('extents', [None, None, None, None])

    ax.set_xlim(
        user_extents[0] if user_extents[0] is not None else np.min(all_extents[:, 0]),
        user_extents[1] if user_extents[1] is not None else np.max(all_extents[:, 1])
    )
    ax.set_ylim(
        user_extents[2] if user_extents[2] is not None else np.min(all_extents[:, 2]),
        user_extents[3] if user_extents[3] is not None else np.max(all_extents[:, 3])
    )

    ax.set_aspect('equal')
    ax.set_xlabel(user_params['xlabel'])
    ax.set_ylabel(user_params['ylabel'])

    slice_no = user_params.get('slice_number', 0)
    dt = user_params.get('dt', 0)
    ax.set_title(f't = {slice_no * dt:.3f} Myr')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    sm = ScalarMappable(norm=norm, cmap=user_params['cmap'])
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label(user_params['cmap_label'])

    out_dir = Path(user_params['output_path']) / user_params['variable']
    out_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{user_params['input_file']}.png"
    out_path = out_dir / fname

    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to {out_path}")

def plot_stitched_data(global_array, global_extent, user_params):
    """
    Plot the stitched global 2D AthenaK slice data.
    
    Parameters:
    - global_array: 2D numpy array from stitch_meshblocks_to_global
    - global_extent: tuple (x_min, x_max, y_min, y_max) corresponding to global_array
    - user_params: dict with keys 'cmap', 'xlabel', 'ylabel', 'figsize', 'cmap_label'
    - norm: matplotlib color normalization instance for consistent colormap scaling
    
    Behavior:
    - Plots the global 2D array with imshow over the specified extent.
    - Uses user colormap and normalization.
    - Sets axis labels and colorbar label.
    """
    
    global_min = global_array.min()
    global_max = global_array.max()
    
    vmin = user_params['clim'][0] if user_params['clim'][0] is not None else global_min
    vmax = user_params['clim'][1] if user_params['clim'][1] is not None else global_max
    norm = set_normalization(user_params, vmin, vmax)
    fig, ax = plt.subplots(figsize=user_params.get('figsize', (10, 8)))
    im = ax.imshow(global_array, origin='lower', extent=global_extent,
                   cmap=user_params['cmap'], norm=norm, aspect='auto')
    
    ax.set_xlabel(user_params['xlabel'])
    ax.set_ylabel(user_params['ylabel'])
    ax.set_aspect('equal')
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    sm = ScalarMappable(norm=norm, cmap=user_params['cmap'])
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label(user_params['cmap_label'])
    
    plt.tight_layout()
    
    out_dir = Path(user_params['output_path']) / user_params['variable']
    out_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{user_params['input_file']}.png"
    out_path = out_dir / fname

    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to {out_path}")

def plot_individual_blocks(data_dict, user_params):
    """
    Plot each block's 2D slice individually using data from extract_athenak_slice_as_dataframe.
    
    Parameters:
    - data_dict: dict output from extract_athenak_slice_as_dataframe with keys:
        'df_quantities', 'df_extents', 'block_shape', etc.
    - user_params: dict for plotting options (cmap, figsize, xlabel, ylabel, axes_scale)
    """
    df_quantities = data_dict['df_quantities']
    df_extents = data_dict['df_extents']
    block_shape = data_dict['block_shape']
    
    num_blocks = df_extents.shape[0]
    cmap = user_params.get("cmap", "viridis")
    figsize = user_params.get("figsize", (6, 6))
    xlabel = user_params.get("xlabel", "x")
    ylabel = user_params.get("ylabel", "y")
    axes_scale = user_params.get("axes_scale", 1.0)
    vmin, vmax = user_params.get('clim', (None, None))
    
    for block in range(num_blocks):
        block_data = df_quantities.loc[block]
        # Convert Series with MultiIndex (row, col) to 2D array
        slice_2d = block_data.unstack(level=-1).values
        
        extent = df_extents.loc[block]
        img_extent = [
            extent['x_min'] * axes_scale,
            extent['x_max'] * axes_scale,
            extent['y_min'] * axes_scale,
            extent['y_max'] * axes_scale
        ]
        # Determine color limits
        actual_vmin = vmin if vmin is not None else block_data.min()
        actual_vmax = vmax if vmax is not None else block_data.max()
        norm = set_normalization(user_params,actual_vmin,actual_vmax)
        plt.figure(figsize=figsize)
        im = plt.imshow(slice_2d, origin='lower', extent=img_extent, cmap=cmap, norm=norm, aspect='auto')
        plt.colorbar(im)
        plt.title(f"Block {block} slice")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        # Save if output_path given
        if 'output_path' in user_params:
            import os
            os.makedirs(user_params['output_path'], exist_ok=True)
            fname = f"{user_params['variable']}_meshblock_{block}.png"
            plt.savefig(os.path.join(user_params['output_path'], fname), dpi=200, bbox_inches='tight')
            print("figure_saved")
            plt.close()
        else:
            plt.show()

