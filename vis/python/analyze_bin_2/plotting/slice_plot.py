import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
import matplotlib.colors as plt_col
from matplotlib.ticker import LogFormatterMathtext
from pathlib import Path
import os
import numpy as np
from utils.helpers import mask_half, build_split_params
from utils.io_utils import extract_slice_number
from data_processing.slice_data import get_stitched_slice_for_variable
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
    print(global_array)
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

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def _make_norm(user_params, data):
    global_min = np.nanmin(data)
    global_max = np.nanmax(data)
    vmin = user_params["clim"][0] if user_params["clim"][0] is not None else global_min
    vmax = user_params["clim"][1] if user_params["clim"][1] is not None else global_max
    return set_normalization(user_params, vmin, vmax), vmin, vmax


def draw_split_panel(
    ax,
    stitch_data_top, stitch_extent_top, user_params_top,
    stitch_data_bot, stitch_extent_bot, user_params_bot,
    add_colorbars=True, add_panel_title=True
):
    file_number = extract_slice_number(user_params_top["input_file"])

    norm_top, vmin_top, vmax_top = _make_norm(user_params_top, stitch_data_top)
    norm_bot, vmin_bot, vmax_bot = _make_norm(user_params_bot, stitch_data_bot)

    cmap_top = plt.get_cmap(user_params_top["cmap"]).copy()
    cmap_bot = plt.get_cmap(user_params_bot["cmap"]).copy()
    cmap_top.set_bad(alpha=0)
    cmap_bot.set_bad(alpha=0)

    top_plot = mask_half(stitch_data_top, stitch_extent_top, keep="top")
    bot_plot = mask_half(stitch_data_bot, stitch_extent_bot, keep="bottom")

    ax.imshow(
        top_plot,
        origin="lower",
        extent=stitch_extent_top,
        cmap=cmap_top,
        norm=norm_top,
        aspect="auto"
    )

    ax.imshow(
        bot_plot,
        origin="lower",
        extent=stitch_extent_bot,
        cmap=cmap_bot,
        norm=norm_bot,
        aspect="auto"
    )

    # ax.axhline(0, color="white", lw=1.0)

    ax.set_xlabel(user_params_top["xlabel"], labelpad=2.0, fontsize=20, fontweight="bold")
    ax.set_ylabel(user_params_top["ylabel"], labelpad=2.0, fontsize=20, fontweight="bold")
    ax.tick_params(axis="both", labelsize=18)
    ax.set_aspect("equal")
    # for tick in ax.get_xticklabels():
    #     tick.set_fontweight("bold")
    # for tick in ax.get_yticklabels():
    #     tick.set_fontweight("bold")
    if add_panel_title:
        ax.set_title(f"t = {file_number} {user_params_top['time_label']}", fontsize=16)

    # if add_colorbars:
    #     divider = make_axes_locatable(ax)
    #     cax_col = divider.append_axes("right", size="4%", pad=0.03)
    #     cax_col.set_axis_off()

    #     cax_top = cax_col.inset_axes([0.20, 0.55, 0.45, 0.40])
    #     cax_bot = cax_col.inset_axes([0.20, 0.00, 0.45, 0.40])

    #     sm_top = ScalarMappable(norm=norm_top, cmap=cmap_top)
    #     sm_top.set_array([])
    #     cbar_top = ax.figure.colorbar(sm_top, cax=cax_top)
    #     cbar_top.ax.set_title(user_params_top["cmap_label"], loc="left", pad=8)
    #     cbar_top.set_ticks([vmin_top, np.sqrt(vmin_top*vmax_top), vmax_top])
    #     cbar_top.formatter = LogFormatterMathtext()
    #     cbar_top.update_ticks()
    #     # cbar_top.set_ticklabels([f"{vmin_top:.2e}", f"{vmax_top:.2e}"])

    #     sm_bot = ScalarMappable(norm=norm_bot, cmap=cmap_bot)
    #     sm_bot.set_array([])
    #     cbar_bot = ax.figure.colorbar(sm_bot, cax=cax_bot)
    #     cbar_bot.ax.set_title(user_params_bot["cmap_label"], loc="left", pad=8)
    #     cbar_bot.set_ticks([vmin_bot,np.sqrt(vmin_bot * vmax_bot), vmax_bot])
    #     cbar_bot.formatter = LogFormatterMathtext()
    #     cbar_bot.update_ticks()
    #     # cbar_bot.set_ticklabels([f"{vmin_bot:.2e}", f"{vmax_bot:.2e}"])
    if add_colorbars:
        divider = make_axes_locatable(ax)
        cax_col = divider.append_axes("right", size="4%", pad=0.03)
        cax_col.set_axis_off()

        cax_top = cax_col.inset_axes([0.20, 0.55, 0.45, 0.45])
        cax_bot = cax_col.inset_axes([0.20, 0.00, 0.45, 0.45])

        sm_top = ScalarMappable(norm=norm_top, cmap=cmap_top)
        sm_top.set_array([])
        cbar_top = ax.figure.colorbar(sm_top, cax=cax_top)
        cbar_top.set_label(user_params_top["cmap_label"], rotation=90, labelpad=10)
        cbar_top.ax.yaxis.label.set_fontsize(24)
        cbar_top.ax.yaxis.label.set_fontweight("bold")
        cbar_top.set_ticks([vmin_top, np.sqrt(vmin_top * vmax_top), vmax_top])
        for t in cbar_top.ax.get_yticklabels():
            t.set_fontsize(18)
            t.set_fontweight("bold")  
        cbar_top.formatter = LogFormatterMathtext()
        cbar_top.update_ticks()

        sm_bot = ScalarMappable(norm=norm_bot, cmap=cmap_bot)
        sm_bot.set_array([])
        cbar_bot = ax.figure.colorbar(sm_bot, cax=cax_bot)
        cbar_bot.set_label(user_params_bot["cmap_label"], rotation=90, labelpad=10, fontweight="bold", fontsize=24)
        cbar_bot.set_ticks([vmin_bot, np.sqrt(vmin_bot * vmax_bot), vmax_bot])
        for t in cbar_bot.ax.get_yticklabels():
            t.set_fontsize(18)
            t.set_fontweight("bold")
        cbar_bot.formatter = LogFormatterMathtext()
        cbar_bot.update_ticks()

def plot_split_stitched_data(
    stitch_data_top, stitch_extent_top, user_params_top,
    stitch_data_bot, stitch_extent_bot, user_params_bot
):
    fig, ax = plt.subplots(figsize=user_params_top.get("figsize", (12, 6)))

    draw_split_panel(
        ax,
        stitch_data_top, stitch_extent_top, user_params_top,
        stitch_data_bot, stitch_extent_bot, user_params_bot,
        add_colorbars=True, add_panel_title=True
    )

    out_dir = Path(user_params_top["output_path"]) / f"{user_params_top['variable']}_{user_params_bot['variable']}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{user_params_top['input_file']}.png"
    out_path = out_dir / fname

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {out_path}")
    
def plot_combined_split_pairs(user_params, save_individual=True):
    file_number = extract_slice_number(user_params["input_file"])
    pairs = user_params["variable"]
    # ncolumns = len(pairs)
    nrows = len(pairs)
    # If only one pair, just make the individual plot and skip the combined figure
    if nrows == 1:
        top_params = build_split_params(user_params, 0, 0)
        bot_params = build_split_params(user_params, 0, 1)

        top_data, top_extent, _ = get_stitched_slice_for_variable(top_params)
        bot_data, bot_extent, _ = get_stitched_slice_for_variable(bot_params)

        plot_split_stitched_data(
            top_data, top_extent, top_params,
            bot_data, bot_extent, bot_params
        )
        return
    figsize_ind = user_params.get("figsize")
    figsize_combined = (figsize_ind[0], 10)
    # figsize_combined = (21,figsize_ind[1])
    # fig, axes = plt.subplots(
    #     nrows=1,
    #     ncols=ncolumns,
    #     figsize = figsize_combined,
    #     squeeze=False
    # )
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=1,
        figsize = figsize_combined,
        squeeze=False
    )
    axes = axes.ravel()
    fig.suptitle(f"t = {file_number} {user_params['time_label']}", fontsize=18, y=0.995)

    for pair_idx, ax in enumerate(axes):
        top_params = build_split_params(user_params, pair_idx, 0)
        bot_params = build_split_params(user_params, pair_idx, 1)

        top_data, top_extent,_ = get_stitched_slice_for_variable(top_params)
        bot_data, bot_extent,_ = get_stitched_slice_for_variable(bot_params)

        # save individual panel too
        if save_individual:
            plot_split_stitched_data(
                top_data, top_extent, top_params,
                bot_data, bot_extent, bot_params
            )

        # draw into the combined multi-panel figure
        draw_split_panel(
            ax,
            top_data, top_extent, top_params,
            bot_data, bot_extent, bot_params,
            add_colorbars=True, add_panel_title= False
        )

        # ax.set_title(
        #     f"{pairs[pair_idx][0]} / {pairs[pair_idx][1]}",
        #     loc="left",
        #     fontsize=14
        # )

    # plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.tight_layout()

    out_dir = Path(user_params["output_path"]) / "combined_pairs"
    out_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{user_params['input_file']}_combined.png"
    out_path = out_dir / fname

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {out_path}")

# def plot_split_stitched_data(
#     stitch_data_top, stitch_extent_top, user_params_top,
#     stitch_data_bot, stitch_extent_bot, user_params_bot
# ):
#     file_number = extract_slice_number(user_params_top['input_file'])
#     global_min_top = np.nanmin(stitch_data_top)
#     global_max_top = np.nanmax(stitch_data_top)
#     global_min_bot = np.nanmin(stitch_data_bot)
#     global_max_bot = np.nanmax(stitch_data_bot)

#     vmin_top = user_params_top["clim"][0] if user_params_top["clim"][0] is not None else global_min_top
#     vmax_top = user_params_top["clim"][1] if user_params_top["clim"][1] is not None else global_max_top
#     vmin_bot = user_params_bot["clim"][0] if user_params_bot["clim"][0] is not None else global_min_bot
#     vmax_bot = user_params_bot["clim"][1] if user_params_bot["clim"][1] is not None else global_max_bot

#     norm_top = set_normalization(user_params_top, vmin_top, vmax_top)
#     norm_bot = set_normalization(user_params_bot, vmin_bot, vmax_bot)
    
#     cmap_top = plt.get_cmap(user_params_top["cmap"]).copy()
#     cmap_bot = plt.get_cmap(user_params_bot["cmap"]).copy()
#     cmap_top.set_bad(alpha=0)
#     cmap_bot.set_bad(alpha=0)

#     top_plot = mask_half(stitch_data_top, stitch_extent_top, keep="top")
#     bot_plot = mask_half(stitch_data_bot, stitch_extent_bot, keep="bottom")
#     fig, ax = plt.subplots(figsize=user_params_top.get("figsize", (12, 6)))

#     im_top = ax.imshow(
#         top_plot,
#         origin="lower",
#         extent=stitch_extent_top,
#         cmap=cmap_top,
#         norm=norm_top,
#         aspect="auto"
#     )

#     im_bot = ax.imshow(
#         bot_plot,
#         origin="lower",
#         extent=stitch_extent_bot,
#         cmap=cmap_bot,
#         norm=norm_bot,
#         aspect="auto"
#     )
#     ax.axhline(0, color="white", lw=1.0)
#     ax.set_xlabel(user_params_top["xlabel"], labelpad=2.0, fontsize=14)
#     ax.set_ylabel(user_params_top["ylabel"], labelpad=2.0, fontsize=14)
#     ax.set_aspect("equal")
#     ax.set_title(f"t = {file_number} {user_params_top['time_label']}", fontsize=16)

#     if "title" in user_params_top:
#         ax.set_title(user_params_top["title"])

#     # One shared right-side column
#     divider = make_axes_locatable(ax)
#     cax_col = divider.append_axes("right", size="4%", pad=0.03)
#     cax_col.set_axis_off()
#     cax_top = cax_col.inset_axes([0.20, 0.55, 0.45, 0.40])
#     cax_bot = cax_col.inset_axes([0.20, 0.00, 0.45, 0.40])

#     sm_top = ScalarMappable(norm=norm_top, cmap=user_params_top["cmap"])
#     sm_top.set_array([])
#     cbar_top = fig.colorbar(sm_top, cax=cax_top)
#     cbar_top.ax.set_title(user_params_top["cmap_label"], loc="left", pad=8)

#     sm_bot = ScalarMappable(norm=norm_bot, cmap=user_params_bot["cmap"])
#     sm_bot.set_array([])
#     cbar_bot = fig.colorbar(sm_bot, cax=cax_bot)
#     cbar_bot.ax.set_title(user_params_bot["cmap_label"], loc="left", pad=8)

#     out_dir = Path(user_params_top["output_path"]) / f"{user_params_top['variable']}_{user_params_bot['variable']}"
#     out_dir.mkdir(parents=True, exist_ok=True)
#     fname = f"{user_params_top['input_file']}.png"
#     out_path = out_dir / fname

#     fig.savefig(out_path, dpi=300, bbox_inches="tight")
#     plt.close(fig)
#     print(f"Saved plot to {out_path}")