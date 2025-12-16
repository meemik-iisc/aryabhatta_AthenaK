import os
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
from PIL import Image, ImageDraw, ImageFont

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

# Prompt for folder containing binary files
input_dir = Path(input("Enter path to binary folder: "))
if not input_dir.is_dir():
    raise NotADirectoryError(f"{input_dir} is not a directory")
#Timestep interval
dt=1.0
#Plotting variables
dens_var={
    'cmap':"coolwarm",
    'norm':None,
    'vmin':None,
    'vmax':None,
    'x1min':None,
    'x1max':None,
    'x2min':None,
    'x2max':None,
}
pres_var={
    'cmap':"viridis",
    'norm':None,
    'vmin':None,
    'vmax':None,
    'x1min':None,
    'x1max':None,
    'x2min':None,
    'x2max':None,
}
velr_var={
    'cmap':"seismic",
    'norm':None,
    'vmin':None,
    'vmax':None,
    'x1min':None,
    'x1max':None,
    'x2min':None,
    'x2max':None,
}
# vely_var={
#     'cmap':"seismic",
#     'norm':None,
#     'vmin':None
#     'vmax':None
#     'x1min':None,
#     'x1max':None,
#     'x2min':None,
#     'x2max':None,
# }
temp_var={
    'cmap':"hot",
    'norm':None,
    'vmin':None,
    'vmax':None,
    'x1min':None,
    'x1max':None,
    'x2min':None,
    'x2max':None,
}

# Prepare output directories
out_root = input_dir.parent / (input_dir.name + '_outputs')
folders = {
    'dens': out_root / 'dens',
    'pres': out_root / 'pres',
    'velr': out_root / 'velr',
    # 'vely': out_root / 'vely',
    'temp': out_root / 'temp',
    'data': out_root / 'data'
}
for folder in folders.values():
    ensure_dir(folder)

# Locate plot_slice.py (must accept args: input_file, quantity, output_file)
plotter = Path(__file__).parent / 'plot_slice_2.py'
if not plotter.exists():
    raise FileNotFoundError(f"Could not find plot_slice.py at {plotter}")

# Find all .bin files in the input directory
bin_files = sorted(input_dir.glob('*.bin'))
if not bin_files:
    print("No .bin files found in the specified directory.")
    exit(0)

# Define the quantities and their output folders
quantities = {
    'dens': folders['dens'],
    'pgas': folders['pres'],
    'velr': folders['velr'],
    # 'vely': folders['vely'],
    'temp': folders['temp']
}

def plot_data(bf,qty,out_file,qtyvar,dim='z'):
    cmd=['python', str(plotter), str(bf), qty, str(out_file),
            f'--dimension={dim}', f"--cmap={qtyvar['cmap']}", "--notex"]
    if qtyvar['norm'] is not None:
        cmd.append( f"--norm={qtyvar['norm']}")
    if qtyvar['vmin'] is not None:
        cmd.append( f"--vmin={qtyvar['vmin']}")
    if qtyvar['vmax'] is not None:
        cmd.append( f"--vmax={qtyvar['vmax']}")
    if qtyvar['x1min'] is not None:
        cmd.append( f"--x1_min={qtyvar['x1min']}")
    if qtyvar['x1max'] is not None:
        cmd.append( f"--x1_max={qtyvar['x1max']}")
    if qtyvar['x2min'] is not None:
        cmd.append( f"--x2_min={qtyvar['x2min']}")
    if qtyvar['x2max'] is not None:
        cmd.append( f"--x2_max={qtyvar['x2max']}")
    subprocess.run(cmd, check=True)


# Process each .bin file and generate plots
for idx,bf in enumerate(bin_files):
    timestep=idx*dt
    basename = bf.stem
    # Generate individual plots for density

    # Plot density
    dens_out = folders['dens'] / f"{basename}_dens.png"
    plot_data(bf,'dens',dens_out,dens_var)

    # Plot pressure (derived)
    pres_out = folders['pres'] / f"{basename}_pgas.png"
    plot_data(bf,'eint',pres_out,pres_var)
    
    

    # Plot radial velocity
    velr_out = folders['velr'] / f"{basename}_velr.png"
    plot_data(bf,'derived:velr',velr_out,velr_var)

    # # Plot velocity y
    # vely_out = folders['vely'] / f"{basename}_vely.png"
    # plot_data(bf,'vely',vely_out,vely_var)
    
    #Plot temperature
    temp_out = folders['temp'] / f"{basename}_temp.png"
    plot_data(bf,'derived:T',temp_out,temp_var)


    # Create combined 4x1 subplot: dens, pgas, velr, temp
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    plot_order = ['dens', 'pgas', 'velr', 'temp']
    titles = ['Density', 'Pressure', 'Radial Velocity', 'Temperature']
    for ax, key, title in zip(axes.flatten(), plot_order, titles):
        img_path = folders[key if key != 'pgas' else 'pres'] / f"{basename}_{key}.png"
        if img_path.exists():
            img = plt.imread(str(img_path))
            ax.imshow(img)
        ax.set_title(title)
        ax.axis('off')
    fig.suptitle(f"t={timestep:.2f}"+r'$t_{ff}$')

    combo_file = folders['data'] / f"{basename}_combined.png"
    fig.tight_layout()
    fig.savefig(combo_file, dpi=300)
    plt.close(fig)

# Function to create movie from image sequence
# def create_movie(img_folder, pattern, output_path, fps=2):
#     images = []
#     for img_file in sorted(img_folder.glob(pattern)):
#         images.append(imageio.imread(str(img_file)))
#     if images:
#         writer = imageio.get_writer(str(output_path), fps=fps)
#         for img in images:
#             writer.append_data(img)
#         writer.close()

# def create_movie(img_folder, pattern, output_path, fps=2, t_step=0.1, t_unit='Myr'):
#     images = []
#     font = ImageFont.load_default()  # For better look, use truetype font
    
#     sorted_imgs = sorted(img_folder.glob(pattern))
    
#     for i, img_file in enumerate(sorted_imgs):
#         t_val = i * t_step
#         timestamp = f"t = {t_val:.1f} {t_unit}"
        
#         # Open image and convert to RGB
#         img = Image.open(img_file).convert("RGB")
#         draw = ImageDraw.Draw(img)
        
#         # Draw time label
#         draw.text((10, 10), timestamp, font=font, fill=(0, 0, 0))
        
#         # Convert PIL image to numpy array
#         img_np = np.array(img)
#         images.append(img_np)

#     # Write MP4 movie
#     if images:
#         imageio.mimsave(output_path, images, fps=fps)

# # Generate movies for each quantity
# for qty, folder in quantities.items():
#     movie_path = out_root / f"{qty}.mp4"
#     create_movie(folder, f"*_{qty}.png", movie_path)

# # Generate movie for combined data plots
# combo_folder = folders['data']
# combo_movie = out_root / 'combined.mp4'
# create_movie(combo_folder, "*_combined.png", combo_movie)

print(f"All slices, combined plots, and movies saved under {out_root}")
