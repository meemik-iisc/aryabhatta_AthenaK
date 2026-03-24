import os
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont

# ============================================================================
# CONFIGURATION - Easy to modify
# ============================================================================
config = {
    # File processing
    'dt': 1.0,                    # Timestep interval (code units)
    'time_label': 'Myr',
    'plot_dimension': 'z',
    
    # Figure sizes
    # 'fig_size_single': (22.0, 24),  # Combined subplot
    'fig_size_single': (13.0, 24),  # Combined subplot
    # 'fig_size_single': (7.0, 15.0),  # Combined subplot
    
    
    # Plot variables (name → settings)
    'plot_vars': {
        'dens': {
            'label': 'Density',
            'folder_key': 'dens',
            'quantity': 'dens',
            'cmap': 'PiYG',
            'norm': 'log',
            'vmin': 1.0e-5, 
            'vmax': 1.0e-2
            # 'vmin': None, 
            # 'vmax': None
        },
        'pres': {
            'label': 'Pressure',
            'folder_key': 'pres',
            'quantity': 'eint',      # Derived from internal energy
            'cmap': 'viridis',
            'norm': 'log',
            'vmin': 1e-7, 
            'vmax': 1e-1
        },
        'velx': {
            'label': 'X Velocity',
            'folder_key': 'velx',
            'quantity': 'velx',
            'cmap': 'seismic',
            'norm': None,            # Linear
            'vmin': -1e4, 
            'vmax': 1e4
        },
        'temp': {
            'label': 'Temperature',
            'folder_key': 'temp',
            'quantity': 'derived:T',
            'cmap': 'coolwarm',
            'norm': 'log',
            'vmin': 1e4, 
            'vmax': 1e9
        },
        't_cool': {
            'label': 'Cooling Time',
            'folder_key': 't_cool',
            'quantity': 'derived:t_cool',
            'cmap': 'turbo',
            'norm': 'log',
            'vmin': 1, 
            'vmax': 1e6
        },
        'tracer': {
            'label': 'Outflow Tracer',
            'folder_key': 'tracer',
            'quantity': 's_00',
            'cmap': 'bwr',
            'norm': None,
            'vmin': 0.0, 
            'vmax': 1.0
            # 'vmin': None, 
            # 'vmax': None
        }
    },
    
    # Combined plot order (top to bottom)
    'plot_order': ['dens', 'pres', 'velx', 'temp', 't_cool', 'tracer']
    # 'plot_order': ['dens', 'pres', 'velx', 'temp', 't_cool']
    
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def ensure_dir(path: Path):
    """Create directory if it doesn't exist."""
    path.mkdir(parents=True, exist_ok=True)

def get_plot_cmd(plotter_path: Path, bin_file: Path, quantity: str, out_file: Path, 
                 var_config: dict, dimension: str = 'z') -> list:
    """Build subprocess command for plot_slice_2.py."""
    cmd = [
        'python', str(plotter_path), str(bin_file), quantity, str(out_file),
        f'--dimension={dimension}', f"--cmap={var_config['cmap']}", "--notex"
    ]
    
    if var_config['norm']:
        cmd.append(f"--norm={var_config['norm']}")
    if var_config['vmin'] is not None:
        cmd.append(f"--vmin={var_config['vmin']}")
    if var_config['vmax'] is not None:
        cmd.append(f"--vmax={var_config['vmax']}")
    
    return cmd

def plot_single_quantity(bin_file: Path, var_name: str, out_root: Path, 
                        plotter: Path, config: dict) -> Path:
    """Generate single plot for one quantity."""
    var_config  = config['plot_vars'][var_name]
    out_dir     = out_root / var_config['folder_key']
    ensure_dir(out_dir) 
    out_file    =out_dir / f"{bin_file.stem}_{var_name}.png"
    
    cmd = get_plot_cmd(plotter, bin_file, var_config['quantity'], out_file, 
                      var_config, config['plot_dimension'])
    
    subprocess.run(cmd, check=True)
    print(f"{var_name}: {out_file.name}")
    return out_file

def create_combined_subplot(out_root: Path, basename: str, timestep: float, 
                          config: dict) -> Path:
    """Create vertical stack of all plots."""
    out_dir = out_root / "combined"
    ensure_dir(out_dir)
    fig, axes = plt.subplots(len(config['plot_order']), 1, figsize=config['fig_size_single'])
    if len(config['plot_order']) == 1:
        axes = [axes]
    
    titles = [config['plot_vars'][k]['label'] for k in config['plot_order']]
    
    for ax, var_name, title in zip(axes.flatten(), config['plot_order'], titles):
        img_path = (out_dir.parent / config['plot_vars'][var_name]['folder_key'] 
                   / f"{basename}_{var_name}.png")
        if img_path.exists():
            img = plt.imread(str(img_path))
            ax.imshow(img)
        ax.set_title(title, fontsize=12)
        ax.axis('off')
    
    fig.suptitle(f"t={timestep:.2f} {config['time_label']}", fontsize=18)
    # fig.tight_layout()
    # Reserve top space FIRST
    fig.subplots_adjust(top=0.98)

    # THEN tight_layout, avoiding suptitle region
    fig.tight_layout(rect=[0, 0, 1, 0.98])  # rect=[left, bottom, right, top]
    # Position suptitle manually ABOVE tight_layout area
    # fig.tight_layout(rect=[0, 0, 1, 0.90])  # rect=[left, bottom, right, top]
    
    combo_file = out_dir / f"{basename}_combined.png"
    fig.savefig(combo_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Combined: {combo_file.name}")
    return combo_file

# Replace the matplotlib PdfPages function with this:

def create_pdf_from_combined(out_root:Path, config: dict):
    """Convert *_combined.png → sorted PDF using img2pdf (lossless, fast)."""
    try:
        import img2pdf
    except ImportError:
        print("Install img2pdf: pip install img2pdf")
        return
    
    combined_dir = out_root / "combined"
    pdf_path = combined_dir / "output.pdf"
    combo_files = sorted(combined_dir.glob("*_combined.png"))  # Natural order
    
    if not combo_files:
        print("No combined PNGs found!")
        return
    
    print(f"Creating PDF ({len(combo_files)} pages) with img2pdf...")
    
    # Convert directly: PNGs → PDF bytes → write file
    pdf_bytes = img2pdf.convert([str(f) for f in combo_files])
    
    with open(pdf_path, "wb") as f:
        f.write(pdf_bytes)
    
    print(f"✓ PDF saved: {pdf_path} ({len(combo_files)} pages)")


# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def main():
    # 1. Get input directory
    input_dir = Path(input("Enter path to binary folder: "))
    if not input_dir.is_dir():
        raise NotADirectoryError(f"{input_dir} is not a directory")
    
    # 2. Setup output directories
    out_root = input_dir.parent / f"{input_dir.name}_outputs"
        
    # 3. Find plotter
    plotter = Path(__file__).parent / 'plot_slice_2.py'
    if not plotter.exists():
        raise FileNotFoundError(f"Missing {plotter}")
    
    # 4. Process all .bin files
    bin_files = sorted(input_dir.glob('*.bin'))
    if not bin_files:
        print("No .bin files found!")
        return
    
    print(f"Processing {len(bin_files)} files...")
    
    for idx, bin_file in enumerate(bin_files):
        timestep = idx * config['dt']
        basename = bin_file.stem
        
        # Generate individual plots
        for var_name in config['plot_order']:
            plot_single_quantity(bin_file, var_name, out_root, 
                               plotter, config)
        
        # Create combined plot
        create_combined_subplot(out_root, basename, timestep, config)
    
    print(f"\nAll done! Output saved to: {out_root}")
    
    #Create combined PDF
    create_pdf_from_combined(out_root, config)

if __name__ == "__main__":
    main()
