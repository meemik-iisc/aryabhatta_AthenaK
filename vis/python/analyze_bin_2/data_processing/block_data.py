import struct
import numpy as np
import pandas as pd
from utils.constants import mu,mp_cgs,kB_cgs,G_cgs,s_Myr,M_bh_cgs,gamma,Msolar_per_kpc_sq_per_yr,\
length_cgs,mass_cgs,time_cgs,pres_cgs,rho_cgs,v_cgs,\
Temp_norm,Velr_scale, dens_scale, pres_scale, tcool_threshold_Myr, t_in_threshold_Myr
from utils.helpers import cool_lambda


def extract_athenak_3D_block(user_params):
    """
    Extract 3D block data of requested variable from AthenaK binary file.

    Returns dict with:
      - df_extents: DataFrame with columns ['x_min','x_max','y_min','y_max','z_min','z_max']
      - blocks: list of per-block dicts:
            * data  : 3D numpy array (nz, ny, nx)
            * x, y, z : 1D coordinate arrays (cell centers)
            * extent  : (x_min,x_max,y_min,y_max,z_min,z_max)
            * ijk     : (block_i, block_j, block_k, block_level)
      - num_blocks: int
      - block_shape: (nz, ny, nx)
    """
    file_path   = user_params["input_folder"].rstrip('/') + '/' + user_params["input_file"]
    variable    = user_params["variable"]
    if variable.startswith('derived:'):
        variable    = variable.split(":", 1)[1].strip()
    else:
        var_name    = variable
        
        
    with open(file_path, "rb") as f:
        line    = f.readline().decode('ascii')
        if line != 'Athena binary output version=1.1\n':
            raise RuntimeError('Unrecognized AthenaK binary file format.')
        for _ in range(3):
            next(f)
        location_size       = int(f.readline().decode('ascii')[19:])
        variable_size       = int(f.readline().decode('ascii')[19:])
        next(f)
        variable_names_base = f.readline().decode('ascii')[12:].split()
        header_offset       = int(f.readline().decode('ascii')[16:])

        

        if var_name not in variable_names_base:
            raise RuntimeError(f'Variable "{var_name}" not found in file variables.')
        var_ind             = variable_names_base.index(var_name)

        location_fmt        = 'f' if location_size == 4 else 'd'
        variable_fmt        = 'f' if variable_size == 4 else 'd'

        input_data          = {}
        start_data          = f.tell() + header_offset
        while f.tell() < start_data:
            line            = f.readline().decode('ascii')
            if line.startswith('#'):
                continue
            if line.startswith('<'):
                section     = line[1:-2]
                input_data[section] = {}
                continue
            if '=' in line:
                k, v        = line.split('=', 1)
                input_data[section][k.strip()] = v.split('#')[0].strip()

        num_ghost           = int(input_data['mesh']['nghost'])
        num_variables_base  = len(variable_names_base)
        file_size           = f.seek(0, 2)
        f.seek(start_data, 0)

        quantities_list     = []
        extents_list        = []
        blocks              = []
        num_blocks_used     = 0
        block_nx_ref = block_ny_ref = block_nz_ref = None

        while f.tell() < file_size:
            # indices and AMR level
            block_indices           = np.array(struct.unpack('@6i', f.read(24))) - num_ghost
            block_i, block_j, block_k, block_level = struct.unpack('@4i', f.read(16))

            block_nx = block_indices[1] - block_indices[0] + 1
            block_ny = block_indices[3] - block_indices[2] + 1
            block_nz = block_indices[5] - block_indices[4] + 1
            cells_per_block = block_nz * block_ny * block_nx
            block_cell_fmt = '=' + str(cells_per_block) + variable_fmt
            variable_data_size = cells_per_block * variable_size

            if block_nx_ref is None:
                block_nx_ref = block_nx
                block_ny_ref = block_ny
                block_nz_ref = block_nz

            # spatial extents
            block_lims = struct.unpack('=6' + location_fmt, f.read(6 * location_size))
            x_min, x_max, y_min, y_max, z_min, z_max = block_lims
            extents_list.append((x_min, x_max, y_min, y_max, z_min, z_max))

            # read this variable's data
            cell_data_start = f.tell()
            f.seek(cell_data_start + var_ind * variable_data_size, 0)
            cell_data = np.array(struct.unpack(block_cell_fmt, f.read(variable_data_size)))
            cell_data = cell_data.reshape((block_nz, block_ny, block_nx))

            if user_params['variable']=="dens":
                cell_data = dens_scale*cell_data
            if user_params['variable']=="eint":
                cell_data = pres_scale*cell_data
            # store per-block numpy representation for fast later use
            x_coords = np.linspace(x_min, x_max, block_nx)
            y_coords = np.linspace(y_min, y_max, block_ny)
            z_coords = np.linspace(z_min, z_max, block_nz)
            blocks.append(
                dict(
                    data    = cell_data,
                    x       = x_coords,
                    y       = y_coords,
                    z       = z_coords,
                    extent  = (x_min, x_max, y_min, y_max, z_min, z_max),
                    ijk     = (block_i, block_j, block_k, block_level),
                )
            )
            num_blocks_used += 1

            # move to next block
            f.seek(cell_data_start + num_variables_base * variable_data_size, 0)

        # extents DataFrame
        df_extents = pd.DataFrame(
            extents_list,
            columns=['x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max']
        )
        return {
            "df_extents"    : df_extents,
            "blocks"        : blocks,
            "num_blocks"    : num_blocks_used,
            "block_shape"   : (block_nz_ref, block_ny_ref, block_nx_ref)  # (nz, ny, nx)
        }
        
def extract_temp_data(user_params):
    rho_params              = user_params.copy()
    rho_params['variable']  = 'dens'
    rho_data_dict           = extract_athenak_3D_block(rho_params)
    rho_blocks              = rho_data_dict['blocks']
    
    # Extract internal energy (eint) as 3D blocks
    eint_params             = user_params.copy()
    eint_params['variable'] = 'eint'
    eint_data_dict          = extract_athenak_3D_block(eint_params)
    eint_blocks             = eint_data_dict['blocks']
    
    # Compute temperature for each block: T = (eint/rho) * mu*mp/(kB*(gamma-1)) * (L/t)^2 / Temp_norm
    temp_blocks = []
    
    for block_idx, (rho_block, eint_block) in enumerate(zip(rho_blocks, eint_blocks)):
        rho_data    = rho_block['data']       # (nz, ny, nx)
        eint_data   = eint_block['data']     # (nz, ny, nx)
        
        # Avoid division by zero
        rho_safe    = np.where(rho_data > 0, rho_data, np.nan)
        
        # Temperature calculation
        temp_data   = (eint_data / rho_safe) * \
                    (mu * mp_cgs / (kB_cgs * (gamma - 1))) * \
                    ((length_cgs / time_cgs) ** 2) / Temp_norm
        
        # Create temperature block dict (same structure as input blocks)
        temp_block      = dict(
            data    = temp_data,                    # (nz, ny, nx)
            x       = rho_block['x'],                  # 1D coordinate array
            y       = rho_block['y'],                  # 1D coordinate array
            z       = rho_block['z'],                  # 1D coordinate array
            extent  = rho_block['extent'],        # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk     = rho_block['ijk'],              # (block_i, block_j, block_k, block_level)
        )
        temp_blocks.append(temp_block)
    
    return {
        "df_extents"    : eint_data_dict['df_extents'],
        "blocks"        : temp_blocks,
        "num_blocks"    : len(temp_blocks),
        "block_shape"   : eint_data_dict['block_shape']  # (nz, ny, nx)
    }

def extract_velr_data(user_params):
    vel_params                  = user_params.copy()
    vel_params.update(variable="velx")
    velx_data_dict              = extract_athenak_3D_block(vel_params)
    velx_blocks                 = velx_data_dict['blocks']
    vel_params.update(variable="vely")
    vely_blocks                 = extract_athenak_3D_block(vel_params)['blocks']
    vel_params.update(variable="velz")
    velz_blocks                 = extract_athenak_3D_block(vel_params)['blocks']
    
    velr_blocks = []
    
    for block_idx, (velx_block, vely_block, velz_block) in enumerate(zip(velx_blocks, vely_blocks, velz_blocks)):
        velx_data = velx_block['data']       # (nz, ny, nx)
        vely_data = vely_block['data']     # (nz, ny, nx)
        velz_data = velz_block['data']
                
        # Temperature calculation
        velr_data = np.sqrt(np.square(velx_data)+np.square(vely_data)+np.square(velz_data))*Velr_scale
        
        # Create velocity block dict (same structure as input blocks)
        velr_block      = dict(
            data    = velr_data,                        # (nz, ny, nx)
            x       = velx_block['x'],                  # 1D coordinate array
            y       = velx_block['y'],                  # 1D coordinate array
            z       = velx_block['z'],                  # 1D coordinate array
            extent  = velx_block['extent'],             # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk     = velx_block['ijk'],                # (block_i, block_j, block_k, block_level)
        )
        velr_blocks.append(velr_block)
    
    return {
        "df_extents"    : velx_data_dict['df_extents'],
        "blocks"        : velr_blocks,
        "num_blocks"    : len(velr_blocks),
        "block_shape"   : velx_data_dict['block_shape']  # (nz, ny, nx)
    }
def extract_cooling_rate_data(user_params):
    eint_params                 = user_params.copy()
    eint_params.update(variable = "eint")
    eint_data_dict              = extract_athenak_3D_block(eint_params)
    eint_blocks                 = eint_data_dict['blocks']
    rho_params                  = user_params.copy()
    rho_params.update(variable  = "dens")
    rho_blocks                  = extract_athenak_3D_block(rho_params)['blocks']
    
    cooling_rate_blocks = []
    
    for block_idx, (eint_block, rho_block) in enumerate(zip(eint_blocks, rho_blocks)):
        eint_data       = eint_block['data']
        rho_data        = rho_block['data']
        
        #Convert to CGS
        pres_data_cgs   = eint_data*(gamma-1)*pres_cgs
        rho_data_cgs    = rho_data*rho_cgs
        
        # Temperature calculation
        temp_data_cgs   = (pres_data_cgs*mu*mp_cgs)/(rho_data_cgs*kB_cgs)
        #Cooling Rate calculation
        lambda_cgs      = cool_lambda(temp_data_cgs)
        #Calculate cooling time
        cooling_rate_cgs = ((rho_data_cgs**2)*lambda_cgs)/((mu**2)*(mp_cgs**2))
        # Create velocity block dict (same structure as input blocks)
        cooling_rate_block = dict(
            data    = cooling_rate_cgs,                 # (nz, ny, nx)
            x       = rho_block['x'],                   # 1D coordinate array
            y       = rho_block['y'],                   # 1D coordinate array
            z       = rho_block['z'],                   # 1D coordinate array
            extent  = rho_block['extent'],              # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk     = rho_block['ijk'],                 # (block_i, block_j, block_k, block_level)
        )
        cooling_rate_blocks.append(cooling_rate_block)
    
    return {
        "df_extents"    : eint_data_dict['df_extents'],
        "blocks"        : cooling_rate_blocks,
        "num_blocks"    : len(cooling_rate_blocks),
        "block_shape"   : eint_data_dict['block_shape']  # (nz, ny, nx)
    }
    
def extract_t_cool_data(user_params):
    eint_params                 = user_params.copy()
    eint_params.update(variable = "eint")
    eint_data_dict              = extract_athenak_3D_block(eint_params)
    eint_blocks                 = eint_data_dict['blocks']
    cooling_rate_blocks         = extract_cooling_rate_data(user_params)['blocks']
    
    tcool_blocks = []
    
    for block_idx, (eint_block, cooling_rate_block) in enumerate(zip(eint_blocks, cooling_rate_blocks)):
        eint_data               = eint_block['data']
        cooling_rate_data_cgs   = cooling_rate_block['data']
        
        #Convert to CGS
        pres_data_cgs           = eint_data*(gamma-1)*pres_cgs
        #Calculate cooling time
        t_cool_cgs              = (gamma*pres_data_cgs)/((gamma-1)*cooling_rate_data_cgs)
        t_cool_Myr              = t_cool_cgs/s_Myr
        # ===== MASK EXTREME VALUES =====
        bad_mask                = t_cool_Myr > tcool_threshold_Myr
        if np.any(bad_mask):
            # print(f"Block {block_idx}: Masking {np.sum(bad_mask)} points with t_cool > {tcool_threshold_Myr} Myr")
            t_cool_Myr[bad_mask] = np.nan  # Set to NaN
        # ================================
        # Create velocity block dict (same structure as input blocks)
        tcool_block     = dict(
            data    = t_cool_Myr,                    # (nz, ny, nx)
            x       = eint_block['x'],                  # 1D coordinate array
            y       = eint_block['y'],                  # 1D coordinate array
            z       = eint_block['z'],                  # 1D coordinate array
            extent  = eint_block['extent'],        # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk     = eint_block['ijk'],              # (block_i, block_j, block_k, block_level)
        )
        tcool_blocks.append(tcool_block)
    
    return {
        "df_extents"    : eint_data_dict['df_extents'],
        "blocks"        : tcool_blocks,
        "num_blocks"    : len(tcool_blocks),
        "block_shape"   : eint_data_dict['block_shape']  # (nz, ny, nx)
    }
def extract_t_free_fall_data(user_params):
    rho_params                  = user_params.copy()
    rho_params.update(variable  = "dens")
    rho_data_dict               = extract_athenak_3D_block(rho_params)
    rho_blocks                  = rho_data_dict['blocks']
    
    nz, ny, nx  = rho_data_dict['block_shape']
    
    x_center    = user_params.get('center_x', 0.0)
    y_center    = user_params.get('center_y', 0.0)
    z_center    = user_params.get('center_z', 0.0)
    
    t_ff_blocks = []
    # Iterate through each block
    for block_idx, rho_block in enumerate(rho_blocks):
        
        rho_data_3d     = rho_block['data']
        x_coords        = rho_block['x']
        y_coords        = rho_block['y']
        z_coords        = rho_block['z']
        x_min, x_max, y_min, y_max, z_min, z_max = rho_block['extent']
        
        # Create 3D coordinate grids
        xx, yy, zz      = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
        xx              = np.transpose(xx, (2, 1, 0))
        yy              = np.transpose(yy, (2, 1, 0))
        zz              = np.transpose(zz, (2, 1, 0))
        
        # Compute radial distance
        rr              = np.sqrt((xx - x_center)**2 + (yy - y_center)**2 + (zz - z_center)**2)
        r_cgs           = rr*length_cgs
        t_free_fall_cgs = np.pi*np.sqrt(r_cgs**3/(2*G_cgs*M_bh_cgs))
        t_free_fall_Myr = t_free_fall_cgs/s_Myr
        t_ff_block      = dict(
            data    = t_free_fall_Myr,                    # (nz, ny, nx)
            x       = rho_block['x'],                  # 1D coordinate array
            y       = rho_block['y'],                  # 1D coordinate array
            z       = rho_block['z'],                  # 1D coordinate array
            extent  = rho_block['extent'],        # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk     = rho_block['ijk'],              # (block_i, block_j, block_k, block_level)
        )
        t_ff_blocks.append(t_ff_block)
    return {
        "df_extents"    : rho_data_dict['df_extents'],
        "blocks"        : t_ff_blocks,
        "num_blocks"    : len(t_ff_blocks),
        "block_shape"   : rho_data_dict['block_shape']  # (nz, ny, nx)
    }

def extract_t_in_data(user_params):
    velr_data_dict  = extract_velr_data(user_params)
    velr_blocks     = velr_data_dict['blocks']
    
    nz, ny, nx      = velr_data_dict['block_shape']
    
    x_center        = user_params.get('center_x', 0.0)
    y_center        = user_params.get('center_y', 0.0)
    z_center        = user_params.get('center_z', 0.0)
    
    t_in_blocks = []
    # Iterate through each block
    for block_idx, velr_block in enumerate(velr_blocks):
        
        velr_data_3d    = velr_block['data']
        velr_data_cgs   = velr_data_3d*v_cgs
        x_coords        = velr_block['x']
        y_coords        = velr_block['y']
        z_coords        = velr_block['z']
        x_min, x_max, y_min, y_max, z_min, z_max = velr_block['extent']
        
        # Create 3D coordinate grids
        xx, yy, zz      = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
        xx              = np.transpose(xx, (2, 1, 0))
        yy              = np.transpose(yy, (2, 1, 0))
        zz              = np.transpose(zz, (2, 1, 0))
        
        # Compute radial distance
        rr              = np.sqrt((xx - x_center)**2 + (yy - y_center)**2 + (zz - z_center)**2)
        r_cgs           = rr*length_cgs
        t_in_cgs        = r_cgs/velr_data_cgs
        t_in_Myr        = t_in_cgs/s_Myr
        # ===== MASK INFINITY VALUES =====
        mask = (t_in_Myr > t_in_threshold_Myr) | np.isinf(t_in_Myr) | np.isnan(t_in_Myr)

        if np.any(mask):
            # print(f"Block {block_idx}: Masking {np.sum(bad_mask)} points with t_cool > {tcool_threshold_Myr} Myr")
            t_in_Myr[mask] = t_in_threshold_Myr  # Set to NaN
        # ================================
        t_in_block      = dict(
            data    = t_in_Myr,                    # (nz, ny, nx)
            x       = velr_block['x'],                  # 1D coordinate array
            y       = velr_block['y'],                  # 1D coordinate array
            z       = velr_block['z'],                  # 1D coordinate array
            extent  = velr_block['extent'],        # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk     = velr_block['ijk'],              # (block_i, block_j, block_k, block_level)
        )
        t_in_blocks.append(t_in_block)
    return {
        "df_extents"    : velr_data_dict['df_extents'],
        "blocks"        : t_in_blocks,
        "num_blocks"    : len(t_in_blocks),
        "block_shape"   : velr_data_dict['block_shape']  # (nz, ny, nx)
    }

def extract_mass_flux_data(user_params):
    rho_params                  = user_params.copy()
    rho_params.update(variable  = "dens")
    rho_data_dict               = extract_athenak_3D_block(rho_params)
    rho_blocks                  = rho_data_dict['blocks']
    velr_blocks                 = extract_velr_data(user_params)['blocks']
    mass_flux_blocks            = []
    for block_idx, (rho_block, velr_block) in enumerate(zip(rho_blocks, velr_blocks)):
        rho_data                = rho_block['data']
        rho_data_cgs            = rho_data*rho_cgs
        velr_data               = velr_block['data']
        velr_data_cgs           = velr_data*v_cgs
        mass_flux_cgs           = rho_data_cgs*velr_data_cgs
        mass_flux_unit          = mass_flux_cgs/Msolar_per_kpc_sq_per_yr
        mass_flux_block         = dict(
            data    = mass_flux_unit,                    # (nz, ny, nx)
            x       = rho_block['x'],                  # 1D coordinate array
            y       = rho_block['y'],                  # 1D coordinate array
            z       = rho_block['z'],                  # 1D coordinate array
            extent  = rho_block['extent'],        # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk     = rho_block['ijk'],              # (block_i, block_j, block_k, block_level)
        )
        mass_flux_blocks.append(mass_flux_block)
    return{
        "df_extents"    : rho_data_dict['df_extents'],
        "blocks"        : mass_flux_blocks,
        "num_blocks"    : len(mass_flux_blocks),
        "block_shape"   : rho_data_dict['block_shape']  # (nz, ny, nx)
    }
    
    
    
def compute_radial_profile(df3D, user_params, weight_type='volume'):
    """
    Helper function: Compute radial profile with flexible weighting scheme.
    
    Parameters:
    -----------
    df3D : dict
        Output from extract_athenak_3D_block containing blocks and metadata
    user_params : dict
        Parameters including center coordinates, num_radial_bins, etc.
    weight_type : str
        'volume' - constant volume weighting (all cells equally weighted)
        'mass' - mass weighting (rho * volume per cell)
    
    Returns:
    --------
    dict with keys: r_centers, profile, std, count, sum
    """
    
    blocks = df3D['blocks']
    nz, ny, nx = df3D['block_shape']
    
    x_center = user_params.get('center_x', 0.0)
    y_center = user_params.get('center_y', 0.0)
    z_center = user_params.get('center_z', 0.0)
    num_bins = user_params.get('num_radial_bins', 100)
    
    # Extract density if mass-weighting is requested
    rho_blocks = None
    if weight_type == 'mass':
        rho_params = user_params.copy()
        rho_params['variable'] = 'dens'
        rho_df3D = extract_athenak_3D_block(rho_params)
        rho_blocks = rho_df3D['blocks']
    
    radii_all = []
    data_all = []
    weights_all = []
    
    # Iterate through each block
    for block_idx, block_dict in enumerate(blocks):
        try:
            data_3d = block_dict['data']
            x_coords = block_dict['x']
            y_coords = block_dict['y']
            z_coords = block_dict['z']
            x_min, x_max, y_min, y_max, z_min, z_max = block_dict['extent']
            
            # Create 3D coordinate grids
            xx, yy, zz = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
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
            
            # Compute weights based on type
            if weight_type == 'volume':
                # Constant volume: all cells weighted equally
                weights = np.full_like(data_3d, cell_volume, dtype=np.float64)
            elif weight_type == 'mass':
                # Mass: rho * volume per cell
                rho_3d = rho_blocks[block_idx]['data']
                weights = rho_3d * cell_volume
            else:
                raise ValueError(f"Unknown weight_type: {weight_type}")
            
            # Flatten and accumulate
            radii_all.append(rr.ravel())
            data_all.append(data_3d.ravel())
            weights_all.append(weights.ravel())
            
        except Exception as e:
            print(f"  ✗ Error processing block {block_idx}: {e}")
            continue
    
    # Concatenate all blocks
    radii = np.concatenate(radii_all)
    data = np.concatenate(data_all)
    weights = np.concatenate(weights_all)
    
    # Remove NaN and invalid data
    mask = np.isfinite(data) & np.isfinite(radii) & np.isfinite(weights) & (radii >= 0) & (weights > 0)
    radii = radii[mask]
    data = data[mask]
    weights = weights[mask]
    
    # Define radial bins
    r_min = np.min(radii)
    r_max = np.max(radii)
    r_bins = np.linspace(r_min, r_max, num_bins + 1)
    r_centers = (r_bins[:-1] + r_bins[1:]) / 2
    
    # Bin and compute weighted averages
    radial_profile = np.empty(num_bins)
    radial_std = np.empty(num_bins)
    radial_count = np.zeros(num_bins, dtype=int)
    radial_sum = np.empty(num_bins)
    
    for i in range(num_bins):
        mask_shell = (radii >= r_bins[i]) & (radii < r_bins[i+1])
        n_cells = np.sum(mask_shell)
        
        if n_cells > 0:
            data_in_shell = data[mask_shell]
            weights_in_shell = weights[mask_shell]
            
            # Weighted average using numpy.average
            radial_profile[i] = np.average(data_in_shell, weights=weights_in_shell)
            radial_sum[i] = np.sum(data_in_shell * weights_in_shell)
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
    
    return {
        "r_centers": r_centers_valid,
        "profile": radial_profile_valid,
        "std": radial_std_valid,
        "count": radial_count_valid,
        "sum": radial_sum_valid
    }
