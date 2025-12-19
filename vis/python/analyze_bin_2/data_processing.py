import struct
import numpy as np
import pandas as pd
from constants import mu,mp_cgs,kB_cgs,L_code,M_code,t_code,gamma,Temp_norm

def extract_athenak_slice(user_params):
    """
    Extract 2D slice data of requested variable from AthenaK binary file and return as pandas DataFrames.
    
    user_params: dict with keys including
        - input_folder: folder path of input file
        - input_file: filename of AthenaK binary dump
        - variable: variable name to extract
        - direction: slice direction, e.g. 'x1', 'x2', 'x3' or 'x', 'y', 'z'
        - location: float coordinate of slice position in slice direction
    
    Returns:
        df_quantities: pandas DataFrame with each block's 2D slice flattened and indexed by block
        df_extents: pandas DataFrame of block extents with columns ['x_min', 'x_max', 'y_min', 'y_max']
        additional metadata (num_blocks, block shape)
    """
    file_path = user_params["input_folder"].rstrip('/') + '/' + user_params["input_file"]
    variable = user_params["variable"]
    if variable.startswith("derived:"):
        variable=variable.split(":", 1)[1].strip()
        if variable == "temp":
            return extract_temp_slice(user_params)
        elif variable == "velr":
            return extract_radial_vel_slice(user_params)
    direction = user_params.get("direction", None)
    if direction in ('x1', '1'):
        dimension = 'x'
    elif direction in ('x2', '2'):
        dimension = 'y'
    elif direction in ('x3', '3'):
        dimension = 'z'
    else:
        dimension = direction  # fallback
    location = user_params.get("location", 0.0)
    
    with open(file_path, "rb") as f:
        line = f.readline().decode('ascii')
        if line != 'Athena binary output version=1.1\n':
            raise RuntimeError('Unrecognized AthenaK binary file format.')
        for _ in range(3):
            next(f)
        location_size = int(f.readline().decode('ascii')[19:])
        variable_size = int(f.readline().decode('ascii')[19:])
        next(f)
        variable_names_base = f.readline().decode('ascii')[12:].split()
        header_offset = int(f.readline().decode('ascii')[16:])
        
        if variable[:8] == 'derived:':
            raise NotImplementedError("Derived variable extraction not implemented here.")
        else:
            var_name = variable
        
        if var_name not in variable_names_base:
            raise RuntimeError(f'Variable "{var_name}" not found in file variables.')
        var_ind = variable_names_base.index(var_name)
        
        location_fmt = 'f' if location_size == 4 else 'd'
        variable_fmt = 'f' if variable_size == 4 else 'd'
        
        input_data = {}
        start_data = f.tell() + header_offset
        while f.tell() < start_data:
            line = f.readline().decode('ascii')
            if line.startswith('#'):
                continue
            if line.startswith('<'):
                section = line[1:-2]
                input_data[section] = {}
                continue
            if '=' in line:
                k, v = line.split('=', 1)
                input_data[section][k.strip()] = v.split('#')[0].strip()
        
        num_ghost = int(input_data['mesh']['nghost'])
        num_variables_base = len(variable_names_base)
        quantities_list = []
        extents_list = []
        file_size = f.seek(0, 2)
        f.seek(start_data, 0)
        
        max_level_calculated = -1
        block_loc_for_level = []
        block_ind_for_level = []
        num_blocks_used = 0
        
        while f.tell() < file_size:
            block_indices = np.array(struct.unpack('@6i', f.read(24))) - num_ghost
            block_i, block_j, block_k, block_level = struct.unpack('@4i', f.read(16))
            
            block_nx = block_indices[1] - block_indices[0] + 1
            block_ny = block_indices[3] - block_indices[2] + 1
            block_nz = block_indices[5] - block_indices[4] + 1
            cells_per_block = block_nz * block_ny * block_nx
            block_cell_fmt = '=' + str(cells_per_block) + variable_fmt
            variable_data_size = cells_per_block * variable_size
            
            if dimension is None:
                if block_nx > 1 and block_ny > 1 and block_nz > 1:
                    dim = 'z'
                elif block_nx > 1 and block_ny > 1:
                    dim = 'z'
                elif block_nx > 1 and block_nz > 1:
                    dim = 'y'
                elif block_ny > 1 and block_nz > 1:
                    dim = 'x'
                else:
                    raise RuntimeError('Input file only contains 1D data.')
            else:
                dim = dimension
            
            if dim == 'x':
                if block_ny == 1 or block_nz == 1:
                    raise RuntimeError('Data in file has no extent in required dims')
                block_nx1 = block_ny
                block_nx2 = block_nz
                slice_block_n = block_nx
                slice_location_min = float(input_data['mesh']['x1min'])
                slice_location_max = float(input_data['mesh']['x1max'])
                slice_root_blocks = int(input_data['mesh']['nx1']) // int(input_data['meshblock']['nx1'])
            elif dim == 'y':
                if block_nx == 1 or block_nz == 1:
                    raise RuntimeError('Data in file has no extent in required dims')
                block_nx1 = block_nx
                block_nx2 = block_nz
                slice_block_n = block_ny
                slice_location_min = float(input_data['mesh']['x2min'])
                slice_location_max = float(input_data['mesh']['x2max'])
                slice_root_blocks = int(input_data['mesh']['nx2']) // int(input_data['meshblock']['nx2'])
            else:  # z
                if block_nx == 1 or block_ny == 1:
                    raise RuntimeError('Data in file has no extent in required dims')
                block_nx1 = block_nx
                block_nx2 = block_ny
                slice_block_n = block_nz
                slice_location_min = float(input_data['mesh']['x3min'])
                slice_location_max = float(input_data['mesh']['x3max'])
                slice_root_blocks = int(input_data['mesh']['nx3']) // int(input_data['meshblock']['nx3'])
            
            slice_normalized_coord = (location - slice_location_min) / (slice_location_max - slice_location_min)
            
            if block_level > max_level_calculated:
                for level in range(max_level_calculated + 1, block_level + 1):
                    if location <= slice_location_min:
                        block_loc_for_level.append(0)
                        block_ind_for_level.append(0)
                    elif location >= slice_location_max:
                        block_loc_for_level.append(slice_root_blocks - 1)
                        block_ind_for_level.append(slice_block_n - 1)
                    else:
                        slice_mesh_n = slice_block_n * slice_root_blocks * 2**level
                        mesh_ind = int(slice_normalized_coord * slice_mesh_n)
                        block_loc_for_level.append(mesh_ind // slice_block_n)
                        block_ind_for_level.append(mesh_ind - slice_block_n * block_loc_for_level[-1])
                max_level_calculated = block_level
            
            if (dim == 'x' and block_i != block_loc_for_level[block_level]) or \
               (dim == 'y' and block_j != block_loc_for_level[block_level]) or \
               (dim == 'z' and block_k != block_loc_for_level[block_level]):
                f.seek(6*location_size + num_variables_base*variable_data_size, 1)
                continue
            
            num_blocks_used += 1
            
            block_lims = struct.unpack('=6'+location_fmt, f.read(6*location_size))
            if dim == 'x':
                extents_list.append((block_lims[2], block_lims[3], block_lims[4], block_lims[5]))
            elif dim == 'y':
                extents_list.append((block_lims[0], block_lims[1], block_lims[4], block_lims[5]))
            else:  # z
                extents_list.append((block_lims[0], block_lims[1], block_lims[2], block_lims[3]))
            
            cell_data_start = f.tell()
            f.seek(cell_data_start + var_ind * variable_data_size, 0)
            cell_data = np.array(struct.unpack(block_cell_fmt, f.read(variable_data_size)))
            cell_data = cell_data.reshape((block_nz, block_ny, block_nx))
            if dim == 'x':
                slice_2d = cell_data[:, :, block_ind_for_level[block_level]]
            elif dim == 'y':
                slice_2d = cell_data[:, block_ind_for_level[block_level], :]
            else:  # z
                slice_2d = cell_data[block_ind_for_level[block_level], :, :]
            
            # Convert the 2D block slice to a DataFrame indexed by coordinates in grid order
            rows, cols = slice_2d.shape
            df_block = pd.DataFrame(slice_2d, index=range(rows), columns=range(cols))
            quantities_list.append(df_block)
            f.seek(cell_data_start + num_variables_base * variable_data_size, 0)
        
        # Create DataFrame for extents
        df_extents = pd.DataFrame(extents_list, columns=['x_min', 'x_max', 'y_min', 'y_max'])
        
        quantities_stacked = [df.stack() for df in quantities_list]
        combined_df = pd.concat(quantities_stacked, keys=range(len(quantities_stacked)), names=['block', 'row', 'col'])
        
        return {
            "df_quantities": combined_df,
            "df_extents": df_extents,
            "num_blocks": num_blocks_used,
            "block_shape": (block_nx2, block_nx1)
        }
        
        
def extract_temp_slice(user_params):
    rho_params = user_params.copy()
    rho_params.update(variable="dens")
    rho_data_df=extract_athenak_slice(rho_params)
    rho_data=rho_data_df['df_quantities']
    pres_params = user_params.copy()
    pres_params.update(variable="eint")
    pres_data_df = extract_athenak_slice(pres_params)
    pres_data = pres_data_df['df_quantities']
    temp_data=(pres_data/rho_data)* (mu*mp_cgs/(kB_cgs*(gamma-1)))*((L_code/t_code)**2)/Temp_norm
    return {
        "df_quantities": temp_data,
        "df_extents": pres_data_df['df_extents'],
        "num_blocks": pres_data_df['num_blocks'],
        "block_shape": pres_data_df['block_shape']
    }

def extract_radial_vel_slice(user_params):
    vel_params = user_params.copy()
    vel_params.update(variable="velx")
    velx_data_df=extract_athenak_slice(vel_params)
    velx_data = velx_data_df['df_quantities']
    vel_params.update(variable="vely")
    vely_data = extract_athenak_slice(vel_params)['df_quantities']
    vel_params.update(variable="velz")
    velz_data = extract_athenak_slice(vel_params)['df_quantities']
    velr_data = np.sqrt(np.square(velx_data)+np.square(vely_data)+np.square(velz_data))
    return {
        "df_quantities": velr_data,
        "df_extents": velx_data_df['df_extents'],
        "num_blocks": velx_data_df['num_blocks'],
        "block_shape": velx_data_df['block_shape']
    }
import numpy as np

def stitch_meshblocks_to_global(data_dict, user_params):
    """
    Stitch multiple meshblocks (2D slice blocks) into one global 2D numpy array.
    
    Parameters:
    - data_dict: dict output from extract_athenak_slice_as_dataframe, containing:
        'df_quantities', 'df_extents', 'block_shape'
    - user_params: dict containing 'axes_scale' (float) and optionally 'fill_value' to fill gaps
    
    Returns:
    - global_array: 2D numpy array combining all blocks spatially
    - global_extent: tuple (x_min, x_max, y_min, y_max) of stitched array in scaled coordinates
    """
    df_quantities = data_dict['df_quantities']
    df_extents = data_dict['df_extents']
    block_nx2, block_nx1 = data_dict['block_shape']
    axes_scale = user_params.get('axes_scale', 1.0)
    fill_value = user_params.get('fill_value', np.nan)  # value to fill empty regions if any
    
    # Scale extents
    scaled_extents = [tuple(axes_scale * x for x in extent) for extent in df_extents.values]
    extents_array = np.array(scaled_extents)
    
    # Determine global extent boundaries
    x_min = extents_array[:, 0].min()
    x_max = extents_array[:, 1].max()
    y_min = extents_array[:, 2].min()
    y_max = extents_array[:, 3].max()
    
    # Determine resolution from first block (assumes uniform resolution)
    res_x = (extents_array[0, 1] - extents_array[0, 0]) / block_nx1
    res_y = (extents_array[0, 3] - extents_array[0, 2]) / block_nx2
    
    # Compute global grid shape
    nx_total = int(round((x_max - x_min) / res_x))
    ny_total = int(round((y_max - y_min) / res_y))
    
    # Initialize global array
    global_array = np.full((ny_total, nx_total), fill_value, dtype=float)
    
    for block_idx, extent in enumerate(extents_array):
        block_data = df_quantities.loc[block_idx].unstack(level=-1).values
        
        # Compute pixel indices for this block within global array
        ix_start = int(round((extent[0] - x_min) / res_x))
        ix_end = ix_start + block_nx1
        iy_start = int(round((extent[2] - y_min) / res_y))
        iy_end = iy_start + block_nx2
        
        # Place block data in global array
        global_array[iy_start:iy_end, ix_start:ix_end] = block_data
    
    global_extent = (x_min, x_max, y_min, y_max)
    return global_array, global_extent

