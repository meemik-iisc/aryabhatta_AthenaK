import struct
import numpy as np
import pandas as pd
from constants import mu,mp_cgs,kB_cgs,s_Myr,length_cgs,mass_cgs,time_cgs,pres_cgs,rho_cgs,gamma,Temp_norm,Velr_scale
from ismcooling import cool_lambda


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
    file_path = user_params["input_folder"].rstrip('/') + '/' + user_params["input_file"]
    variable = user_params["variable"]
    if variable.startswith('derived:'):
        variable=variable.split(":", 1)[1].strip()
    else:
        var_name = variable
        
        
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
        file_size = f.seek(0, 2)
        f.seek(start_data, 0)

        quantities_list = []
        extents_list = []
        blocks = []
        num_blocks_used = 0
        block_nx_ref = block_ny_ref = block_nz_ref = None

        while f.tell() < file_size:
            # indices and AMR level
            block_indices = np.array(struct.unpack('@6i', f.read(24))) - num_ghost
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

            # store per-block numpy representation for fast later use
            x_coords = np.linspace(x_min, x_max, block_nx)
            y_coords = np.linspace(y_min, y_max, block_ny)
            z_coords = np.linspace(z_min, z_max, block_nz)
            blocks.append(
                dict(
                    data=cell_data,
                    x=x_coords,
                    y=y_coords,
                    z=z_coords,
                    extent=(x_min, x_max, y_min, y_max, z_min, z_max),
                    ijk=(block_i, block_j, block_k, block_level),
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
            "df_extents": df_extents,
            "blocks": blocks,
            "num_blocks": num_blocks_used,
            "block_shape": (block_nz_ref, block_ny_ref, block_nx_ref)  # (nz, ny, nx)
        }
        
def extract_temp_data(user_params):
    rho_params = user_params.copy()
    rho_params['variable'] = 'dens'
    rho_data_dict = extract_athenak_3D_block(rho_params)
    rho_blocks = rho_data_dict['blocks']
    
    # Extract internal energy (eint) as 3D blocks
    eint_params = user_params.copy()
    eint_params['variable'] = 'eint'
    eint_data_dict = extract_athenak_3D_block(eint_params)
    eint_blocks = eint_data_dict['blocks']
    
    # Compute temperature for each block: T = (eint/rho) * mu*mp/(kB*(gamma-1)) * (L/t)^2 / Temp_norm
    temp_blocks = []
    
    for block_idx, (rho_block, eint_block) in enumerate(zip(rho_blocks, eint_blocks)):
        rho_data = rho_block['data']       # (nz, ny, nx)
        eint_data = eint_block['data']     # (nz, ny, nx)
        
        # Avoid division by zero
        rho_safe = np.where(rho_data > 0, rho_data, np.nan)
        
        # Temperature calculation
        temp_data = (eint_data / rho_safe) * \
                    (mu * mp_cgs / (kB_cgs * (gamma - 1))) * \
                    ((length_cgs / time_cgs) ** 2) / Temp_norm
        
        # Create temperature block dict (same structure as input blocks)
        temp_block = dict(
            data=temp_data,                    # (nz, ny, nx)
            x=rho_block['x'],                  # 1D coordinate array
            y=rho_block['y'],                  # 1D coordinate array
            z=rho_block['z'],                  # 1D coordinate array
            extent=rho_block['extent'],        # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk=rho_block['ijk'],              # (block_i, block_j, block_k, block_level)
        )
        temp_blocks.append(temp_block)
    
    return {
        "df_extents": eint_data_dict['df_extents'],
        "blocks": temp_blocks,
        "num_blocks": len(temp_blocks),
        "block_shape": eint_data_dict['block_shape']  # (nz, ny, nx)
    }

def extract_velr_data(user_params):
    vel_params = user_params.copy()
    vel_params.update(variable="velx")
    velx_data_dict=extract_athenak_3D_block(vel_params)
    velx_blocks = velx_data_dict['blocks']
    vel_params.update(variable="vely")
    vely_blocks = extract_athenak_3D_block(vel_params)['blocks']
    vel_params.update(variable="velz")
    velz_blocks = extract_athenak_3D_block(vel_params)['blocks']
    
    velr_blocks = []
    
    for block_idx, (velx_block, vely_block, velz_block) in enumerate(zip(velx_blocks, vely_blocks, velz_blocks)):
        velx_data = velx_block['data']       # (nz, ny, nx)
        vely_data = vely_block['data']     # (nz, ny, nx)
        velz_data = velz_block['data']
                
        # Temperature calculation
        velr_data = np.sqrt(np.square(velx_data)+np.square(vely_data)+np.square(velz_data))*Velr_scale
        
        # Create velocity block dict (same structure as input blocks)
        velr_block = dict(
            data=velr_data,                    # (nz, ny, nx)
            x=velx_block['x'],                  # 1D coordinate array
            y=velx_block['y'],                  # 1D coordinate array
            z=velx_block['z'],                  # 1D coordinate array
            extent=velx_block['extent'],        # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk=velx_block['ijk'],              # (block_i, block_j, block_k, block_level)
        )
        velr_blocks.append(velr_block)
    
    return {
        "df_extents": velx_data_dict['df_extents'],
        "blocks": velr_blocks,
        "num_blocks": len(velr_blocks),
        "block_shape": velx_data_dict['block_shape']  # (nz, ny, nx)
    }
def extract_cooling_rate_data(user_params):
    eint_params = user_params.copy()
    eint_params.update(variable="eint")
    eint_data_dict=extract_athenak_3D_block(eint_params)
    eint_blocks = eint_data_dict['blocks']
    rho_params = user_params.copy()
    rho_params.update(variable="dens")
    rho_blocks = extract_athenak_3D_block(rho_params)['blocks']
    
    cooling_rate_blocks = []
    
    for block_idx, (eint_block, rho_block) in enumerate(zip(eint_blocks, rho_blocks)):
        eint_data = eint_block['data']
        rho_data = rho_block['data']
        
        #Convert to CGS
        pres_data_cgs = eint_data*(gamma-1)*pres_cgs
        rho_data_cgs = rho_data*rho_cgs
        
        # Temperature calculation
        temp_data_cgs = (pres_data_cgs*mu*mp_cgs)/(rho_data_cgs*kB_cgs)
        #Cooling Rate calculation
        lambda_cgs = cool_lambda(temp_data_cgs)
        #Calculate cooling time
        cooling_rate_cgs = ((rho_data_cgs**2)*lambda_cgs)/((mu**2)*(mp_cgs**2))
        # Create velocity block dict (same structure as input blocks)
        cooling_rate_block = dict(
            data=cooling_rate_cgs,                    # (nz, ny, nx)
            x=rho_block['x'],                  # 1D coordinate array
            y=rho_block['y'],                  # 1D coordinate array
            z=rho_block['z'],                  # 1D coordinate array
            extent=rho_block['extent'],        # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk=rho_block['ijk'],              # (block_i, block_j, block_k, block_level)
        )
        cooling_rate_blocks.append(cooling_rate_block)
    
    return {
        "df_extents": eint_data_dict['df_extents'],
        "blocks": cooling_rate_blocks,
        "num_blocks": len(cooling_rate_blocks),
        "block_shape": eint_data_dict['block_shape']  # (nz, ny, nx)
    }
    
def extract_cool_time_data(user_params):
    eint_params = user_params.copy()
    eint_params.update(variable="eint")
    eint_data_dict=extract_athenak_3D_block(eint_params)
    eint_blocks = eint_data_dict['blocks']
    rho_params = user_params.copy()
    rho_params.update(variable="dens")
    rho_blocks = extract_athenak_3D_block(rho_params)['blocks']
    
    tcool_blocks = []
    
    for block_idx, (eint_block, rho_block) in enumerate(zip(eint_blocks, rho_blocks)):
        eint_data = eint_block['data']
        rho_data = rho_block['data']
        
        #Convert to CGS
        pres_data_cgs = eint_data*(gamma-1)*pres_cgs
        rho_data_cgs = rho_data*rho_cgs
        
        # Temperature calculation
        temp_data_cgs = (pres_data_cgs*mu*mp_cgs)/(rho_data_cgs*kB_cgs)
        #Cooling Rate calculation
        lambda_cgs = cool_lambda(temp_data_cgs)
        #Calculate cooling time
        t_cool_cgs = (gamma*pres_data_cgs*(mu**2)*(mp_cgs**2))/((gamma-1)*(rho_data_cgs**2)*lambda_cgs)
        t_cool_Myr = t_cool_cgs/s_Myr
        # Create velocity block dict (same structure as input blocks)
        tcool_block = dict(
            data=t_cool_Myr,                    # (nz, ny, nx)
            x=rho_block['x'],                  # 1D coordinate array
            y=rho_block['y'],                  # 1D coordinate array
            z=rho_block['z'],                  # 1D coordinate array
            extent=rho_block['extent'],        # (x_min, x_max, y_min, y_max, z_min, z_max)
            ijk=rho_block['ijk'],              # (block_i, block_j, block_k, block_level)
        )
        tcool_blocks.append(tcool_block)
    
    return {
        "df_extents": eint_data_dict['df_extents'],
        "blocks": tcool_blocks,
        "num_blocks": len(tcool_blocks),
        "block_shape": eint_data_dict['block_shape']  # (nz, ny, nx)
    }
