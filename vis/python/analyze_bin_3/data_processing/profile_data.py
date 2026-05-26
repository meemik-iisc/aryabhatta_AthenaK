import struct
import numpy as np
import pandas as pd
from utils.constants import mu,mp_cgs,kB_cgs,length_cgs,mass_cgs,time_cgs,gamma,Temp_norm,rho_cgs,pres_cgs,s_Myr, Velr_scale
from utils.helpers import cool_lambda, mask_half

def get_zh_data(data_dict, user_params):
    """
    Returns concatenated 1D vertical cut (ZH) profile from AthenaK blocks along y=user_params['profile_slice'].
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

    return x_filtered, data_filtered