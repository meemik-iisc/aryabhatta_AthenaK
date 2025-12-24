# analysis.py
import os
import sys
import numpy as np

# Ensure parent directory (where bin_convert.py lives) is on sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)
    
import bin_convert
from io_utils import ensure_dir, extract_slice_number
from data_processing.slice_data import extract_athenak_slice, stitch_meshblocks_to_global
from data_processing.block_data import extract_athenak_3D_block, extract_temp_data, extract_velr_data, extract_cooling_rate_data, extract_cool_time_data
from plotting.slice_plot import plot_athenak_combined, plot_individual_blocks, plot_stitched_data 
from plotting.plot_1d_profiles import plot_zh,plot_rc,plot_x_profile, plot_spherical_volume_weighted_avg_profile, plot_spherical_mass_weighted_avg_profile
from plotting.streamlines import plot_streamlines_from_dataframes

def run(analysis_type, user_params):
    ensure_dir(user_params['output_path'])
    def handle():
        if analysis_type=="slice":
            file_data_2d = extract_athenak_slice(user_params)
            stitch_data, stitch_extent = stitch_meshblocks_to_global(file_data_2d,user_params)
            plot_stitched_data(stitch_data, stitch_extent, user_params)
            # print(file_data_2d['df_quantities'].head)
            # plot_athenak_combined_2(file_data_2d,user_params)
            # plot_individual_meshblocks(df_2d,user_params)       
            # plot_slice(df_2d, user_params)1

        elif analysis_type=="profiles":
            if user_params['profile_variable']=='avg':
                
                variable = user_params["variable"]
                if variable.startswith('derived:'):
                    variable=variable.split(":", 1)[1].strip()
                    if variable == "temp":
                        df3D = extract_temp_data(user_params)
                        plot_spherical_volume_weighted_avg_profile(df3D, user_params)
                    elif variable == "velr":
                        df3D = extract_velr_data(user_params)
                        plot_spherical_mass_weighted_avg_profile(df3D, user_params)
                    elif variable == "cooling_rate":
                        df3D = extract_cooling_rate_data(user_params)
                        plot_spherical_volume_weighted_avg_profile(df3D,user_params)
                    elif variable == "tcool":
                        df3D = extract_cool_time_data(user_params)
                        plot_spherical_mass_weighted_avg_profile(df3D,user_params)
                else:
                    df3D = extract_athenak_3D_block(user_params)
                    plot_spherical_volume_weighted_avg_profile(df3D, user_params)
                
            else:
                df = extract_athenak_slice(user_params)
                if user_params['profile_variable']=='zh':
                    plot_zh(df, user_params)
                elif user_params['profile_variable']=='rc':
                    plot_rc(df,user_params)
                elif user_params['profile_variable']=='x':
                    plot_x_profile(df,user_params)
            
        else:
            p_u = user_params.copy(); p_u['variable']="velx"
            p_v = user_params.copy(); p_v['variable']="velz"
            df_u = extract_athenak_slice(p_u)
            df_v = extract_athenak_slice(p_v)
            plot_streamlines_from_dataframes(df_u, df_v, user_params)

    if user_params['loop_bin']:
        for fn in user_params['input_files']:
            fp = os.path.join(user_params['input_folder'], fn)
            data = bin_convert.read_binary(fp)
            user_params.update(input_file=fn,
                               bin_path=fp,
                               slice_number=extract_slice_number(fn))
            handle()
    else:
        # data = bin_convert.read_binary(user_params['bin_path'])
        # print(list(np.array(data['mb_logical'])[:,-1]))
        user_params.update(slice_number=extract_slice_number(user_params['input_file']))
        handle()
