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
from data_processing import extract_athenak_slice, stitch_meshblocks_to_global
from plotting.slice_plot import plot_athenak_combined, plot_individual_blocks, plot_stitched_data 
from plotting.plot_1d_profiles import plot_zh,plot_rc
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
            df = extract_athenak_slice(user_params)
            if user_params['profile_variable']=='zh':
                plot_zh(df, user_params)
            elif user_params['profile_variable']=='rc':
                plot_rc(df,user_params)
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
