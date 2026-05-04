# menu.py
import sys, os
from config import default_dict
from params import build_user_params, parse_int_list, parse_str_list, parse_str_pair_list, parse_float_pair_list
from utils.io_utils import ensure_dir

try:
    from simple_term_menu import TerminalMenu
    MENU_AVAILABLE = True
except ImportError:
    print("simple-term-menu not installed. Install with: pip install simple-term-menu")
    MENU_AVAILABLE = False

def select_analysis():
    if MENU_AVAILABLE:
        options = ["Plot 2D slice", "Plot split slice", "Plot 1D profiles","Plot Comparative Profile","Plot streamlines","Quit"]
        menu = TerminalMenu(options, title="Choose analysis:")
        idx = menu.show()
        if idx is None:
            return "quit"
        return ["slice","split_slice","profiles","comparative","streamlines","quit"][idx]
    else:
        print("\n1) Plot 2D slice \n2) Plot split slice\n3) Plot 1D profiles\n4) Plot Comparative Profile\n5) Plot streamlines\n6) Quit")
        while True:
            c = input("Choice? ").strip()
            return {"1":"slice","2":"split_slice","3":"profiles","4":"comparative","5":"streamlines","6":"quit"}.get(c)

def build_params_for(analysis_type):
    """Return a meta-dict of prompts for the chosen analysis."""
    base = {
        "input_folder": {"prompt":"Input folder","default":default_dict['input_folder'],"type":str},
        "loop_bin":     {"prompt":"Loop over all .bin files? (True/False)","default":default_dict['loop_bin'],"type":bool},
        "dt":           {"prompt":"Time step","default":default_dict['dt'],"type":float},
        "time_label":   {"prompt":"Time label","default":default_dict['time_label'],"type":str},
        "direction":    {"prompt":"Slice direction: ['x1','x2','x3', or None], if None default to 'x2'","default":default_dict['direction'],"type":str},
        "location":     {"prompt":"Slice location","default":default_dict['location'],"type":int},
    }
    if analysis_type=="slice":
        base.update({
            "variable":     {"prompt":"Variable to be plotted. check file_data['var_names'] for valid names","default":default_dict['variable'],"type":str},
            "cmap":         {"prompt":"Color map","default":default_dict['cmap'],"type":str},
            "cmap_label":   {"prompt":"Color label","default":default_dict['cmap_label'],"type":str},
            "norm":         {"prompt":"Color map normalization:[None=linear,'log'=logarithmic]","default":default_dict['norm'],"type":str},
            "clim":         {"prompt":"Colorbar limits: [vmin,vmax]","default":default_dict['clim'],"type":parse_int_list},
            "extents":      {"prompt":"Plot extents: [x1min,x1max,x2min,x2max]","default":default_dict['extents'],"type":parse_int_list},
            "figsize":      {"prompt":"Figure size","default":default_dict['figsize'],"type":tuple},
            "xlabel":       {"prompt":"X label","default":default_dict['xlabel'],"type":str},
            "ylabel":       {"prompt":"Y label","default":default_dict['ylabel'],"type":str},
            "axes_scale":   {"prompt":"Scale Axis","default":default_dict['axes_scale'],"type":str},
        })
    elif analysis_type=="split_slice":
        base.update({
            "variable":     {"prompt":"List of variables [top, bottom] to be plotted.","default":default_dict['variable'],"type":parse_str_pair_list},
            "cmap":         {"prompt":"List of color map [top, bottom]","default":default_dict['cmap'],"type":parse_str_pair_list},
            "cmap_label":   {"prompt":"List of color label","default":default_dict['cmap_label'],"type":parse_str_pair_list},
            "norm":         {"prompt":"List of color map normalization:[None=linear,'log'=logarithmic]","default":default_dict['norm'],"type":parse_str_pair_list},
            "clim":         {"prompt":"List of colorbar limits: [vmin_top,vmax_top,vmin_bot,vmax_bot]","default":default_dict['clim'],"type":parse_float_pair_list},
            "extents":      {"prompt":"Plot extents: [x1min,x1max,x2min,x2max]","default":default_dict['extents'],"type":parse_float_pair_list},
            "figsize":      {"prompt":"Figure size","default":default_dict['figsize'],"type":tuple},
            "xlabel":       {"prompt":"X label","default":default_dict['xlabel'],"type":str},
            "ylabel":       {"prompt":"Y label","default":default_dict['ylabel'],"type":str},
            "axes_scale":   {"prompt":"Scale Axis","default":default_dict['axes_scale'],"type":str},
        })
    elif analysis_type=="profiles":
        base.update({
            "variable":         {"prompt":"Variable to be plotted. check file_data['var_names'] for valid names","default":default_dict['variable'],"type":str},
            "profile_variable": {"prompt":"Profile variable [r_c or z_h or x or avg]","default":default_dict['profile_variable'],"type":str},
            "profile_slice":    {"prompt":"Select axis of 1D slice","default":default_dict['profile_slice'],"type":int},
            "axes_scale":       {"prompt":"Scale Axis","default":default_dict['axes_scale'],"type":str},
            "cmap_label":       {"prompt":"Color label","default":default_dict['cmap_label'],"type":str},
            "norm":             {"prompt":"Y axis normalization:[None=linear,'log'=logarithmic]","default":default_dict['norm'],"type":str},
            "clim":             {"prompt":"Colorbar limits: [vmin,vmax]","default":default_dict['clim'],"type":parse_int_list},
            "color":            {"prompt":"Plot color","default":default_dict['color'],"type":str},
            "xlabel":           {"prompt":"X label","default":default_dict['xlabel'],"type":str},
            "axes_scale":       {"prompt":"Scale Axis","default":default_dict['axes_scale'],"type":str},
            })
    elif analysis_type=="comparative":
        base.update({
            "rad_input_folder": {"prompt":"Radiative Cooling Input folder","default":default_dict['rad_input_folder'],"type":str},
            "variable":         {"prompt":"Variable to be plotted. check file_data['var_names'] for valid names","default":default_dict['variable'],"type":str},
            "profile_variable": {"prompt":"Profile variable [r_c or z_h or x or avg]","default":default_dict['profile_variable'],"type":str},
            "profile_slice":    {"prompt":"Select axis of 1D slice","default":default_dict['profile_slice'],"type":int},
            "axes_scale":       {"prompt":"Scale Axis","default":default_dict['axes_scale'],"type":str},
            "cmap_label":       {"prompt":"Color label","default":default_dict['cmap_label'],"type":str},
            "norm":             {"prompt":"Y axis normalization:[None=linear,'log'=logarithmic]","default":default_dict['norm'],"type":str},
            "clim":             {"prompt":"Colorbar limits: [vmin,vmax]","default":default_dict['clim'],"type":parse_int_list},
            "color":            {"prompt":"Plot color","default":default_dict['color'],"type":str},
            "xlabel":           {"prompt":"X label","default":default_dict['xlabel'],"type":str},
            "axes_scale":       {"prompt":"Scale Axis","default":default_dict['axes_scale'],"type":str},
            })
    else:  # streamlines
        base.update({
            "cmap":                 {"prompt":"Color map","default":default_dict['cmap'],"type":str},
            "norm":                 {"prompt":"Color map normalization:[None=linear,'log'=logarithmic]","default":default_dict['norm'],"type":str},
            "streamline_density":   {"prompt":"Streamline density","default":default_dict['streamline_density'],"type":float},
            "extents":              {"prompt":"Plot extents: [x1min,x1max,x2min,x2max]","default":default_dict['extents'],"type":parse_int_list},
            "figsize":              {"prompt":"Figure size","default":default_dict['figsize'],"type":tuple}
        })
    return build_user_params(base)
