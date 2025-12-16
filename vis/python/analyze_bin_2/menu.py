# menu.py
import sys, os
from config import default_dict
from params import build_user_params, parse_int_list
from io_utils import ensure_dir

try:
    from simple_term_menu import TerminalMenu
    MENU_AVAILABLE = True
except ImportError:
    print("simple-term-menu not installed. Install with: pip install simple-term-menu")
    MENU_AVAILABLE = False

def select_analysis():
    if MENU_AVAILABLE:
        options = ["Plot 2D slice","Plot 1D profiles","Plot streamlines","Quit"]
        menu = TerminalMenu(options, title="Choose analysis:")
        idx = menu.show()
        if idx is None:
            return "quit"
        return ["slice","profiles","streamlines","quit"][idx]
    else:
        print("\n1) Plot 2D slice\n2) Plot 1D profiles\n3) Plot streamlines\n4) Quit")
        while True:
            c = input("Choice? ").strip()
            return {"1":"slice","2":"profiles","3":"streamlines","4":"quit"}.get(c)

def build_params_for(analysis_type):
    """Return a meta-dict of prompts for the chosen analysis."""
    base = {
        "input_folder": {"prompt":"Input folder","default":default_dict['input_folder'],"type":str},
        "loop_bin":     {"prompt":"Loop over all .bin files? (True/False)","default":default_dict['loop_bin'],"type":bool},
        "dt":           {"prompt":"Time step","default":default_dict['dt'],"type":float},
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
    elif analysis_type=="profiles":
        base.update({
            "variable":         {"prompt":"Variable to be plotted. check file_data['var_names'] for valid names","default":default_dict['variable'],"type":str},
            "profile_variable": {"prompt":"Profile variable [r_c or z_h or x]","default":default_dict['profile_variable'],"type":str},
            "profile_slice":    {"prompt":"Select axis of 1D slice","default":default_dict['profile_slice'],"type":int},
            "axes_scale":       {"prompt":"Scale Axis","default":default_dict['axes_scale'],"type":str},
            "cmap_label":   {"prompt":"Color label","default":default_dict['cmap_label'],"type":str},
            "clim":         {"prompt":"Colorbar limits: [vmin,vmax]","default":default_dict['clim'],"type":parse_int_list},
            "color":   {"prompt":"Plot color","default":default_dict['color'],"type":str},
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
