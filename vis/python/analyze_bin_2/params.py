import ast
from config import default_dict
from pathlib import Path
from utils.io_utils import extract_slice_number, sort_by_slice_number
def parse_int_list(s: str) -> list[int]:
    s = s.strip()
    try:
        val = ast.literal_eval(s)
        if isinstance(val, list) and all(isinstance(x, int) for x in val):
            return val
    except Exception:
        pass
    return [int(x) for x in s.split(',') if x.strip()]

def prompt_params(meta: dict) -> dict:
    filled = {}
    for key, info in meta.items():
        prompt = f"{info['prompt']} [{info['default']}]: "
        val = input(prompt).strip()
        if not val:
            filled[key] = info["default"]
        else:
            try:
                filled[key] = info["type"](val)
            except Exception:
                filled[key] = val  # fallback to raw
        print(f" → {key} = {filled[key]!r}")
    return filled

def build_user_params(meta: dict) -> dict:
    p = prompt_params(meta)

    inp = Path(p['input_folder'])
    if not p['loop_bin']:
        fn = input(f"Input file name [{default_dict['input_file']}]: ").strip() or default_dict['input_file']
        p['input_file'] = fn
        print(f" → input_file = {p['input_file']!r}")

    if p['loop_bin']:
        input_files = [f.name for f in inp.glob("*.bin")]
        input_files_sorted = sort_by_slice_number(input_files)
        p['input_files'] = input_files_sorted
    else:
        p['bin_path'] = inp / p['input_file']
        p['slice_number'] = extract_slice_number(p['input_file'])

    import os
    p['output_path'] = os.path.join(os.path.dirname(p['input_folder']), "output_folder")
    return p

def rad_build_params(user_params_dict: dict)->dict:
    rad_params_dict = user_params_dict.copy()
    rad_input_folder = user_params_dict['rad_input_folder']
    rad_input_path = Path(rad_input_folder)
    rad_params_dict.update(input_folder = rad_input_folder)
    if not user_params_dict['loop_bin']:
        fn = input(f"Input Radiative file name [{default_dict['rad_input_file']}]: ").strip() or default_dict['rad_input_file']
        rad_params_dict['input_file'] = fn
        print(f" → rad_input_file = {user_params_dict['input_file']!r}")

    if user_params_dict['loop_bin']:
        rad_input_files = [f.name for f in rad_input_path.glob("*.bin")]
        rad_input_files_sorted = sort_by_slice_number(rad_input_files)
        rad_params_dict['input_files'] = rad_input_files_sorted
        slice_number = user_params_dict['slice_number']
        rad_params_dict['input_file'] = rad_input_files_sorted[slice_number]
        rad_params_dict['bin_path'] = rad_input_path/rad_params_dict['input_file']
    else:
        rad_params_dict['bin_path'] = rad_input_path / rad_params_dict['input_file']

    import os
    rad_params_dict['output_path'] = os.path.join(os.path.dirname(rad_params_dict['input_folder']), "output_folder")
    return rad_params_dict