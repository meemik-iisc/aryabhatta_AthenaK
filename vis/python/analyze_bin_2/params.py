import ast
from config import default_dict

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
    from pathlib import Path
    inp = Path(p['input_folder'])
    if not p['loop_bin']:
        fn = input(f"Input file name [{default_dict['input_file']}]: ").strip() or default_dict['input_file']
        p['input_file'] = fn
        print(f" → input_file = {p['input_file']!r}")

    if p['loop_bin']:
        p['input_files'] = [f.name for f in inp.glob("*.bin")]
    else:
        p['bin_path'] = inp / p['input_file']
        from io_utils import extract_slice_number
        p['slice_number'] = extract_slice_number(p['input_file'])

    import os
    p['output_path'] = os.path.join(os.path.dirname(p['input_folder']), "output_folder")
    return p
