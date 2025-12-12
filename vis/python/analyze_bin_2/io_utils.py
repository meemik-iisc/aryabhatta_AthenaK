from pathlib import Path
import re

def ensure_dir(path: str):
    Path(path).mkdir(parents=True, exist_ok=True)

def extract_slice_number(filename: str) -> int:
    match = re.search(r'\.(\d+)\.bin$', filename.strip())
    if match:
        return int(match.group(1))
    raise ValueError(f"Could not extract slice number from filename: {filename}")
