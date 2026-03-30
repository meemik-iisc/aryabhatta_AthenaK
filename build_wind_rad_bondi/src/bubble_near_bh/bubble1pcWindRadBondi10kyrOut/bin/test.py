import numpy as np
import os
# import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
fname = os.path.join(script_dir, 'wind_rad_bondi_accr.hydro_w.00000.bin')

data = np.fromfile(fname, dtype=np.float64)

# EXACT single meshblock: skip global header + first block's 6 loc doubles
skip_global = 5375 // 8 + 1  # 672
skip_loc = 6  # First block locations
skip_total = skip_global + skip_loc
print(f"Skipping {skip_total} doubles")

data = data[skip_total:]

nx1=nx2=nx3=128; nvar=6
ntotal = nx1*nx2*nx3*nvar
assert len(data) >= ntotal, f"File too small: {len(data):,} < {ntotal:,}"

u = data[:ntotal].reshape(nx3,nx2,nx1,nvar, order='F')
s00 = u[:,:,:,5]

print(f"Data loaded: {ntotal:,} doubles")
print(f"s_00 min/max: {s00.min():.8f}/{s00.max():.8f}")
print(f"frac s_00==1: {(s00 == 1.0).mean():.4f}")
print(f"s_00==0 frac: {(s00 == 0.0).mean():.4f}")