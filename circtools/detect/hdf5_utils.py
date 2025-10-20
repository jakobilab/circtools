# utils_hdf5.py
import h5py, numpy as np, pandas as pd, os

def write_dataframe_to_hdf5(df, h5_path, group_name, overwrite=True, verbose=True):
    """Writes or updates a dataframe into the specified HDF5 group."""
    os.makedirs(os.path.dirname(h5_path), exist_ok=True)
    mode = "a" if os.path.exists(h5_path) else "w"
    with h5py.File(h5_path, mode) as h5:
        if group_name in h5 and overwrite:
            del h5[group_name]
        grp = h5.create_group(group_name)
        for col in df.columns:
            arr = np.array(df[col].astype(str).values, dtype="S")
            grp.create_dataset(col, data=arr)
    if verbose:
        print(f"[HDF5] ✅ Wrote {group_name} ({len(df)} rows, {len(df.columns)} cols) to {h5_path}")
