"""
tas_nexus_reader.py
===================
Reader function for TAS NeXus HDF5 files written by tas_nexus_writer.py.

Returns
-------
df       : pandas DataFrame  — one row per scan point, all scalar channels
psd      : numpy array       — shape (n_points, 128, 128), float64
                               None if the dataset is all-zero / absent
scan_info: dict with keys:
    num_points      int     — total number of scan points
    scanning_axis   str     — e.g. "s2"
    start_time      str     — ISO-8601
    end_time        str     — ISO-8601
    title           str
    command         str
    filename        str
    facility        str
    source          str
    instrument_name str
    local_contact   str
    experiment_id   str
    proposal_no     str
    users           list[str]
    sample_name     str
    sample_type     str
    sample_mosaic   float   — minutes of arc
    sample_v1       np.ndarray (3,)
    sample_v2       np.ndarray (3,)
    unit_cell       np.ndarray (6,)  [a,b,c,alpha,beta,gamma]
    ub_matrix       np.ndarray (3,3)
    mono_crystal    str
    ana_crystal     str
    sense           str
    distance_vs_mono     float  m
    distance_mono_sample float  m
    distance_sample_ana  float  m
    distance_ana_det     float  m
    psd_present     bool    — True if non-zero PSD data was found
"""

import h5py
import numpy as np
import pandas as pd

# Mirrors _SCALAR_PATHS from the writer — maps DataFrame column → HDF5 path
_SCALAR_PATHS = {
    "vs_left":   "/entry/instrument/source/vs_left/value",
    "vs_right":  "/entry/instrument/source/vs_right/value",
    "m1":        "/entry/instrument/monochromator/m1/value",
    "m2":        "/entry/instrument/monochromator/m2/value",
    "monovf":    "/entry/instrument/monochromator/monovf/value",
    "monohf":    "/entry/instrument/monochromator/monohf/value",
    "monotilt":  "/entry/instrument/monochromator/monotilt/value",
    "monotrans": "/entry/instrument/monochromator/monotrans/value",
    "ei":        "/entry/instrument/monochromator/ei/value",
    "a1":        "/entry/instrument/analyzer/a1/value",
    "a2":        "/entry/instrument/analyzer/a2/value",
    "anavf":     "/entry/instrument/analyzer/anavf/value",
    "anahf":     "/entry/instrument/analyzer/anahf/value",
    "anatilt":   "/entry/instrument/analyzer/anatilt/value",
    "anatrans":  "/entry/instrument/analyzer/anatrans/value",
    "ef":        "/entry/instrument/analyzer/ef/value",
    "ps_left":   "/entry/instrument/slits/ps_left/value",
    "ps_right":  "/entry/instrument/slits/ps_right/value",
    "ps_top":    "/entry/instrument/slits/ps_top/value",
    "ps_bottom": "/entry/instrument/slits/ps_bottom/value",
    "pa_left":   "/entry/instrument/slits/pa_left/value",
    "pa_right":  "/entry/instrument/slits/pa_right/value",
    "pa_top":    "/entry/instrument/slits/pa_top/value",
    "pa_bottom": "/entry/instrument/slits/pa_bottom/value",
    "col_1_motor": "/entry/instrument/collimators/col_1_motor/value",
    "col_2_motor": "/entry/instrument/collimators/col_2_motor/value",
    "col_3_motor": "/entry/instrument/collimators/col_3_motor/value",
    "col_4_motor": "/entry/instrument/collimators/col_4_motor/value",
    "coll_alpha1": "/entry/instrument/collimators/coll_alpha1/value",
    "coll_alpha2": "/entry/instrument/collimators/coll_alpha2/value",
    "coll_alpha3": "/entry/instrument/collimators/coll_alpha3/value",
    "coll_alpha4": "/entry/instrument/collimators/coll_alpha4/value",
    "counts":    "/entry/instrument/det_group/detector/data",
    "monitor":   "/entry/instrument/monitor/data",
    "s1":  "/entry/sample/s1/value",
    "s2":  "/entry/sample/s2/value",
    "sgu": "/entry/sample/sgu/value",
    "sgl": "/entry/sample/sgl/value",
    "stu": "/entry/sample/stu/value",
    "stl": "/entry/sample/stl/value",
    "qh":  "/entry/sample/qh/value",
    "qk":  "/entry/sample/qk/value",
    "ql":  "/entry/sample/ql/value",
    "en":  "/entry/virtual_motors/en/value",
    "sample_temp1":    "/entry/sample/sample_env/sample_temp1/value",
    "sample_temp2":    "/entry/sample/sample_env/sample_temp2/value",
    "sample_temp3":    "/entry/sample/sample_env/sample_temp3/value",
    "sample_temp4":    "/entry/sample/sample_env/sample_temp4/value",
    "sample_mfield":   "/entry/sample/sample_env/sample_mfield/value",
    "sample_efield":   "/entry/sample/sample_env/sample_efield/value",
    "sample_pressure": "/entry/sample/sample_env/sample_pressure/value",
    "cryo_he_pressure":"/entry/sample/sample_env/cryo_he_pressure/value",
    "cryo_needlevalve":"/entry/sample/sample_env/cryo_needlevalve/value",
}

# HDF5 paths for scalar metadata (non-array, non-scan-axis fields)
_META_STR_PATHS = {
    "start_time":      "/entry/start_time",
    "end_time":        "/entry/end_time",
    "title":           "/entry/metadata/title",
    "command":         "/entry/metadata/command",
    "filename":        "/entry/metadata/filename",
    "facility":        "/entry/metadata/facility",
    "source":          "/entry/metadata/source",
    "instrument_name": "/entry/metadata/instrument_name",
    "local_contact":   "/entry/metadata/local_contact",
    "experiment_id":   "/entry/metadata/experiment_id",
    "proposal_no":     "/entry/metadata/proposal_no",
    "sample_name":     "/entry/metadata/sample_name",
    "sample_type":     "/entry/metadata/sample_type",
    "scanning_axis":   "/entry/metadata/scanning_axis",
    "mono_crystal":    "/entry/instrument/monochromator/mono_crystal",
    "ana_crystal":     "/entry/instrument/analyzer/ana_crystal",
    "sense":           "/entry/instrument/sense",
}

_META_FLOAT_PATHS = {
    "sample_mosaic":        "/entry/sample/sample_mosaic",
    "distance_vs_mono":     "/entry/instrument/distance_vs_mono",
    "distance_mono_sample": "/entry/instrument/distance_mono_sample",
    "distance_sample_ana":  "/entry/instrument/distance_sample_ana",
    "distance_ana_det":     "/entry/instrument/distance_ana_det",
}

_META_ARRAY_PATHS = {
    "sample_v1":  "/entry/sample/sample_v1",
    "sample_v2":  "/entry/sample/sample_v2",
    "unit_cell":  "/entry/sample/unit_cell",
    "ub_matrix":  "/entry/sample/ub_matrix",
}


def _read_str(f, path, default=""):
    """Safely read a scalar string dataset; decode bytes if needed."""
    if path not in f:
        return default
    val = f[path][()]
    if isinstance(val, bytes):
        return val.decode("utf-8")
    if isinstance(val, np.ndarray):
        # string array (e.g. users list)
        return [v.decode("utf-8") if isinstance(v, bytes) else str(v)
                for v in val.flat]
    return str(val)


def _read_float(f, path, default=0.0):
    """Safely read a scalar float dataset."""
    if path not in f:
        return default
    return float(f[path][()])


def _read_array(f, path):
    """Safely read an array dataset; returns None if absent."""
    if path not in f:
        return None
    return f[path][()].copy()


def load_from_hdf(filename, load_psd=True):
    """
    Read a TAS NeXus HDF5 file into a DataFrame, PSD array, and scan_info dict.

    Parameters
    ----------
    filename : str
        Path to the .tas.nxs.h5 file.
    load_psd : bool
        If True (default), load the full PSD array (n_points, 128, 128).
        Set to False to skip loading PSD data (faster for scalar-only work).

    Returns
    -------
    df : pandas DataFrame
        One row per scan point. Columns match the writer's _SCALAR_PATHS keys.
        Only columns that exist in the file and have non-uniform data OR are
        recognised detector/motor channels are included.

    psd : numpy array or None
        Shape (n_points, 128, 128), dtype float64.
        None if load_psd=False or the PSD dataset is absent in the file.

    scan_info : dict
        Scalar metadata and array quantities — see module docstring for keys.
    """
    with h5py.File(filename, "r") as f:

        # ── 1. Scalar scan data → DataFrame ──────────────────────────────────
        data = {}
        for col, path in _SCALAR_PATHS.items():
            if path in f:
                arr = f[path][()].astype(np.float64)
                data[col] = arr

        df = pd.DataFrame(data)
        n_points = len(df)

        # ── 2. PSD ───────────────────────────────────────────────────────────
        psd = None
        psd_path = "/entry/instrument/det_group/psd/data"
        psd_present = False
        if load_psd and psd_path in f:
            psd = f[psd_path][()].astype(np.float64)
            psd_present = bool(np.any(psd != 0))
            if not psd_present:
                psd = None   # return None for all-zero placeholder

        # ── 3. Metadata ───────────────────────────────────────────────────────
        scan_info = {}

        # string fields
        for key, path in _META_STR_PATHS.items():
            scan_info[key] = _read_str(f, path)

        # users is a string array — override the scalar read above
        users_path = "/entry/metadata/user"
        if users_path in f:
            raw = f[users_path][()]
            if raw.ndim == 0:
                # scalar stored as single string
                val = raw.item()
                scan_info["users"] = [val.decode() if isinstance(val, bytes) else str(val)]
            else:
                scan_info["users"] = [
                    v.decode("utf-8") if isinstance(v, bytes) else str(v)
                    for v in raw
                ]
        else:
            scan_info["users"] = []

        # float fields
        for key, path in _META_FLOAT_PATHS.items():
            scan_info[key] = _read_float(f, path)

        # array fields
        for key, path in _META_ARRAY_PATHS.items():
            arr = _read_array(f, path)
            scan_info[key] = arr  # None if absent

        # derived / computed fields
        scan_info["num_points"]   = n_points
        scan_info["psd_present"]  = psd_present

        # scanning_axis already read from metadata; fall back to NXdata @axes
        if not scan_info.get("scanning_axis"):
            nxdata_path = "/entry/data"
            if nxdata_path in f:
                axes_attr = f[nxdata_path].attrs.get("axes", "")
                scan_info["scanning_axis"] = (
                    axes_attr.decode() if isinstance(axes_attr, bytes)
                    else str(axes_attr)
                )

    return df, psd, scan_info


def export_to_text(df, psd, scan_info, output_filename):
    """
    Exports TAS DataFrame, metadata, and PSD array to a formatted text file.
    
    Parameters
    ----------
    df : pandas DataFrame
        DataFrame returned by load_from_hdf.
    psd : numpy array or None
        PSD array returned by load_from_hdf.
    scan_info : dict
        Metadata dictionary returned by load_from_hdf.
    output_filename : str
        Path for the destination text file.
    """
    # 1. Gather scanning axes
    scanning_axis_str = scan_info.get("scanning_axis", "")
    # Handle possible comma/space separated multiple scanning axes
    scan_axes = [ax.strip() for ax in scanning_axis_str.replace(",", " ").split() if ax.strip()]
    if not scan_axes:
        scan_axes = ["s2"] # fallback if completely empty
        
    # 2. Define the exact column order requirements
    fixed_start = ["qh", "qk", "ql", "en", "ei", "ef"]
    fixed_middle = ["m1", "m2", "s1", "s2", "a1", "a2"]
    fixed_end = ["counts", "monitor"] # counts maps to 'detector' in your header
    
    # 3. Sort the remaining columns based on instrument component order
    # component priorities matching beamline path: source -> mono -> slits/collimators -> sample -> sample_env -> ana -> slits -> detector
    component_order = ["vs_", "mono", "m1", "m2", "ei", "ps_", "pa_", "col_", "s1", "s2", "sg", "st", "qh", "qk", "ql", "en", "a1", "a2", "ana", "ef", "sample_", "cryo_", "counts", "monitor"]
    
    def get_component_priority(col_name):
        for idx, prefix in enumerate(component_order):
            if col_name.startswith(prefix):
                return idx
        return len(component_order)

    remaining_cols = [c for c in df.columns if c not in scan_axes + fixed_start + fixed_middle + fixed_end]
    remaining_cols.sort(key=get_component_priority)
    
    # 4. Assemble final column order for the text file
    final_columns = list(scan_axes) + fixed_start + fixed_middle + fixed_end + remaining_cols
    # Ensure we only include columns actually present in the DataFrame
    final_columns = [c for c in final_columns if c in df.columns]
    
    # Create the export DataFrame copy and map 'counts' to 'detector' for the header print
    export_df = df[final_columns].copy()
    export_df.insert(0, "Pt.", range(1, len(export_df) + 1))
    export_df = export_df.rename(columns={"counts": "detector"})

    with open(output_filename, "w", encoding="utf-8") as out:
        # 5. Write Metadata Header
        # Mapping requested keys to the actual keys inside your scan_info dictionary
        metadata_mapping = [
            ("TAS_NeXus_Version", "tas_nexus_version"),
            ("software_version", "software_version"),
            ("filename", "filename"),
            #
            ("facility", "facility"),
            ("source", "source"),
            ("instrument", "instrument_name"),
            ("experiment_id", "experiment_id"),
            ("proposal", "proposal_no"),
            ("user(s)", "users"),
            ("local_contact", "local_contact"),
            #
            ("mono_crystal", "mono_crystal"),
            ("ana_crystal", "ana_crystal"),
            ("sense", "sense"), 
            ("distance_vs_mono", "distance_vs_mono"),
            ("distance_mono_sample", "distance_mono_sample"),
            ("distance_sample_ana", "distance_sample_ana"),
            ("distance_ana_det", "distance_ana_det"),
            #
            ("sample_name", "sample_name"),   
            ("sample_type", "sample_type"),
            ("sample_mosaic", "sample_mosaic"),
            ("sample_v1", "sample_v1"),
            ("sample_v2", "sample_v2"),
            ("unit_cell", "unit_cell"),
            ("ub_matrix", "ub_matrix"),
            #
            ("scan_no", "scan_no"), # placeholder fallback
            ("scantitle", "title"),
            ("command", "command"),
            ("start_time", "start_time"),
            ("end_time", "end_time"),
            ("scanning_axis", "scanning_axis")
        ]

        
        for header_label, info_key in metadata_mapping:
            val = scan_info.get(info_key, "")
            if header_label == "ub_matrix":
                val = str(val).replace("\n", ", ")  # single-line format for the matrix

            out.write(f"# {header_label} = {val}\n")
            
        out.write("# \n") # spacer line
        
        # 6. Write Column Headers and Scalar Data
        # Format the DataFrame to text with tab-separation for tidy columns
        df_string = export_df.to_string(index=False, index_names=False, col_space=16, justify="right")
        # Prepend '#' to the header line
        lines = df_string.splitlines()
        lines[0] = "# " + lines[0]
        
        for line in lines:
            out.write(line + "\n")
            
        # 7. Write PSD Data Slices if available
        if psd is not None and scan_info.get("psd_present", False):
            out.write("\n################PSD Data################\n")
            for slice_idx in range(psd.shape[0]):
                out.write(f"#SLICE_{slice_idx + 1:03d}\n")
                # Save 2D matrix frame slice-by-slice
                np.savetxt(out, psd[slice_idx], fmt="%12.1f")
                out.write("\n")

# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import os, sys

    test_file = "test_batch.tas.nxs.h5"
    if not os.path.exists(test_file):
        print(f"Test file '{test_file}' not found — run tas_nexus_writer_v1.py first.")
        sys.exit(1)

    df, psd, info = load_from_hdf(test_file)

    print(f"\n{'='*55}")
    print(f"File : {test_file}")
    print(f"{'='*55}")

    print(f"\n── DataFrame ──────────────────────────────────────────")
    print(f"  shape       : {df.shape}  ({info['num_points']} points × {len(df.columns)} columns)")
    print(f"  columns     : {list(df.columns)}")
    print(f"  s2 range    : {df['s2'].min():.3f} → {df['s2'].max():.3f} degree")
    print(f"  counts peak : {df['counts'].max():.1f}  (at s2={df.loc[df['counts'].idxmax(),'s2']:.3f})")
    print(f"  first row   :\n{df.head(1).T.to_string()}")

    print(f"\n── PSD ────────────────────────────────────────────────")
    if psd is not None:
        print(f"  shape       : {psd.shape}")
        print(f"  dtype       : {psd.dtype}")
        print(f"  frame[0] sum: {psd[0].sum():.0f}")
        print(f"  total sum   : {psd.sum():.0f}")
    else:
        print("  Not present or all-zero.")

    print(f"\n── scan_info ──────────────────────────────────────────")
    for k, v in info.items():
        if isinstance(v, np.ndarray):
            print(f"  {k:25s}: array{v.shape} = {v.ravel()[:6]}{'...' if v.size > 6 else ''}")
        elif isinstance(v, list):
            print(f"  {k:25s}: {v}")
        else:
            print(f"  {k:25s}: {v}")

    # Round-trip check: values written == values read
    print(f"\n── Round-trip value checks ────────────────────────────")
    checks = {
        "s2 first"  : (df["s2"].iloc[0],   -65.0,   1e-6),
        "s2 last"   : (df["s2"].iloc[-1],  -61.0,   1e-6),
        "ei mean"   : (df["ei"].mean(),     14.6,    1e-6),
        "monitor[0]": (df["monitor"].iloc[0], 10000.0, 1e-6),
        "n_points"  : (info["num_points"],  101,     0),
        "dist_vs_mono": (info["distance_vs_mono"], 3.5, 1e-9),
        "dist_m_s"  : (info["distance_mono_sample"], 2.0, 1e-9),
        "mosaic"    : (info["sample_mosaic"], 0.3,   1e-9),
        "scanning_axis": (info["scanning_axis"], "s2", None),
        "sense"     : (info["sense"],        "+-+",  None),
        "mono_xtal" : (info["mono_crystal"], "PG",   None),
    }
    all_ok = True
    for label, (got, expected, tol) in checks.items():
        if tol is None:
            ok = str(got) == str(expected)
        elif tol == 0:
            ok = got == expected
        else:
            ok = abs(float(got) - float(expected)) <= tol
        tag = "✓" if ok else "✗"
        if not ok:
            all_ok = False
        print(f"  {tag}  {label:20s}  got={got!r}  expected={expected!r}")

    print(f"\n{'All round-trip checks passed.' if all_ok else 'Some checks FAILED.'}")

    # Example execution within the "__main__" block:
    txt_output_file = "test_batch_exported.txt"
    export_to_text(df, psd, info, txt_output_file)
    print(f"Successfully exported data to: {txt_output_file}")
