"""
tas_nexus_writer.py
===================
Two functions for writing TAS (Triple-Axis Spectrometer) scan data
to NeXus-compliant HDF5 files, matching the structure of v2e_psd.py.

Assumptions
-----------
- The pandas DataFrame / Series has one row per scan point for all
  scalar motor/detector/environment columns.
- PSD data is passed SEPARATELY as a numpy array:
    - save_to_hdf            : shape (n_points, 128, 128)
    - save_to_hdf_point_by_point : shape (128, 128), one frame per call
- save_to_hdf_point_by_point receives n_points upfront at file creation
  so datasets are pre-allocated (no dynamic resize needed).
- Metadata is passed as a separate dict to both functions.
- Both functions produce identical HDF5 structure.

DataFrame / Series column layout
----------------------------------
Motors (sample):   s1, s2, sgu, sgl, stu, stl
Virtual (Q/E):     qh, qk, ql, ei, ef, en
Mono motors:       m1, m2, monovf, monohf, monotilt, monotrans
Ana motors:        a1, a2, anavf, anahf, anatilt, anatrans
Source:            vs_left, vs_right
Slits:             ps_left, ps_right, ps_top, ps_bottom,
                   pa_left, pa_right, pa_top, pa_bottom
Collimators:       col_1_motor, col_2_motor, col_3_motor, col_4_motor
Sample env:        sample_temp1..4, sample_mfield, sample_efield,
                   sample_pressure, cryo_he_pressure, cryo_needlevalve
Detectors:         counts, monitor

Metadata dict keys (all optional, sensible defaults used if absent)
--------------------------------------------------------------------
facility, instrument_name, software_version, tas_nexus_version,
experiment_id, proposal_no, sample_name, users (list[str]),
title, start_time, end_time,
unit_cell     (list of 6 floats [a,b,c,alpha,beta,gamma]),
ub_matrix     (3×3 array-like),
distance_mono_sample (float, metres),
distance_sample_ana  (float, metres),
distance_ana_det     (float, metres)
"""

import os
from datetime import datetime, timezone

import h5py
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Constants & helpers
# ---------------------------------------------------------------------------

_STR_DT = h5py.string_dtype("utf-8")

_MOTOR_UNITS = {
    "s1": "degree",  "s2": "degree",
    "sgu": "degree", "sgl": "degree",
    "stu": "degree", "stl": "degree",
    "qh": "rlu",     "qk": "rlu",     "ql": "rlu",
    "ei": "meV",     "ef": "meV",     "en": "meV",
    "m1": "degree",  "m2": "degree",
    "monovf": "degree", "monohf": "degree",
    "monotilt": "degree", "monotrans": "mm",
    "a1": "degree",  "a2": "degree",
    "anavf": "degree", "anahf": "degree",
    "anatilt": "degree", "anatrans": "mm",
    "vs_left": "mm",   "vs_right": "mm",
    "ps_left": "mm",   "ps_right": "mm",
    "ps_top": "mm",    "ps_bottom": "mm",
    "pa_left": "mm",   "pa_right": "mm",
    "pa_top": "mm",    "pa_bottom": "mm",
    "col_1_motor": "arcmin", "col_2_motor": "arcmin",
    "col_3_motor": "arcmin", "col_4_motor": "arcmin",
    "coll_alpha1": "arcmin", "coll_alpha2": "arcmin",
    "coll_alpha3": "arcmin", "coll_alpha4": "arcmin",
    "sample_temp1": "K",   "sample_temp2": "K",
    "sample_temp3": "K",   "sample_temp4": "K",
    "sample_mfield": "T",  "sample_efield": "V",
    "sample_pressure": "GPa",
    "cryo_he_pressure": "mBar", "cryo_needlevalve": "%",
    "counts": "counts",    "monitor": "counts",
}

# HDF5 path of the "value" (or "data") dataset for every scalar column
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


def _now_iso():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _meta(meta, key, default=""):
    return meta.get(key, default) if meta else default


def _col_arr(source, col, n, dtype=np.float64):
    """Extract a column from a DataFrame or Series as a 1-D numpy array.
    Returns zeros of length n if the column is absent."""
    if isinstance(source, pd.Series):
        return np.array([source[col]], dtype=dtype) if col in source.index \
               else np.zeros(1, dtype=dtype)
    return source[col].to_numpy(dtype=dtype) if col in source.columns \
           else np.zeros(n, dtype=dtype)


# ---------------------------------------------------------------------------
# Skeleton builder  (called once when a new file is created)
# ---------------------------------------------------------------------------

def _build_skeleton(f, meta, n_points):
    """
    Create every group, static dataset, and pre-allocated scan dataset
    for a file that will hold exactly n_points scan points.

    All resizable/scan datasets are created with shape=(n_points,) and
    filled with zeros; actual data is written afterwards by the callers.
    """
    n   = n_points
    str_dt = _STR_DT

    def positioner(parent, name):
        """Create an NXpositioner group with a pre-allocated value array."""
        g = parent.require_group(name)
        g.attrs["NX_class"] = "NXpositioner"
        g.create_dataset("value", shape=(n,), dtype=np.float64)
        g["value"].attrs["units"] = _MOTOR_UNITS.get(name, "degree")
        return g

    def sensor(parent, name, unit):
        """Create an NXsensor group with a pre-allocated value array."""
        g = parent.require_group(name)
        g.attrs["NX_class"] = "NXsensor"
        g.create_dataset("value", shape=(n,), dtype=np.float64)
        g["value"].attrs["units"] = unit
        return g

    # ── /entry ──────────────────────────────────────────────────────────────
    entry = f.require_group("entry")
    entry.attrs["NX_class"] = "NXentry"
    entry.attrs["title"]    = _meta(meta, "title", "TAS scan")
    entry.create_dataset("definition", data="NXtas",                    dtype=str_dt)
    entry.create_dataset("start_time", data=_meta(meta,"start_time",_now_iso()), dtype=str_dt)
    entry.create_dataset("end_time",   data=_meta(meta,"end_time",""),  dtype=str_dt)

    # ── instrument ──────────────────────────────────────────────────────────
    inst = entry.require_group("instrument")
    inst.attrs["NX_class"] = "NXinstrument"

    # source
    src = inst.require_group("source")
    src.attrs["NX_class"] = "NXsource"
    for name in ["vs_left", "vs_right"]:
        positioner(src, name)

    # monochromator
    mono = inst.require_group("monochromator")
    mono.attrs["NX_class"] = "NXmonochromator"
    mono.create_dataset("distance", data=float(_meta(meta,"distance_mono_sample",2.0)))
    mono["distance"].attrs["units"] = "m"
    mono.create_dataset("mono_crystal", data=_meta(meta,"mono_crystal","PG"), dtype=str_dt)
    for name in ["m1","m2","monovf","monohf","monotilt","monotrans","ei"]:
        positioner(mono, name)

    # analyzer
    ana = inst.require_group("analyzer")
    ana.attrs["NX_class"] = "NXmonochromator"
    ana.create_dataset("distance", data=float(_meta(meta,"distance_sample_ana",1.5)))
    ana["distance"].attrs["units"] = "m"
    ana.create_dataset("ana_crystal", data=_meta(meta,"ana_crystal","PG"), dtype=str_dt)  # NEW
    for name in ["a1","a2","anavf","anahf","anatilt","anatrans","ef"]:
        positioner(ana, name)

    # det_group → single detector
    dg = inst.require_group("det_group")
    dg.attrs["NX_class"] = "NXdetector_group"
    sd = dg.require_group("detector")
    sd.attrs["NX_class"] = "NXdetector"
    sd.create_dataset("data", shape=(n,), dtype=np.float64)
    sd["data"].attrs["units"] = "counts"
    sd.create_dataset("distance", data=float(_meta(meta,"distance_ana_det",0.5)))
    sd["distance"].attrs["units"] = "m"

    # det_group → PSD  — shape (n_points, 128, 128), passed separately
    psd = dg.require_group("psd")
    psd.attrs["NX_class"] = "NXdetector"
    psd.create_dataset("data", shape=(n, 128, 128), dtype=np.float64)
    psd["data"].attrs["units"] = "counts"
    psd.create_dataset("x_pixel_offset", data=np.arange(128))
    psd["x_pixel_offset"].attrs["units"] = "pixel"
    psd.create_dataset("y_pixel_offset", data=np.arange(128))
    psd["y_pixel_offset"].attrs["units"] = "pixel"

    # monitor
    mon = inst.require_group("monitor")
    mon.attrs["NX_class"] = "NXmonitor"
    mon.create_dataset("data", shape=(n,), dtype=np.float64)
    mon["data"].attrs["units"] = "counts"

    # slits
    slits = inst.require_group("slits")
    slits.attrs["NX_class"] = "NXcollection"
    for name in ["ps_left","ps_right","ps_top","ps_bottom",
                 "pa_left","pa_right","pa_top","pa_bottom"]:
        positioner(slits, name)

    # collimators
    colli = inst.require_group("collimators")
    colli.attrs["NX_class"] = "NXcollection"
    for name in ["col_1_motor","col_2_motor","col_3_motor","col_4_motor"]:
        positioner(colli, name)
    for name in ["coll_alpha1","coll_alpha2","coll_alpha3","coll_alpha4"]:
        positioner(colli, name)
    
    # --- instrument-level geometry & configuration ---
    inst.create_dataset("distance_vs_mono",     data=float(_meta(meta,"distance_vs_mono",    0.0)))
    inst["distance_vs_mono"].attrs["units"]     = "m"
    inst.create_dataset("distance_mono_sample", data=float(_meta(meta,"distance_mono_sample",2.0)))
    inst["distance_mono_sample"].attrs["units"] = "m"
    inst.create_dataset("distance_sample_ana",  data=float(_meta(meta,"distance_sample_ana", 1.5)))
    inst["distance_sample_ana"].attrs["units"]  = "m"
    inst.create_dataset("distance_ana_det",     data=float(_meta(meta,"distance_ana_det",    0.5)))
    inst["distance_ana_det"].attrs["units"]     = "m"
    inst.create_dataset("sense", data=_meta(meta,"sense","+-+"), dtype=str_dt)

    # ── sample ──────────────────────────────────────────────────────────────
    sample = entry.require_group("sample")
    sample.attrs["NX_class"] = "NXsample"
    sample.create_dataset("distance", data=float(_meta(meta,"distance_mono_sample",1.5)))
    sample["distance"].attrs["units"] = "m"
    uc = _meta(meta, "unit_cell", [3.82, 3.82, 11.68, 90, 90, 90])
    sample.create_dataset("unit_cell", data=np.array(uc, dtype=np.float64))
    ub = _meta(meta, "ub_matrix",
               [[0.1254,0.0021,0.0],[0.0018,0.1189,0.0],[0.0,0.0,0.0765]])
    ub_ds = sample.create_dataset("ub_matrix", data=np.array(ub, dtype=np.float64))
    ub_ds.attrs["description"]   = "Orientation matrix transforming HKL to Q"
    ub_ds.attrs["interpretation"] = "matrix"

     # --- NEW: sample descriptors ---
    sample.create_dataset("sample_mosaic", data=float(_meta(meta,"sample_mosaic", 0.0)))
    sample["sample_mosaic"].attrs["units"] = "minutes of arc"
    sample.create_dataset("sample_v1", data=np.array(_meta(meta,"sample_v1",[1,0,0]), dtype=np.float64))
    sample.create_dataset("sample_v2", data=np.array(_meta(meta,"sample_v2",[0,1,0]), dtype=np.float64))

    for name in ["s1","s2","sgu","sgl","stu","stl"]:
        positioner(sample, name)
    for name in ["qh","qk","ql"]:
        positioner(sample, name)

    # sample environment
    se = sample.require_group("sample_env")
    se.attrs["NX_class"] = "NXenvironment"
    se.create_dataset("name",        data="Cryostat CF16",       dtype=str_dt)
    se.create_dataset("type",        data="Cryostat",            dtype=str_dt)
    se.create_dataset("description", data="Closed-cycle cryostat", dtype=str_dt)
    for name, unit in [
        ("sample_temp1","K"), ("sample_temp2","K"),
        ("sample_temp3","K"), ("sample_temp4","K"),
        ("sample_mfield","T"), ("sample_efield","V"),
        ("sample_pressure","GPa"),
        ("cryo_he_pressure","mBar"), ("cryo_needlevalve","%"),
    ]:
        sensor(se, name, unit)

    # ── virtual motors ───────────────────────────────────────────────────────
    vm = entry.require_group("virtual_motors")
    vm.attrs["NX_class"] = "NXcollection"
    vm["qh"] = sample["qh"]          # hard links into sample
    vm["qk"] = sample["qk"]
    vm["ql"] = sample["ql"]
    vm["ei"] = inst["monochromator/ei"]
    vm["ef"] = inst["analyzer/ef"]
    en_g = vm.require_group("en")
    en_g.attrs["NX_class"] = "NXpositioner"
    en_g.create_dataset("value", shape=(n,), dtype=np.float64)
    en_g["value"].attrs["units"] = "meV"


    # ── metadata ─────────────────────────────────────────────────────────────
    md = entry.require_group("metadata")
    md.create_dataset("facility",            data=_meta(meta,"facility","Facility_Name"),        dtype=str_dt)
    md.create_dataset("source",              data=_meta(meta,"source",""),                       dtype=str_dt)  # e.g. "HIFR Reactor"
    md.create_dataset("instrument_name",     data=_meta(meta,"instrument_name","TAS"),           dtype=str_dt)
    md.create_dataset("software version",    data=_meta(meta,"software_version","0.9.1"),        dtype=str_dt)
    md.create_dataset("TAS_NeXus_Version",   data=_meta(meta,"tas_nexus_version","0.9.1"),       dtype=str_dt)
    md.create_dataset("experiment_id",       data=_meta(meta,"experiment_id",""),                dtype=str_dt)
    md.create_dataset("proposal_no",         data=_meta(meta,"proposal_no",""),                  dtype=str_dt)
    md.create_dataset("local_contact",       data=_meta(meta,"local_contact",""),                dtype=str_dt)
    md.create_dataset("sample_name",         data=_meta(meta,"sample_name",""),                  dtype=str_dt)
    md.create_dataset("sample_type",         data=_meta(meta,"sample_type","crystal"),           dtype=str_dt)

    md.create_dataset("user",                data=_meta(meta,"users",["Unknown"]),               dtype=str_dt)
    md.create_dataset("title",               data=_meta(meta,"title",""),                        dtype=str_dt)
    md.create_dataset("command",             data=_meta(meta,"command",""),                      dtype=str_dt)
    md.create_dataset("filename",            data=_meta(meta,"filename",""),                     dtype=str_dt)

    md.create_dataset("scanning_axis",       data=_meta(meta,"scanning_axis","s2"),                dtype=str_dt)

    # ── NXdata (default plot view — single detector 1-D) ────────────────────
    nxdata = entry.require_group("data")
    nxdata.attrs["NX_class"]       = "NXdata"
    nxdata.attrs["signal"]         = "data"
    nxdata.attrs["axes"]           = "s2"
    nxdata.attrs["interpretation"] = "spectrum"
    nxdata.attrs["s2_indices"]     = 0
    nxdata.attrs["en_indices"]     = 0
    nxdata.attrs["qh_indices"]     = 0

    nxdata["data"]    = f["/entry/instrument/det_group/detector/data"]
    nxdata["monitor"] = f["/entry/instrument/monitor/data"]
    nxdata["s2"]  = f["/entry/sample/s2/value"]
    nxdata["qh"]  = f["/entry/virtual_motors/qh/value"]
    nxdata["qk"]  = f["/entry/virtual_motors/qk/value"]
    nxdata["ql"]  = f["/entry/virtual_motors/ql/value"]
    nxdata["en"]  = f["/entry/virtual_motors/en/value"]
    nxdata["ei"]  = f["/entry/virtual_motors/ei/value"]
    nxdata["ef"]  = f["/entry/virtual_motors/ef/value"]
    nxdata["qh"].attrs["long_name"] = "QH Reciprocal Coordinate"
    nxdata["qk"].attrs["long_name"] = "QK Reciprocal Coordinate"
    nxdata["ql"].attrs["long_name"] = "QL Reciprocal Coordinate"
    nxdata["en"].attrs["long_name"] = "Energy Transfer"
    nxdata["ei"].attrs["long_name"] = "Incident Energy"
    nxdata["ef"].attrs["long_name"] = "Final Energy"

    f.attrs["default"]     = "entry"
    entry.attrs["default"] = "data"


# ---------------------------------------------------------------------------
# Internal write helpers
# ---------------------------------------------------------------------------

def _write_scalars(f, source, start, n):
    """Write all scalar columns from source (DataFrame or Series) at [start:start+n]."""
    end = start + n
    for col, path in _SCALAR_PATHS.items():
        f[path][start:end] = _col_arr(source, col, n)


def _write_psd(f, psd_data, start, n):
    """
    Write PSD frames into /entry/instrument/det_group/psd/data.

    psd_data : numpy array, shape (n, 128, 128)  — for batch
                                  (128, 128)      — for a single point
                             or None (skipped, zeros remain)
    """
    if psd_data is None:
        return
    arr = np.asarray(psd_data, dtype=np.float64)
    if arr.ndim == 2:                     # single point → add scan axis
        arr = arr[np.newaxis, ...]        # (1, 128, 128)
    f["/entry/instrument/det_group/psd/data"][start:start + n] = arr


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def save_to_hdf(df, filename, psd=None, meta=None, overwrite=True):
    """
    Convert an entire pandas DataFrame to a NeXus HDF5 file in one shot.

    Parameters
    ----------
    df        : pandas DataFrame — one row per scan point, scalar columns only.
    filename  : output .h5 filepath.
    psd       : numpy array, shape (n_points, 128, 128), or None.
                PSD frames for every scan point, passed separately from df.
    meta      : dict of metadata (see module docstring).  Optional.
    overwrite : if True (default), silently replace an existing file.

    Returns
    -------
    filename (str)
    """
    if overwrite and os.path.exists(filename):
        os.remove(filename)

    n = len(df)
    h5py.get_config().track_order = True
    with h5py.File(filename, "w", track_order=True) as f:
        _build_skeleton(f, meta, n_points=n)
        _write_scalars(f, df, start=0, n=n)
        _write_psd(f, psd, start=0, n=n)

        # Stamp end time if not supplied in metadata
        if not _meta(meta, "end_time"):
            f["/entry/end_time"][()] = _now_iso().encode()

    print(f"[save_to_hdf] Wrote {n} points → {filename}")
    return filename


def save_to_hdf_point_by_point(data_row, filename, psd_frame=None,
                                meta=None, n_points=None):
    """
    Append a single new scan point to a NeXus HDF5 file.

    On the first call (file does not yet exist) the file is created with
    all datasets pre-allocated to n_points.  Subsequent calls write into
    the next available slot without reopening the structure.

    Parameters
    ----------
    data_row  : pandas Series — one row of scalar motor/detector values.
    filename  : output .h5 filepath (created on first call, appended after).
    psd_frame : numpy array, shape (128, 128), or None.
                The PSD frame for this single scan point.
    meta      : dict of metadata — only used when creating the file.
    n_points  : total number of scan points expected in the scan.
                REQUIRED on the first call (file creation); ignored after.

    Returns
    -------
    current_point (int) — 1-based index of the point just written.

    Raises
    ------
    ValueError  if the file does not exist and n_points is not provided.
    """
    is_new = not os.path.exists(filename)
    h5py.get_config().track_order = True

    if is_new:
        if n_points is None:
            raise ValueError(
                "n_points must be provided when creating a new file."
            )
        with h5py.File(filename, "w", track_order=True) as f:
            _build_skeleton(f, meta, n_points=n_points)

    with h5py.File(filename, "a", track_order=True) as f:
        # Use the detector dataset length as the authoritative written count
        n_written = int(
            np.count_nonzero(
                np.ones(f["/entry/instrument/det_group/detector/data"].shape[0])
            )
        )
        # A cleaner index: track via a small scalar attribute
        idx = f["entry"].attrs.get("_next_point_index", 0)

        _write_scalars(f, data_row, start=idx, n=1)
        _write_psd(f, psd_frame, start=idx, n=1)

        f["entry"].attrs["_next_point_index"] = idx + 1
        current_point = idx + 1

    print(f"[save_to_hdf_point_by_point] Point {current_point}/{n_points or '?'} → {filename}")
    return current_point


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    rng = np.random.default_rng(42)
    n   = 101
    motor_pos = np.linspace(-65, -61, n)

    # ── Scalar DataFrame (no psd_data column) ───────────────────────────────
    df = pd.DataFrame({
        "s1": np.zeros(n),      "s2": motor_pos,
        "sgu": np.zeros(n),     "sgl": np.zeros(n),
        "stu": np.zeros(n),     "stl": np.zeros(n),
        "qh": np.zeros(n),      "qk": np.zeros(n),
        "ql": np.linspace(1.4, 1.6, n),
        "ei": np.full(n, 14.6), "ef": np.full(n, 14.6), "en": np.zeros(n),
        "m1": np.full(n,45.0),  "m2": np.full(n,45.0),
        "monovf": np.full(n,45.0), "monohf": np.full(n,45.0),
        "monotilt": np.full(n,45.0), "monotrans": np.full(n,45.0),
        "a1": np.full(n,60.0),  "a2": np.full(n,60.0),
        "anavf": np.full(n,60.0), "anahf": np.full(n,60.0),
        "anatilt": np.full(n,60.0), "anatrans": np.full(n,60.0),
        "vs_left": np.full(n,10.0), "vs_right": np.full(n,10.0),
        "ps_left": np.full(n,10.0), "ps_right": np.full(n,10.0),
        "ps_top":  np.full(n,10.0), "ps_bottom": np.full(n,10.0),
        "pa_left": np.full(n,10.0), "pa_right":  np.full(n,10.0),
        "pa_top":  np.full(n,10.0), "pa_bottom": np.full(n,10.0),
        "col_1_motor": np.ones(n),  "col_2_motor": np.ones(n),
        "col_3_motor": np.ones(n),  "col_4_motor": np.ones(n),
        "sample_temp1": np.full(n,295.0), "sample_temp2": np.full(n,295.0),
        "sample_temp3": np.full(n,295.0), "sample_temp4": np.full(n,295.0),
        "sample_mfield": np.zeros(n),     "sample_efield": np.zeros(n),
        "sample_pressure": np.zeros(n),
        "cryo_he_pressure": np.full(n,10.0),
        "cryo_needlevalve": np.full(n,7.0),
        "counts": (1000*np.exp(-np.power(motor_pos+63,2)/(2*0.5**2))
                   + rng.normal(0,20,n)),
        "monitor": np.full(n,10000.0),
    })

    # ── PSD array passed separately: (n_points, 128, 128) ───────────────────
    psd_all = rng.poisson(100, (n, 128, 128)).astype(np.float64)

    meta = {
    "facility":             "ANSTO",
    "source":               "OPAL Reactor",
    "instrument_name":      "TAIPAN",
    "experiment_id":        "exp1234",
    "proposal_no":          "P21234",
    "users":                ["Andrew Brown", "Alex Green"],
    "local_contact":        "Dr. Jane Doe",

    "mono_crystal":         "PG",
    "ana_crystal":          "PG",
    "sense":                "+-+",
    "distance_vs_mono":     3.5,       # NEW — source-to-mono distance
    "distance_mono_sample": 2.0,
    "distance_sample_ana":  1.5,
    "distance_ana_det":     0.5,

    "title":                "s2 scan of Bragg (0 0 1.5)",
    "command":              "scan s2 -65 -61 0.04 mon 10000",
    "filename":             "TAIPAN_#1001235.tas.nxs.h5",

    "sample_name":          "YBCO",
    "sample_type":          "crystal",
    "sample_mosaic":        0.3,
    "sample_v1":            [1, 0, 0],
    "sample_v2":            [0, 0, 1],
    "unit_cell":            [3.82, 3.82, 11.68, 90, 90, 90],
    "ub_matrix":  [[0.1254,0.0021,0.0],[0.0018,0.1189,0.0],[0.0,0.0,0.0765]],

    "scanning_axis":        "s2",
}

    # Test 1: batch
    f1 = "test_batch.tas.nxs.h5"
    save_to_hdf(df, f1, psd=psd_all, meta=meta)

    # Test 2: point-by-point (n_points known upfront)
    f2 = "test_pbp.tas.nxs.h5"
    if os.path.exists(f2):
        os.remove(f2)
    for i, row in df.iterrows():
        save_to_hdf_point_by_point(
            row, f2,
            psd_frame=psd_all[i],
            meta=meta,
            n_points=n,          # only used on first call
        )

    # Verify identical shapes and data agreement
    print("\nShape & data comparison (batch vs point-by-point):")
    check_paths = [
        "/entry/instrument/det_group/detector/data",
        "/entry/instrument/det_group/psd/data",
        "/entry/sample/s2/value",
        "/entry/virtual_motors/en/value",
        "/entry/virtual_motors/qh/value",
    ]
    all_ok = True
    with h5py.File(f1,"r") as a, h5py.File(f2,"r") as b:
        for p in check_paths:
            sa, sb  = a[p].shape, b[p].shape
            match_s = sa == sb
            match_d = np.allclose(a[p][()], b[p][()]) if match_s else False
            tag     = "✓" if (match_s and match_d) else "✗"
            print(f"  {tag}  {p}  shapes {sa} | data_equal={match_d}")
            if not (match_s and match_d):
                all_ok = False
    print("\nAll checks passed." if all_ok else "\nSome checks FAILED.")
