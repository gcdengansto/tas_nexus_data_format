# TAS NeXus Writer

A Python module for writing Triple-Axis Spectrometer (TAS) scan data to
NeXus-compliant HDF5 files. Two write modes are provided — a batch converter
for complete datasets and a point-by-point writer for real-time data
acquisition — both producing identical file structures.

---

## Requirements

```
h5py >= 3.0
numpy
pandas
```

Install with:

```bash
pip install h5py numpy pandas
```

---

## Quick Start

```python
import numpy as np
import pandas as pd
from tas_nexus_writer import save_to_hdf, save_to_hdf_point_by_point

# Metadata shared by both functions
meta = {
    "facility":        "ANSTO",
    "instrument_name": "TAIPAN",
    "sample_name":     "YBCO",
    "users":           ["Alice Smith", "Bob Jones"],
    "title":           "s2 scan of Bragg (0 0 1.5)",
    "experiment_id":   "exp1234",
    "proposal_no":     "P21234",
    "unit_cell":       [3.82, 3.82, 11.68, 90, 90, 90],
}

# ── Mode 1: convert a complete DataFrame ─────────────────────────────────────
save_to_hdf(df, "scan.tas.nxs.h5", psd=psd_array, meta=meta)

# ── Mode 2: append one point at a time during acquisition ────────────────────
for i, row in df.iterrows():
    save_to_hdf_point_by_point(
        row, "scan.tas.nxs.h5",
        psd_frame=psd_array[i],
        meta=meta,
        n_points=101,   # required on first call only
    )
```

---

## Functions

### `save_to_hdf(df, filename, psd=None, meta=None, overwrite=True)`

Converts a complete pandas DataFrame to a NeXus HDF5 file in a single
operation. All datasets are created at their final size before any data is
written, making this the most efficient option when the full scan is already
in memory.

| Parameter  | Type | Description |
|------------|------|-------------|
| `df` | `pd.DataFrame` | One row per scan point. Scalar motor, detector, and environment columns (see Column Reference below). No PSD column — PSD is passed separately. |
| `filename` | `str` | Output file path. Conventionally ends in `.tas.nxs.h5`. |
| `psd` | `np.ndarray` or `None` | PSD frames with shape `(n_points, 128, 128)`. Pass `None` if no PSD detector is present; the dataset will be filled with zeros. |
| `meta` | `dict` or `None` | Metadata dictionary (see Metadata Reference below). All keys are optional. |
| `overwrite` | `bool` | If `True` (default), silently replace an existing file. |

Returns the output `filename` as a string.

---

### `save_to_hdf_point_by_point(data_row, filename, psd_frame=None, meta=None, n_points=None)`

Appends a single scan point to an HDF5 file. On the first call the file is
created with all datasets pre-allocated to `n_points`. Subsequent calls write
into the next available slot. The file can be read by analysis software at
any point during the scan — it is a valid NeXus file from the moment of
creation.

| Parameter   | Type | Description |
|-------------|------|-------------|
| `data_row`  | `pd.Series` | One row of scalar values (same columns as the DataFrame used by `save_to_hdf`). |
| `filename`  | `str` | Output file path. Created on the first call; appended to on subsequent calls. |
| `psd_frame` | `np.ndarray` or `None` | Single PSD frame with shape `(128, 128)`. |
| `meta`      | `dict` or `None` | Metadata — only used when the file is first created. |
| `n_points`  | `int` | Total number of scan points expected. **Required on the first call.** Ignored on subsequent calls. |

Returns the current 1-based point index (int).

Raises `ValueError` if the file does not yet exist and `n_points` is not
provided.

---

## DataFrame Column Reference

All columns are optional — missing columns are silently filled with zeros.
PSD data is **never** a DataFrame column; it is always passed separately via
the `psd` / `psd_frame` parameter.

### Sample motors

| Column | Units | Description |
|--------|-------|-------------|
| `s1` | degree | Sample rotation (chi / omega) |
| `s2` | degree | Sample rotation (2-theta / phi) — typically the scan axis |
| `sgu` | degree | Sample goniometer upper |
| `sgl` | degree | Sample goniometer lower |
| `stu` | degree | Sample tilt upper |
| `stl` | degree | Sample tilt lower |

### Reciprocal space & energy (virtual motors)

| Column | Units | Description |
|--------|-------|-------------|
| `qh` | rlu | H reciprocal lattice coordinate |
| `qk` | rlu | K reciprocal lattice coordinate |
| `ql` | rlu | L reciprocal lattice coordinate |
| `ei` | meV | Incident energy |
| `ef` | meV | Final energy |
| `en` | meV | Energy transfer (Ei − Ef) |

### Monochromator motors

| Column | Units | Description |
|--------|-------|-------------|
| `m1` | degree | Monochromator rotation 1 |
| `m2` | degree | Monochromator rotation 2 |
| `monovf` | degree | Monochromator vertical focus |
| `monohf` | degree | Monochromator horizontal focus |
| `monotilt` | degree | Monochromator tilt |
| `monotrans` | mm | Monochromator translation |

### Analyzer motors

| Column | Units | Description |
|--------|-------|-------------|
| `a1` | degree | Analyzer rotation 1 |
| `a2` | degree | Analyzer rotation 2 |
| `anavf` | degree | Analyzer vertical focus |
| `anahf` | degree | Analyzer horizontal focus |
| `anatilt` | degree | Analyzer tilt |
| `anatrans` | mm | Analyzer translation |

### Slits

| Column | Units |
|--------|-------|
| `ps_left`, `ps_right`, `ps_top`, `ps_bottom` | mm |
| `pa_left`, `pa_right`, `pa_top`, `pa_bottom` | mm |

### Collimators

| Column | Units |
|--------|-------|
| `col_1_motor` … `col_4_motor` | arcmin |
| `col_alpha1` … `col_alpha4`   | min |

### Source

| Column | Units |
|--------|-------|
| `vs_left`, `vs_right` | mm |

### Detectors

| Column | Units | Description |
|--------|-------|-------------|
| `counts` | counts | Single-detector integrated counts per point |
| `monitor` | counts | Monitor counts per point |

### Sample environment

| Column | Units | Description |
|--------|-------|-------------|
| `sample_temp1` … `sample_temp4` | K | Temperature sensors |
| `sample_mfield` | T | Magnetic field |
| `sample_efield` | V | Electric field |
| `sample_pressure` | GPa | Pressure |
| `cryo_he_pressure` | mBar | Cryostat He pressure |
| `cryo_needlevalve` | % | Cryostat needle valve position |

---

## Metadata Reference

All keys are optional. Sensible defaults are used when a key is absent.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `facility` | str | `"Facility_Name"` | Facility name |
| `source` | str | `"Reactor or Spallation Source"` | Source name |
| `instrument_name` | str | `"TAS"` | Instrument name |
| `software_version` | str | `"0.9.1"` | DAQ software version |
| `tas_nexus_version` | str | `"0.9.1"` | NeXus writer version |
| `experiment_id` | str | `""` | Experiment identifier |
| `proposal_no` | str | `""` | Proposal number |
| `users` | list[str] | `["Unknown"]` | List of user names |
| `local_contact` list[str] | `["Unknown"]` |  local contact scientist name |
| `Monochromator` | str | `PG` | Monochromator Crystal |
| `Aanlyzer` | str | `PG` | Aanlyzer Crystal |
| `distance_vs_mono` | float | `2.0` | Monochromator–sample distance (m) |
| `distance_mono_sample` | float | `2.0` | Monochromator–sample distance (m) |
| `distance_sample_ana` | float | `1.5` | Sample–analyzer distance (m) |
| `distance_ana_det` | float | `0.5` | Analyzer–detector distance (m) |

| `sample_name` | str | `""` | Sample name |
| `sample_type` | str | `""` | Sample type: crystal or powder |
| `sample_mosaic` | float | `""` | Sample mosaic |
| `sample_v1` | list[float] | `[1 0 0]` | Sample vector 1 |
| `sample_v2` | list[float] | `[0 0 1]` | Sample vector 2 |
| `unit_cell` | list[float] | YBCO defaults | `[a, b, c, α, β, γ]` in Å and degrees |
| `ub_matrix` | array-like (3×3) | identity-like | UB orientation matrix (HKL → Q) |

| `title` | str | `""` | Scan title |
| `command` | str | `""` | Scan command |
| `start_time` | str | current UTC | ISO-8601 start timestamp |
| `end_time` | str | current UTC | ISO-8601 end timestamp |

| `scanning_axis` | str | `""` | Scan command |



---

## NeXus File Structure

The file follows the `NXtas` application definition. The tree below shows
every group and the most important datasets. Scan-axis datasets (marked `[n]`)
have length equal to the number of scan points. PSD data (marked `[n,128,128]`)
has one 128×128 frame per point.

```
/ (root)
│   @default = "entry"
│
└── entry/                          NXentry
    │   @default = "data"
    │   @title
    │   definition = "NXtas"
    │   start_time, end_time
    │
    ├── instrument/                 NXinstrument
    │   ├── source/                 NXsource
    │   │   ├── vs_left/value       [n]  mm
    │   │   └── vs_right/value      [n]  mm
    │   │
    │   ├── monochromator/          NXmonochromator
    │   │   ├── distance                 m
    │   │   ├── m1, m2, monovf,
    │   │   │   monohf, monotilt,
    │   │   │   monotrans / value   [n]  degree / mm
    │   │   └── ei / value          [n]  meV
    │   │
    │   ├── analyzer/               NXmonochromator
    │   │   ├── distance                 m
    │   │   ├── a1, a2, anavf,
    │   │   │   anahf, anatilt,
    │   │   │   anatrans / value    [n]  degree / mm
    │   │   └── ef / value          [n]  meV
    │   │
    │   ├── det_group/              NXdetector_group
    │   │   ├── detector/           NXdetector
    │   │   │   ├── data            [n]  counts   ← single detector
    │   │   │   └── distance             m
    │   │   └── psd/                NXdetector
    │   │       ├── data            [n, 128, 128]  counts  ← PSD
    │   │       ├── x_pixel_offset  [128]  pixel
    │   │       └── y_pixel_offset  [128]  pixel
    │   │
    │   ├── monitor/                NXmonitor
    │   │   └── data                [n]  counts
    │   │
    │   ├── slits/                  NXcollection
    │   │   └── ps_left … pa_bottom / value  [n]  mm
    │   │
    │   ├── collimators/            NXcollection
    │   │   └── col_1_motor … col_4_motor / value  [n]  arcmin
    │   │
    │   └── (filters, attenuator, misc_motors — static)
    │
    ├── sample/                     NXsample
    │   ├── distance                     m
    │   ├── unit_cell                    [a,b,c,α,β,γ]
    │   ├── ub_matrix                    3×3
    │   ├── s1, s2, sgu, sgl,
    │   │   stu, stl / value        [n]  degree
    │   ├── qh, qk, ql / value      [n]  rlu
    │   └── sample_env/             NXenvironment
    │       ├── name, type, description
    │       ├── sample_temp1…4 / value  [n]  K
    │       ├── sample_mfield / value   [n]  T
    │       ├── sample_efield / value   [n]  V
    │       ├── sample_pressure / value [n]  GPa
    │       ├── cryo_he_pressure / value [n] mBar
    │       └── cryo_needlevalve / value [n] %
    │
    ├── virtual_motors/             NXcollection
    │   ├── qh → sample/qh          (hard link)
    │   ├── qk → sample/qk          (hard link)
    │   ├── ql → sample/ql          (hard link)
    │   ├── ei → monochromator/ei   (hard link)
    │   ├── ef → analyzer/ef        (hard link)
    │   └── en / value              [n]  meV
    │
    ├── metadata/                   NXcollection
    │   ├── facility, instrument_name
    │   ├── software version, TAS_NeXus_Version
    │   ├── experiment_id, proposal_no
    │   ├── sample_name, user[], title
    │   └── filename
    │
    └── data/                       NXdata  ← default plot view
        │   @signal = "data"
        │   @axes  = "s2"
        │   @interpretation = "spectrum"
        ├── data  → detector/data   (hard link)
        ├── monitor → monitor/data  (hard link)
        ├── s2  → sample/s2/value   (hard link)
        ├── qh, qk, ql              (hard links via virtual_motors)
        └── en, ei, ef              (hard links via virtual_motors)
```

### Key design points

**Virtual motors** — `qh`, `qk`, `ql`, `ei`, and `ef` are stored once under
`sample/` or `instrument/` and exposed again under `virtual_motors/` via HDF5
hard links. The `NXdata` group links to the `virtual_motors/` copies so that
reciprocal-space and energy axes are immediately available for plotting without
traversing the instrument tree.

**PSD shape** — The position-sensitive detector dataset has shape
`(n_points, 128, 128)` because a full 2D frame is collected at every motor
position during a scan. This is physically correct and distinct from a
single-frame `(128, 128)` snapshot.

**Default plot view** — The `NXdata` group at `/entry/data` is the NeXus
default plot target (`@default = "data"` on the entry). It links the single
detector counts as the signal and `s2` as the primary axis, which allows
NeXus-aware viewers (e.g. NeXpy, Mantid, HDFView with NeXus plugin) to
display the scan immediately on file open.

**Identical structure from both write modes** — `save_to_hdf` and
`save_to_hdf_point_by_point` call the same internal `_build_skeleton` and
`_write_scalars` / `_write_psd` helpers, guaranteeing that files produced
during live acquisition are structurally identical to files converted from a
complete DataFrame after the fact.
