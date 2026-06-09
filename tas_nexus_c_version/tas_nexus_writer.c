/* Request POSIX.1-2008 extensions (gmtime_r) before any system header */
#define _POSIX_C_SOURCE 200809L

/**
 * @file tas_nexus_writer.c
 * @brief TAS NeXus HDF5 writer — C implementation
 *
 * Mirrors tas_nexus_writer_v1.py structure exactly:
 *   /entry/
 *     instrument/
 *       source/        vs_left, vs_right
 *       monochromator/ m1, m2, monovf, monohf, monotilt, monotrans, ei
 *       analyzer/      a1, a2, anavf, anahf, anatilt, anatrans, ef
 *       det_group/
 *         detector/    data[n], distance
 *         psd/         data[n,128,128], x_pixel_offset, y_pixel_offset
 *       monitor/       data[n]
 *       slits/         ps_*, pa_*
 *       collimators/   col_*_motor, coll_alpha*
 *       distance_*, sense
 *     sample/
 *       s1, s2, sgu, sgl, stu, stl
 *       qh, qk, ql
 *       unit_cell, ub_matrix, sample_mosaic, sample_v1, sample_v2
 *       sample_env/  temp1..4, mfield, efield, pressure, cryo_*
 *     virtual_motors/  qh, qk, ql (hard-links), ei, ef (hard-links), en
 *     metadata/        facility, instrument_name, …
 *     data/            NXdata signal=data axes=s2
 */

#include "tas_nexus_writer.h"

#include <hdf5.h>

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* =========================================================================
 * Compile-time settings
 * ========================================================================= */

#ifndef TAS_ERROR_BUF_SIZE
#  define TAS_ERROR_BUF_SIZE 512
#endif

/* =========================================================================
 * Thread-local error buffer
 * ========================================================================= */

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#  define TAS_THREAD_LOCAL _Thread_local
#elif defined(__GNUC__)
#  define TAS_THREAD_LOCAL __thread
#elif defined(_MSC_VER)
#  define TAS_THREAD_LOCAL __declspec(thread)
#else
#  define TAS_THREAD_LOCAL  /* best-effort: shared, not thread-safe */
#endif

static TAS_THREAD_LOCAL char s_errbuf[TAS_ERROR_BUF_SIZE];

static void set_error(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(s_errbuf, sizeof(s_errbuf), fmt, ap);
    va_end(ap);
}

const char *tas_last_error(void) { return s_errbuf; }

/* =========================================================================
 * Column metadata tables  (must stay in sync with tas_col_id_t enum)
 * ========================================================================= */

typedef struct {
    const char *name;          /* column name, e.g. "s2"            */
    const char *hdf5_path;     /* absolute path of the value dataset */
    const char *unit;          /* physical unit string               */
} col_meta_t;

static const col_meta_t COL_META[TAS_N_COLS] = {
    /* Sample */
    [TAS_COL_S1]   = {"s1",  "/entry/sample/s1/value",  "degree"},
    [TAS_COL_S2]   = {"s2",  "/entry/sample/s2/value",  "degree"},
    [TAS_COL_SGU]  = {"sgu", "/entry/sample/sgu/value", "degree"},
    [TAS_COL_SGL]  = {"sgl", "/entry/sample/sgl/value", "degree"},
    [TAS_COL_STU]  = {"stu", "/entry/sample/stu/value", "degree"},
    [TAS_COL_STL]  = {"stl", "/entry/sample/stl/value", "degree"},
    /* Virtual Q/E */
    [TAS_COL_QH]   = {"qh",  "/entry/sample/qh/value",  "rlu"},
    [TAS_COL_QK]   = {"qk",  "/entry/sample/qk/value",  "rlu"},
    [TAS_COL_QL]   = {"ql",  "/entry/sample/ql/value",  "rlu"},
    [TAS_COL_EI]   = {"ei",  "/entry/instrument/monochromator/ei/value", "meV"},
    [TAS_COL_EF]   = {"ef",  "/entry/instrument/analyzer/ef/value",      "meV"},
    [TAS_COL_EN]   = {"en",  "/entry/virtual_motors/en/value",           "meV"},
    /* Monochromator */
    [TAS_COL_M1]       = {"m1",       "/entry/instrument/monochromator/m1/value",       "degree"},
    [TAS_COL_M2]       = {"m2",       "/entry/instrument/monochromator/m2/value",       "degree"},
    [TAS_COL_MONOVF]   = {"monovf",   "/entry/instrument/monochromator/monovf/value",   "degree"},
    [TAS_COL_MONOHF]   = {"monohf",   "/entry/instrument/monochromator/monohf/value",   "degree"},
    [TAS_COL_MONOTILT] = {"monotilt", "/entry/instrument/monochromator/monotilt/value", "degree"},
    [TAS_COL_MONOTRANS]= {"monotrans","/entry/instrument/monochromator/monotrans/value","mm"},
    /* Analyser */
    [TAS_COL_A1]       = {"a1",       "/entry/instrument/analyzer/a1/value",       "degree"},
    [TAS_COL_A2]       = {"a2",       "/entry/instrument/analyzer/a2/value",       "degree"},
    [TAS_COL_ANAVF]    = {"anavf",    "/entry/instrument/analyzer/anavf/value",    "degree"},
    [TAS_COL_ANAHF]    = {"anahf",    "/entry/instrument/analyzer/anahf/value",    "degree"},
    [TAS_COL_ANATILT]  = {"anatilt",  "/entry/instrument/analyzer/anatilt/value",  "degree"},
    [TAS_COL_ANATRANS] = {"anatrans", "/entry/instrument/analyzer/anatrans/value", "mm"},
    /* Source */
    [TAS_COL_VS_LEFT]  = {"vs_left",  "/entry/instrument/source/vs_left/value",  "mm"},
    [TAS_COL_VS_RIGHT] = {"vs_right", "/entry/instrument/source/vs_right/value", "mm"},
    /* Pre-sample slits */
    [TAS_COL_PS_LEFT]  = {"ps_left",  "/entry/instrument/slits/ps_left/value",  "mm"},
    [TAS_COL_PS_RIGHT] = {"ps_right", "/entry/instrument/slits/ps_right/value", "mm"},
    [TAS_COL_PS_TOP]   = {"ps_top",   "/entry/instrument/slits/ps_top/value",   "mm"},
    [TAS_COL_PS_BOTTOM]= {"ps_bottom","/entry/instrument/slits/ps_bottom/value","mm"},
    /* Post-analyser slits */
    [TAS_COL_PA_LEFT]  = {"pa_left",  "/entry/instrument/slits/pa_left/value",  "mm"},
    [TAS_COL_PA_RIGHT] = {"pa_right", "/entry/instrument/slits/pa_right/value", "mm"},
    [TAS_COL_PA_TOP]   = {"pa_top",   "/entry/instrument/slits/pa_top/value",   "mm"},
    [TAS_COL_PA_BOTTOM]= {"pa_bottom","/entry/instrument/slits/pa_bottom/value","mm"},
    /* Collimator motors */
    [TAS_COL_COL1] = {"col_1_motor","/entry/instrument/collimators/col_1_motor/value","arcmin"},
    [TAS_COL_COL2] = {"col_2_motor","/entry/instrument/collimators/col_2_motor/value","arcmin"},
    [TAS_COL_COL3] = {"col_3_motor","/entry/instrument/collimators/col_3_motor/value","arcmin"},
    [TAS_COL_COL4] = {"col_4_motor","/entry/instrument/collimators/col_4_motor/value","arcmin"},
    /* Collimator angular values */
    [TAS_COL_COLL_ALPHA1]={"coll_alpha1","/entry/instrument/collimators/coll_alpha1/value","arcmin"},
    [TAS_COL_COLL_ALPHA2]={"coll_alpha2","/entry/instrument/collimators/coll_alpha2/value","arcmin"},
    [TAS_COL_COLL_ALPHA3]={"coll_alpha3","/entry/instrument/collimators/coll_alpha3/value","arcmin"},
    [TAS_COL_COLL_ALPHA4]={"coll_alpha4","/entry/instrument/collimators/coll_alpha4/value","arcmin"},
    /* Sample environment */
    [TAS_COL_TEMP1]      = {"sample_temp1",    "/entry/sample/sample_env/sample_temp1/value",    "K"},
    [TAS_COL_TEMP2]      = {"sample_temp2",    "/entry/sample/sample_env/sample_temp2/value",    "K"},
    [TAS_COL_TEMP3]      = {"sample_temp3",    "/entry/sample/sample_env/sample_temp3/value",    "K"},
    [TAS_COL_TEMP4]      = {"sample_temp4",    "/entry/sample/sample_env/sample_temp4/value",    "K"},
    [TAS_COL_MFIELD]     = {"sample_mfield",   "/entry/sample/sample_env/sample_mfield/value",   "T"},
    [TAS_COL_EFIELD]     = {"sample_efield",   "/entry/sample/sample_env/sample_efield/value",   "V"},
    [TAS_COL_PRESSURE]   = {"sample_pressure", "/entry/sample/sample_env/sample_pressure/value", "GPa"},
    [TAS_COL_CRYO_HE]    = {"cryo_he_pressure","/entry/sample/sample_env/cryo_he_pressure/value","mBar"},
    [TAS_COL_CRYO_NEEDLE]= {"cryo_needlevalve","/entry/sample/sample_env/cryo_needlevalve/value","%"},
    /* Detectors */
    [TAS_COL_COUNTS]  = {"counts",  "/entry/instrument/det_group/detector/data", "counts"},
    [TAS_COL_MONITOR] = {"monitor", "/entry/instrument/monitor/data",             "counts"},
};

const char *tas_col_name(tas_col_id_t col)
{
    if (col < 0 || col >= TAS_N_COLS) return NULL;
    return COL_META[col].name;
}

int tas_col_from_name(const char *name)
{
    if (!name) return -1;
    for (int i = 0; i < TAS_N_COLS; i++)
        if (COL_META[i].name && strcmp(COL_META[i].name, name) == 0)
            return i;
    return -1;
}

/* =========================================================================
 * Internal HDF5 helpers
 * ========================================================================= */

/* String datatype for UTF-8 scalar strings */
static hid_t make_str_type(void)
{
    hid_t t = H5Tcopy(H5T_C_S1);
    H5Tset_size(t, H5T_VARIABLE);
    H5Tset_cset(t, H5T_CSET_UTF8);
    H5Tset_strpad(t, H5T_STR_NULLTERM);
    return t;
}

/* Write a fixed scalar UTF-8 string dataset */
static herr_t write_str_ds(hid_t loc, const char *name, const char *value)
{
    if (!value) value = "";
    hid_t str_t = make_str_type();
    hid_t sp    = H5Screate(H5S_SCALAR);
    hid_t ds    = H5Dcreate2(loc, name, str_t, sp,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t ret  = H5Dwrite(ds, str_t, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    H5Dclose(ds); H5Sclose(sp); H5Tclose(str_t);
    return ret;
}

/* Write a variable-length string array dataset (NULL-terminated char** src) */
static herr_t write_str_array_ds(hid_t loc, const char *name,
                                  const char * const *values, hsize_t count)
{
    hid_t str_t = make_str_type();
    hid_t sp    = H5Screate_simple(1, &count, NULL);
    hid_t ds    = H5Dcreate2(loc, name, str_t, sp,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t ret  = H5Dwrite(ds, str_t, H5S_ALL, H5S_ALL, H5P_DEFAULT, values);
    H5Dclose(ds); H5Sclose(sp); H5Tclose(str_t);
    return ret;
}

/* Write a scalar attribute string on an object */
static herr_t write_attr_str(hid_t obj, const char *name, const char *value)
{
    if (!value) value = "";
    hid_t str_t = make_str_type();
    hid_t sp    = H5Screate(H5S_SCALAR);
    hid_t at    = H5Acreate2(obj, name, str_t, sp,
                              H5P_DEFAULT, H5P_DEFAULT);
    herr_t ret  = H5Awrite(at, str_t, &value);
    H5Aclose(at); H5Sclose(sp); H5Tclose(str_t);
    return ret;
}

/* Write a scalar double dataset */
static herr_t write_scalar_double(hid_t loc, const char *name, double v)
{
    hid_t sp  = H5Screate(H5S_SCALAR);
    hid_t ds  = H5Dcreate2(loc, name, H5T_NATIVE_DOUBLE, sp,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t r  = H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
    H5Dclose(ds); H5Sclose(sp);
    return r;
}

/* Write a double attribute on an object */
static void write_attr_double(hid_t obj, const char *name, double v)
{
    hid_t sp = H5Screate(H5S_SCALAR);
    hid_t at = H5Acreate2(obj, name, H5T_NATIVE_DOUBLE, sp,
                           H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(at, H5T_NATIVE_DOUBLE, &v);
    H5Aclose(at); H5Sclose(sp);
}

/* Write a scalar int attribute */
static void write_attr_int(hid_t obj, const char *name, int v)
{
    hid_t sp = H5Screate(H5S_SCALAR);
    hid_t at = H5Acreate2(obj, name, H5T_NATIVE_INT, sp,
                           H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(at, H5T_NATIVE_INT, &v);
    H5Aclose(at); H5Sclose(sp);
}

/* Write a 1-D double array dataset */
static herr_t write_double_array(hid_t loc, const char *name,
                                  const double *data, hsize_t len)
{
    hid_t sp = H5Screate_simple(1, &len, NULL);
    hid_t ds = H5Dcreate2(loc, name, H5T_NATIVE_DOUBLE, sp,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t r = H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds); H5Sclose(sp);
    return r;
}

/* Pre-allocate a 1-D double dataset of length n filled with zeros */
static hid_t create_preallocated_1d(hid_t loc, const char *name, hsize_t n)
{
    hid_t sp = H5Screate_simple(1, &n, NULL);
    hid_t ds = H5Dcreate2(loc, name, H5T_NATIVE_DOUBLE, sp,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /* HDF5 zero-fills by default for fixed-size datasets */
    H5Sclose(sp);
    return ds;   /* caller must H5Dclose */
}

/* Pre-allocate a 3-D double dataset [n][NY][NX] filled with zeros */
static hid_t create_preallocated_3d(hid_t loc, const char *name,
                                     hsize_t n, hsize_t ny, hsize_t nx)
{
    hsize_t dims[3] = {n, ny, nx};
    hid_t sp = H5Screate_simple(3, dims, NULL);
    hid_t ds = H5Dcreate2(loc, name, H5T_NATIVE_DOUBLE, sp,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(sp);
    return ds;
}

/* Create an NXpositioner group with a pre-allocated value[n] dataset */
static hid_t make_positioner(hid_t parent, const char *name, hsize_t n)
{
    hid_t g  = H5Gcreate2(parent, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(g, "NX_class", "NXpositioner");
    hid_t ds = create_preallocated_1d(g, "value", n);

    const char *unit = "degree";
    for (int i = 0; i < TAS_N_COLS; i++)
        if (COL_META[i].name && strcmp(COL_META[i].name, name) == 0) {
            unit = COL_META[i].unit; break;
        }
    write_attr_str(ds, "units", unit);
    H5Dclose(ds);
    return g;
}

/* Create an NXsensor group with a pre-allocated value[n] dataset */
static hid_t make_sensor(hid_t parent, const char *name,
                          const char *unit, hsize_t n)
{
    hid_t g  = H5Gcreate2(parent, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(g, "NX_class", "NXsensor");
    hid_t ds = create_preallocated_1d(g, "value", n);
    write_attr_str(ds, "units", unit);
    H5Dclose(ds);
    return g;
}

/* =========================================================================
 * ISO-8601 timestamp
 * ========================================================================= */

static void now_iso(char *buf, size_t bufsz)
{
    time_t t = time(NULL);
    struct tm tm_utc;
#ifdef _WIN32
    gmtime_s(&tm_utc, &t);
#else
    gmtime_r(&t, &tm_utc);
#endif
    strftime(buf, bufsz, "%Y-%m-%dT%H:%M:%SZ", &tm_utc);
}

/* =========================================================================
 * Metadata helpers
 * ========================================================================= */

#define META_STR(meta, field, def) \
    ((meta) && (meta)->field && (meta)->field[0] ? (meta)->field : (def))

#define META_DBL(meta, field, def) \
    ((meta) && (meta)->field != 0.0 ? (meta)->field : (def))

/* =========================================================================
 * Default values
 * ========================================================================= */

static const double DEFAULT_UNIT_CELL[6]  = {3.82, 3.82, 11.68, 90.0, 90.0, 90.0};
static const double DEFAULT_UB[3][3] = {
    {0.1254, 0.0021, 0.0},
    {0.0018, 0.1189, 0.0},
    {0.0,    0.0,    0.0765}
};
static const double DEFAULT_V1[3] = {1.0, 0.0, 0.0};
static const double DEFAULT_V2[3] = {0.0, 1.0, 0.0};

/* =========================================================================
 * Skeleton builder — mirrors Python _build_skeleton()
 * ========================================================================= */

static int build_skeleton(hid_t f, const tas_metadata_t *meta, hsize_t n)
{
    char ts_now[32];
    now_iso(ts_now, sizeof(ts_now));

    /* ── /entry ──────────────────────────────────────────────────────── */
    hid_t entry = H5Gcreate2(f, "entry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(entry, "NX_class", "NXentry");
    write_attr_str(entry, "title", META_STR(meta, title, "TAS scan"));
    write_attr_str(entry, "default", "data");

    write_str_ds(entry, "definition", "NXtas");
    write_str_ds(entry, "start_time",
                 META_STR(meta, start_time, ts_now));
    write_str_ds(entry, "end_time",
                 META_STR(meta, end_time, ""));

    /* ── instrument ───────────────────────────────────────────────────── */
    hid_t inst = H5Gcreate2(entry, "instrument",
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(inst, "NX_class", "NXinstrument");

    /* source */
    hid_t src = H5Gcreate2(inst, "source",
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(src, "NX_class", "NXsource");
    H5Gclose(make_positioner(src, "vs_left",  n));
    H5Gclose(make_positioner(src, "vs_right", n));
    H5Gclose(src);

    /* monochromator */
    hid_t mono = H5Gcreate2(inst, "monochromator",
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(mono, "NX_class", "NXmonochromator");
    {
        double d = meta && meta->distance_mono_sample != 0.0
                 ? meta->distance_mono_sample : 2.0;
        hid_t ds = create_preallocated_1d(mono, "distance", 1);
        double tmp = d;
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
        write_attr_str(ds, "units", "m");
        H5Dclose(ds);
        /* scalar scalar overwrite via write_scalar_double would be simpler,
           but create_preallocated_1d keeps shape consistent; use a true scalar: */
        H5Ldelete(mono, "distance", H5P_DEFAULT);
        write_scalar_double(mono, "distance", d);
        {
            hid_t dst = H5Dopen2(mono, "distance", H5P_DEFAULT);
            write_attr_str(dst, "units", "m");
            H5Dclose(dst);
        }
    }
    write_str_ds(mono, "mono_crystal", META_STR(meta, mono_crystal, "PG"));
    const char *mono_pos[] = {"m1","m2","monovf","monohf","monotilt","monotrans","ei"};
    for (int i = 0; i < 7; i++) H5Gclose(make_positioner(mono, mono_pos[i], n));
    H5Gclose(mono);

    /* analyzer */
    hid_t ana = H5Gcreate2(inst, "analyzer",
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(ana, "NX_class", "NXmonochromator");
    {
        double d = meta && meta->distance_sample_ana != 0.0
                 ? meta->distance_sample_ana : 1.5;
        write_scalar_double(ana, "distance", d);
        hid_t dst = H5Dopen2(ana, "distance", H5P_DEFAULT);
        write_attr_str(dst, "units", "m");
        H5Dclose(dst);
    }
    write_str_ds(ana, "ana_crystal", META_STR(meta, ana_crystal, "PG"));
    const char *ana_pos[] = {"a1","a2","anavf","anahf","anatilt","anatrans","ef"};
    for (int i = 0; i < 7; i++) H5Gclose(make_positioner(ana, ana_pos[i], n));
    H5Gclose(ana);

    /* det_group → detector */
    hid_t dg = H5Gcreate2(inst, "det_group",
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(dg, "NX_class", "NXdetector_group");
    {
        hid_t sd = H5Gcreate2(dg, "detector",
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write_attr_str(sd, "NX_class", "NXdetector");
        {
            hid_t ds = create_preallocated_1d(sd, "data", n);
            write_attr_str(ds, "units", "counts");
            H5Dclose(ds);
        }
        {
            double d = meta && meta->distance_ana_det != 0.0
                     ? meta->distance_ana_det : 0.5;
            write_scalar_double(sd, "distance", d);
            hid_t ds = H5Dopen2(sd, "distance", H5P_DEFAULT);
            write_attr_str(ds, "units", "m");
            H5Dclose(ds);
        }
        H5Gclose(sd);
    }
    /* det_group → psd */
    {
        hid_t psd = H5Gcreate2(dg, "psd",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write_attr_str(psd, "NX_class", "NXdetector");
        {
            hid_t ds = create_preallocated_3d(psd, "data",
                                               n, TAS_PSD_NY, TAS_PSD_NX);
            write_attr_str(ds, "units", "counts");
            H5Dclose(ds);
        }
        /* pixel offsets */
        {
            double *px = (double *)malloc(TAS_PSD_NX * sizeof(double));
            for (int i = 0; i < TAS_PSD_NX; i++) px[i] = (double)i;
            write_double_array(psd, "x_pixel_offset", px, TAS_PSD_NX);
            hid_t ds = H5Dopen2(psd, "x_pixel_offset", H5P_DEFAULT);
            write_attr_str(ds, "units", "pixel"); H5Dclose(ds);
            write_double_array(psd, "y_pixel_offset", px, TAS_PSD_NY);
            ds = H5Dopen2(psd, "y_pixel_offset", H5P_DEFAULT);
            write_attr_str(ds, "units", "pixel"); H5Dclose(ds);
            free(px);
        }
        H5Gclose(psd);
    }
    H5Gclose(dg);

    /* monitor */
    {
        hid_t mon = H5Gcreate2(inst, "monitor",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write_attr_str(mon, "NX_class", "NXmonitor");
        hid_t ds = create_preallocated_1d(mon, "data", n);
        write_attr_str(ds, "units", "counts");
        H5Dclose(ds);
        H5Gclose(mon);
    }

    /* slits */
    {
        hid_t slits = H5Gcreate2(inst, "slits",
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write_attr_str(slits, "NX_class", "NXcollection");
        const char *snames[] = {
            "ps_left","ps_right","ps_top","ps_bottom",
            "pa_left","pa_right","pa_top","pa_bottom"
        };
        for (int i = 0; i < 8; i++) H5Gclose(make_positioner(slits, snames[i], n));
        H5Gclose(slits);
    }

    /* collimators */
    {
        hid_t coll = H5Gcreate2(inst, "collimators",
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write_attr_str(coll, "NX_class", "NXcollection");
        const char *cnames[] = {
            "col_1_motor","col_2_motor","col_3_motor","col_4_motor",
            "coll_alpha1","coll_alpha2","coll_alpha3","coll_alpha4"
        };
        for (int i = 0; i < 8; i++) H5Gclose(make_positioner(coll, cnames[i], n));
        H5Gclose(coll);
    }

    /* instrument-level geometry */
    {
        struct { const char *name; double val; } dists[] = {
            {"distance_vs_mono",     meta && meta->distance_vs_mono     ? meta->distance_vs_mono     : 0.0},
            {"distance_mono_sample", meta && meta->distance_mono_sample ? meta->distance_mono_sample : 2.0},
            {"distance_sample_ana",  meta && meta->distance_sample_ana  ? meta->distance_sample_ana  : 1.5},
            {"distance_ana_det",     meta && meta->distance_ana_det     ? meta->distance_ana_det     : 0.5},
        };
        for (int i = 0; i < 4; i++) {
            write_scalar_double(inst, dists[i].name, dists[i].val);
            hid_t ds = H5Dopen2(inst, dists[i].name, H5P_DEFAULT);
            write_attr_str(ds, "units", "m");
            H5Dclose(ds);
        }
        write_str_ds(inst, "sense", META_STR(meta, sense, "+-+"));
    }
    H5Gclose(inst);

    /* ── sample ───────────────────────────────────────────────────────── */
    hid_t sample = H5Gcreate2(entry, "sample",
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(sample, "NX_class", "NXsample");
    {
        double d = meta && meta->distance_mono_sample ? meta->distance_mono_sample : 1.5;
        write_scalar_double(sample, "distance", d);
        hid_t ds = H5Dopen2(sample, "distance", H5P_DEFAULT);
        write_attr_str(ds, "units", "m"); H5Dclose(ds);
    }

    /* unit_cell */
    {
        const double *uc = (meta && meta->unit_cell[0] != 0.0)
                         ? meta->unit_cell : DEFAULT_UNIT_CELL;
        hsize_t six = 6;
        hid_t sp = H5Screate_simple(1, &six, NULL);
        hid_t ds = H5Dcreate2(sample, "unit_cell", H5T_NATIVE_DOUBLE, sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, uc);
        H5Dclose(ds); H5Sclose(sp);
    }

    /* ub_matrix */
    {
        const double (*ub)[3] = (meta && meta->ub_matrix[0][0] != 0.0)
                               ? (const double (*)[3])meta->ub_matrix
                               : (const double (*)[3])DEFAULT_UB;
        hsize_t dims[2] = {3, 3};
        hid_t sp = H5Screate_simple(2, dims, NULL);
        hid_t ds = H5Dcreate2(sample, "ub_matrix", H5T_NATIVE_DOUBLE, sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ub);
        write_attr_str(ds, "description",   "Orientation matrix transforming HKL to Q");
        write_attr_str(ds, "interpretation","matrix");
        H5Dclose(ds); H5Sclose(sp);
    }

    /* sample_mosaic, sample_v1, sample_v2 */
    {
        double mosaic = (meta) ? meta->sample_mosaic : 0.0;
        write_scalar_double(sample, "sample_mosaic", mosaic);
        hid_t ds = H5Dopen2(sample, "sample_mosaic", H5P_DEFAULT);
        write_attr_str(ds, "units", "minutes of arc"); H5Dclose(ds);

        const double *v1 = (meta && meta->sample_v1[0] != 0.0) ? meta->sample_v1 : DEFAULT_V1;
        const double *v2 = (meta && meta->sample_v2[1] != 0.0) ? meta->sample_v2 : DEFAULT_V2;
        hsize_t three = 3;
        hid_t sp = H5Screate_simple(1, &three, NULL);
        ds = H5Dcreate2(sample, "sample_v1", H5T_NATIVE_DOUBLE, sp,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v1);
        H5Dclose(ds);
        ds = H5Dcreate2(sample, "sample_v2", H5T_NATIVE_DOUBLE, sp,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v2);
        H5Dclose(ds);
        H5Sclose(sp);
    }

    /* sample positioners */
    const char *samp_pos[] = {"s1","s2","sgu","sgl","stu","stl","qh","qk","ql"};
    for (int i = 0; i < 9; i++) H5Gclose(make_positioner(sample, samp_pos[i], n));

    /* sample environment */
    {
        hid_t se = H5Gcreate2(sample, "sample_env",
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write_attr_str(se, "NX_class", "NXenvironment");
        write_str_ds(se, "name",        "Cryostat CF16");
        write_str_ds(se, "type",        "Cryostat");
        write_str_ds(se, "description", "Closed-cycle cryostat");
        struct { const char *name; const char *unit; } sensors[] = {
            {"sample_temp1","K"},   {"sample_temp2","K"},
            {"sample_temp3","K"},   {"sample_temp4","K"},
            {"sample_mfield","T"},  {"sample_efield","V"},
            {"sample_pressure","GPa"},
            {"cryo_he_pressure","mBar"},
            {"cryo_needlevalve","%"},
        };
        for (int i = 0; i < 9; i++)
            H5Gclose(make_sensor(se, sensors[i].name, sensors[i].unit, n));
        H5Gclose(se);
    }
    H5Gclose(sample);

    /* ── virtual_motors ───────────────────────────────────────────────── */
    hid_t vm = H5Gcreate2(entry, "virtual_motors",
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attr_str(vm, "NX_class", "NXcollection");

    /* hard-links into sample */
    H5Lcreate_hard(f, "/entry/sample/qh", vm, "qh", H5P_DEFAULT, H5P_DEFAULT);
    H5Lcreate_hard(f, "/entry/sample/qk", vm, "qk", H5P_DEFAULT, H5P_DEFAULT);
    H5Lcreate_hard(f, "/entry/sample/ql", vm, "ql", H5P_DEFAULT, H5P_DEFAULT);
    /* hard-links into instrument */
    H5Lcreate_hard(f, "/entry/instrument/monochromator/ei", vm, "ei", H5P_DEFAULT, H5P_DEFAULT);
    H5Lcreate_hard(f, "/entry/instrument/analyzer/ef",      vm, "ef", H5P_DEFAULT, H5P_DEFAULT);

    /* en: its own group */
    {
        hid_t en_g = H5Gcreate2(vm, "en",
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write_attr_str(en_g, "NX_class", "NXpositioner");
        hid_t ds = create_preallocated_1d(en_g, "value", n);
        write_attr_str(ds, "units", "meV");
        H5Dclose(ds);
        H5Gclose(en_g);
    }
    H5Gclose(vm);

    /* ── metadata ─────────────────────────────────────────────────────── */
    {
        hid_t md = H5Gcreate2(entry, "metadata",
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write_str_ds(md, "facility",         META_STR(meta, facility,         "Facility_Name"));
        write_str_ds(md, "source",           META_STR(meta, source,           ""));
        write_str_ds(md, "instrument_name",  META_STR(meta, instrument_name,  "TAS"));
        write_str_ds(md, "software version", META_STR(meta, software_version, "0.9.1"));
        write_str_ds(md, "TAS_NeXus_Version",META_STR(meta, tas_nexus_version,"0.9.1"));
        write_str_ds(md, "experiment_id",    META_STR(meta, experiment_id,    ""));
        write_str_ds(md, "proposal_no",      META_STR(meta, proposal_no,      ""));
        write_str_ds(md, "local_contact",    META_STR(meta, local_contact,    ""));
        write_str_ds(md, "sample_name",      META_STR(meta, sample_name,      ""));
        write_str_ds(md, "sample_type",      META_STR(meta, sample_type,      "crystal"));
        write_str_ds(md, "scan_no",          META_STR(meta, scan_no,          "0"));
        write_str_ds(md, "title",            META_STR(meta, title,            ""));
        write_str_ds(md, "command",          META_STR(meta, command,          ""));
        write_str_ds(md, "filename",         META_STR(meta, filename,         ""));
        write_str_ds(md, "scanning_axis",    META_STR(meta, scanning_axis,    "s2"));

        /* users — variable-length string array */
        if (meta && meta->users) {
            hsize_t nu = 0;
            while (meta->users[nu]) nu++;
            if (nu > 0)
                write_str_array_ds(md, "user", meta->users, nu);
            else
                write_str_ds(md, "user", "Unknown");
        } else {
            write_str_ds(md, "user", "Unknown");
        }
        H5Gclose(md);
    }

    /* ── NXdata (default plot view) ───────────────────────────────────── */
    {
        hid_t nxdata = H5Gcreate2(entry, "data",
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write_attr_str(nxdata, "NX_class",       "NXdata");
        write_attr_str(nxdata, "signal",          "data");
        write_attr_str(nxdata, "axes",            "s2");
        write_attr_str(nxdata, "interpretation",  "spectrum");
        write_attr_int(nxdata, "s2_indices",  0);
        write_attr_int(nxdata, "en_indices",  0);
        write_attr_int(nxdata, "qh_indices",  0);

        /* hard-links into the real datasets */
        H5Lcreate_hard(f, "/entry/instrument/det_group/detector/data",
                        nxdata, "data",    H5P_DEFAULT, H5P_DEFAULT);
        H5Lcreate_hard(f, "/entry/instrument/monitor/data",
                        nxdata, "monitor", H5P_DEFAULT, H5P_DEFAULT);
        H5Lcreate_hard(f, "/entry/sample/s2/value",
                        nxdata, "s2",      H5P_DEFAULT, H5P_DEFAULT);
        H5Lcreate_hard(f, "/entry/virtual_motors/qh/value",
                        nxdata, "qh",      H5P_DEFAULT, H5P_DEFAULT);
        H5Lcreate_hard(f, "/entry/virtual_motors/qk/value",
                        nxdata, "qk",      H5P_DEFAULT, H5P_DEFAULT);
        H5Lcreate_hard(f, "/entry/virtual_motors/ql/value",
                        nxdata, "ql",      H5P_DEFAULT, H5P_DEFAULT);
        H5Lcreate_hard(f, "/entry/virtual_motors/en/value",
                        nxdata, "en",      H5P_DEFAULT, H5P_DEFAULT);
        H5Lcreate_hard(f, "/entry/instrument/monochromator/ei/value",
                        nxdata, "ei",      H5P_DEFAULT, H5P_DEFAULT);
        H5Lcreate_hard(f, "/entry/instrument/analyzer/ef/value",
                        nxdata, "ef",      H5P_DEFAULT, H5P_DEFAULT);

        /* long_name attributes */
        struct { const char *ds; const char *longname; } lnames[] = {
            {"qh","QH Reciprocal Coordinate"}, {"qk","QK Reciprocal Coordinate"},
            {"ql","QL Reciprocal Coordinate"}, {"en","Energy Transfer"},
            {"ei","Incident Energy"},           {"ef","Final Energy"},
        };
        for (int i = 0; i < 6; i++) {
            hid_t ds = H5Dopen2(nxdata, lnames[i].ds, H5P_DEFAULT);
            write_attr_str(ds, "long_name", lnames[i].longname);
            H5Dclose(ds);
        }
        H5Gclose(nxdata);
    }

    /* root default attribute */
    write_attr_str(f, "default", "entry");

    H5Gclose(entry);
    return TAS_OK;
}

/* =========================================================================
 * Write scalar data (mirrors Python _write_scalars)
 * ========================================================================= */

/*
 * Write the array `data[n]` into dataset at `path`, starting at offset `start`.
 * If data is NULL, writes zeros (dataset was pre-zeroed, so this is a no-op).
 */
static int write_scalar_slice(hid_t f, const char *path,
                               const double *data, hsize_t start, hsize_t n)
{
    if (!data) return TAS_OK;   /* dataset already zeroed */

    hid_t ds  = H5Dopen2(f, path, H5P_DEFAULT);
    if (ds < 0) { set_error("Dataset not found: %s", path); return TAS_ERR_HDF5; }

    hid_t fsp = H5Dget_space(ds);
    hsize_t off = start, cnt = n;
    H5Sselect_hyperslab(fsp, H5S_SELECT_SET, &off, NULL, &cnt, NULL);

    hid_t msp = H5Screate_simple(1, &cnt, NULL);
    herr_t r  = H5Dwrite(ds, H5T_NATIVE_DOUBLE, msp, fsp, H5P_DEFAULT, data);
    H5Sclose(msp); H5Sclose(fsp); H5Dclose(ds);
    return (r < 0) ? TAS_ERR_HDF5 : TAS_OK;
}

/* Write all scalar columns from a lookup array[TAS_N_COLS] */
static int write_all_scalars_from_array(hid_t f,
                                         const double values[TAS_N_COLS],
                                         hsize_t start)
{
    for (int c = 0; c < TAS_N_COLS; c++) {
        double v = values[c];
        int rc = write_scalar_slice(f, COL_META[c].hdf5_path, &v, start, 1);
        if (rc != TAS_OK) return rc;
    }
    return TAS_OK;
}

/* Write a batch of scalar columns from tas_scan_data_t */
static int write_batch_scalars(hid_t f, const tas_scan_data_t *scan)
{
    for (size_t ci = 0; ci < scan->n_columns; ci++) {
        const tas_scalar_col_t *col = &scan->columns[ci];
        if (col->column < 0 || col->column >= TAS_N_COLS) {
            set_error("Invalid column ID %d", col->column);
            return TAS_ERR_ARG;
        }
        if (!col->data) {
            set_error("NULL data pointer for column %d (%s)",
                      col->column, COL_META[col->column].name);
            return TAS_ERR_NULL_DATA;
        }
        int rc = write_scalar_slice(f, COL_META[col->column].hdf5_path,
                                    col->data, 0, (hsize_t)scan->n_points);
        if (rc != TAS_OK) return rc;
    }
    return TAS_OK;
}

/* =========================================================================
 * Write PSD data (mirrors Python _write_psd)
 * ========================================================================= */

static int write_psd_slice(hid_t f, const double *psd,
                            hsize_t start, hsize_t n_frames)
{
    if (!psd) return TAS_OK;

    hid_t ds  = H5Dopen2(f, "/entry/instrument/det_group/psd/data", H5P_DEFAULT);
    if (ds < 0) { set_error("PSD dataset not found"); return TAS_ERR_HDF5; }

    hid_t fsp = H5Dget_space(ds);
    hsize_t offset[3] = {start, 0, 0};
    hsize_t count[3]  = {n_frames, TAS_PSD_NY, TAS_PSD_NX};
    H5Sselect_hyperslab(fsp, H5S_SELECT_SET, offset, NULL, count, NULL);

    hid_t msp = H5Screate_simple(3, count, NULL);
    herr_t r  = H5Dwrite(ds, H5T_NATIVE_DOUBLE, msp, fsp, H5P_DEFAULT, psd);
    H5Sclose(msp); H5Sclose(fsp); H5Dclose(ds);
    return (r < 0) ? TAS_ERR_HDF5 : TAS_OK;
}

/* =========================================================================
 * Stamp end_time into /entry/end_time
 * ========================================================================= */

static void stamp_end_time(hid_t f)
{
    char ts[32];
    now_iso(ts, sizeof(ts));
    hid_t ds    = H5Dopen2(f, "/entry/end_time", H5P_DEFAULT);
    if (ds < 0) return;
    hid_t str_t = make_str_type();
    const char *p = ts;
    H5Dwrite(ds, str_t, H5S_ALL, H5S_ALL, H5P_DEFAULT, &p);
    H5Tclose(str_t);
    H5Dclose(ds);
}

/* =========================================================================
 * HDF5 file creation helper
 * ========================================================================= */

static hid_t open_new_file(const char *filepath)
{
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    hid_t fcpl = H5Pcreate(H5P_FILE_CREATE);
    H5Pset_link_creation_order(fcpl,
        H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
    H5Pset_attr_creation_order(fcpl,
        H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
    hid_t fid = H5Fcreate(filepath, H5F_ACC_TRUNC, fcpl, fapl);
    H5Pclose(fcpl); H5Pclose(fapl);
    return fid;
}

/* =========================================================================
 * Public API — Batch
 * ========================================================================= */

int tas_save_hdf(const char            *filepath,
                 const tas_scan_data_t *scan,
                 const tas_metadata_t  *meta,
                 int                    flags)
{
    if (!filepath || !scan) {
        set_error("tas_save_hdf: filepath and scan must not be NULL");
        return TAS_ERR_ARG;
    }
    if (scan->n_points == 0) {
        set_error("tas_save_hdf: scan.n_points must be >= 1");
        return TAS_ERR_ARG;
    }

    /* Check overwrite */
    if (!(flags & TAS_OVERWRITE)) {
        FILE *fp = fopen(filepath, "r");
        if (fp) { fclose(fp); set_error("File exists: %s", filepath); return TAS_ERR_EXISTS; }
    }

    /* Suppress HDF5 default error output during file ops (we handle errors) */
    H5E_auto2_t old_func; void *old_data;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_data);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    hid_t fid = open_new_file(filepath);
    if (fid < 0) {
        H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
        set_error("Failed to create file: %s", filepath);
        return TAS_ERR_HDF5;
    }

    int rc = build_skeleton(fid, meta, (hsize_t)scan->n_points);
    if (rc == TAS_OK) rc = write_batch_scalars(fid, scan);
    if (rc == TAS_OK && scan->psd_data)
        rc = write_psd_slice(fid, scan->psd_data, 0, (hsize_t)scan->n_points);

    /* Stamp end_time if not provided */
    if (rc == TAS_OK && !(meta && meta->end_time && meta->end_time[0]))
        stamp_end_time(fid);

    H5Fclose(fid);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_data);

    if (rc == TAS_OK)
        fprintf(stdout, "[tas_save_hdf] Wrote %zu points -> %s\n",
                scan->n_points, filepath);
    return rc;
}

/* =========================================================================
 * Public API — Point-by-point handle
 * ========================================================================= */

struct tas_pbp_handle {
    hid_t   fid;          /* open HDF5 file id          */
    size_t  n_points;     /* total pre-allocated points */
    size_t  next_index;   /* next slot to write (0-based) */
    char    filepath[4096];
};

tas_pbp_handle_t *tas_pbp_open(const char           *filepath,
                                size_t                n_points,
                                const tas_metadata_t *meta)
{
    if (!filepath) { set_error("tas_pbp_open: filepath is NULL"); return NULL; }
    if (n_points == 0) { set_error("tas_pbp_open: n_points must be >= 1"); return NULL; }

    H5E_auto2_t old_func; void *old_data;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_data);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    hid_t fid = open_new_file(filepath);
    if (fid < 0) {
        H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
        set_error("Failed to create file: %s", filepath);
        return NULL;
    }

    if (build_skeleton(fid, meta, (hsize_t)n_points) != TAS_OK) {
        H5Fclose(fid);
        H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
        return NULL;
    }

    H5Eset_auto2(H5E_DEFAULT, old_func, old_data);

    tas_pbp_handle_t *h = (tas_pbp_handle_t *)calloc(1, sizeof(*h));
    if (!h) { H5Fclose(fid); set_error("Out of memory"); return NULL; }

    h->fid        = fid;
    h->n_points   = n_points;
    h->next_index = 0;
    snprintf(h->filepath, sizeof(h->filepath), "%s", filepath);
    return h;
}

int tas_pbp_write_point(tas_pbp_handle_t  *handle,
                        const tas_point_t *point)
{
    if (!handle || !point) {
        set_error("tas_pbp_write_point: NULL argument");
        return TAS_ERR_ARG;
    }
    if (handle->next_index >= handle->n_points) {
        set_error("tas_pbp_write_point: exceeded pre-allocated n_points=%zu",
                  handle->n_points);
        return TAS_ERR_ARG;
    }

    hsize_t idx = (hsize_t)handle->next_index;

    int rc = write_all_scalars_from_array(handle->fid, point->values, idx);
    if (rc == TAS_OK && point->psd_frame)
        rc = write_psd_slice(handle->fid, point->psd_frame, idx, 1);

    if (rc == TAS_OK) {
        handle->next_index++;
        int current = (int)handle->next_index;
        fprintf(stdout, "[tas_pbp_write_point] Point %d/%zu -> %s\n",
                current, handle->n_points, handle->filepath);
        return current;
    }
    return rc;
}

int tas_pbp_close(tas_pbp_handle_t *handle)
{
    if (!handle) { set_error("tas_pbp_close: NULL handle"); return TAS_ERR_ARG; }
    stamp_end_time(handle->fid);
    H5Fclose(handle->fid);
    free(handle);
    return TAS_OK;
}
