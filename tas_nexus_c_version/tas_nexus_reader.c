/* Request POSIX.1-2008 extensions before any system header */
#define _POSIX_C_SOURCE 200809L

/**
 * @file tas_nexus_reader.c
 * @brief TAS NeXus HDF5 reader — C implementation
 *
 * Mirrors tas_nexus_reader.py:
 *   - load_from_hdf()   → tas_read_hdf()
 *   - export_to_text()  → tas_export_text()
 *   - self-test summary → tas_print_summary()
 */

#include "tas_nexus_reader.h"
#include "tas_nexus_writer.h"   /* COL_META via extern, TAS_N_COLS, TAS_PSD_* */

#include <hdf5.h>

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* =========================================================================
 * Error buffer  (same pattern as the writer)
 * ========================================================================= */

#ifndef TAS_ERROR_BUF_SIZE
#  define TAS_ERROR_BUF_SIZE 512
#endif

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#  define TAS_TL _Thread_local
#elif defined(__GNUC__)
#  define TAS_TL __thread
#elif defined(_MSC_VER)
#  define TAS_TL __declspec(thread)
#else
#  define TAS_TL
#endif

static TAS_TL char s_errbuf[TAS_ERROR_BUF_SIZE];

static void set_error(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(s_errbuf, sizeof(s_errbuf), fmt, ap);
    va_end(ap);
}

/* The writer already defines tas_last_error(); we reuse it via the shared
   translation unit.  If the reader is compiled standalone without the writer,
   uncomment the following definition: */
/*
const char *tas_last_error(void) { return s_errbuf; }
*/

/* =========================================================================
 * HDF5 path tables  (mirrors Python _SCALAR_PATHS, _META_* dicts)
 *
 * These are intentionally defined here (not in the header) so that changing
 * the on-disk layout only requires updating this file.
 * ========================================================================= */

/* The writer's COL_META table gives us the hdf5_path for each column ID.
   We access it via an extern declaration of the internal array. */

/* col_meta_t is defined in tas_nexus_writer.c; COL_META is a file-scope array
   there.  We use a matching struct here and declare it extern so the linker
   resolves it from the writer object. */
typedef struct { const char *name; const char *hdf5_path; const char *unit; } col_meta_t;
extern const col_meta_t COL_META[];   /* defined in tas_nexus_writer.c */

/* ── Metadata paths ──────────────────────────────────────────────────────── */

typedef struct { const char *key; const char *path; } str_path_t;
typedef struct { const char *key; const char *path; } flt_path_t;

static const str_path_t META_STR[] = {
    {"start_time",      "/entry/start_time"},
    {"end_time",        "/entry/end_time"},
    {"title",           "/entry/metadata/title"},
    {"command",         "/entry/metadata/command"},
    {"filename",        "/entry/metadata/filename"},
    {"facility",        "/entry/metadata/facility"},
    {"source",          "/entry/metadata/source"},
    {"instrument_name", "/entry/metadata/instrument_name"},
    {"local_contact",   "/entry/metadata/local_contact"},
    {"experiment_id",   "/entry/metadata/experiment_id"},
    {"proposal_no",     "/entry/metadata/proposal_no"},
    {"sample_name",     "/entry/metadata/sample_name"},
    {"sample_type",     "/entry/metadata/sample_type"},
    {"scanning_axis",   "/entry/metadata/scanning_axis"},
    {"scan_no",         "/entry/metadata/scan_no"},
    {"mono_crystal",    "/entry/instrument/monochromator/mono_crystal"},
    {"ana_crystal",     "/entry/instrument/analyzer/ana_crystal"},
    {"sense",           "/entry/instrument/sense"},
    {NULL, NULL}
};

static const flt_path_t META_FLT[] = {
    {"sample_mosaic",        "/entry/sample/sample_mosaic"},
    {"distance_vs_mono",     "/entry/instrument/distance_vs_mono"},
    {"distance_mono_sample", "/entry/instrument/distance_mono_sample"},
    {"distance_sample_ana",  "/entry/instrument/distance_sample_ana"},
    {"distance_ana_det",     "/entry/instrument/distance_ana_det"},
    {NULL, NULL}
};

/* =========================================================================
 * Internal HDF5 helpers
 * ========================================================================= */

/* Returns 1 if a dataset or group at 'path' exists in file f */
static int h5_exists(hid_t f, const char *path)
{
    return (H5Lexists(f, path, H5P_DEFAULT) > 0) ? 1 : 0;
}

/**
 * Read a scalar string dataset at 'path' into a newly malloc'd C string.
 * Caller must free().  Returns strdup(default_val) if absent.
 */
static char *h5_read_str(hid_t f, const char *path, const char *default_val)
{
    if (!h5_exists(f, path))
        return strdup(default_val ? default_val : "");

    hid_t ds  = H5Dopen2(f, path, H5P_DEFAULT);
    hid_t typ = H5Dget_type(ds);
    char *result = NULL;

    if (H5Tis_variable_str(typ)) {
        /* Variable-length UTF-8 string */
        char *buf = NULL;
        hid_t mem = H5Tcopy(H5T_C_S1);
        H5Tset_size(mem, H5T_VARIABLE);
        H5Tset_cset(mem, H5T_CSET_UTF8);
        H5Tset_strpad(mem, H5T_STR_NULLTERM);
        H5Dread(ds, mem, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf);
        result = strdup(buf ? buf : "");
        H5free_memory(buf);
        H5Tclose(mem);
    } else {
        /* Fixed-length string */
        size_t sz = H5Tget_size(typ);
        char *buf = (char *)calloc(sz + 1, 1);
        H5Dread(ds, typ, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        result = strdup(buf);
        free(buf);
    }

    H5Tclose(typ);
    H5Dclose(ds);
    return result ? result : strdup("");
}

/**
 * Read a scalar double dataset at 'path'.
 * Returns default_val if absent or on error.
 */
static double h5_read_double(hid_t f, const char *path, double default_val)
{
    if (!h5_exists(f, path)) return default_val;
    hid_t ds = H5Dopen2(f, path, H5P_DEFAULT);
    double v = default_val;
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
    H5Dclose(ds);
    return v;
}

/**
 * Read a fixed-size double array dataset at 'path' into dst[count].
 * Returns 0 on success, -1 if absent.
 */
static int h5_read_double_array(hid_t f, const char *path,
                                 double *dst, size_t count)
{
    if (!h5_exists(f, path)) return -1;
    hid_t ds = H5Dopen2(f, path, H5P_DEFAULT);

    /* Flatten to 1D regardless of dataset shape */
    hsize_t total = (hsize_t)count;
    hid_t msp = H5Screate_simple(1, &total, NULL);
    H5Dread(ds, H5T_NATIVE_DOUBLE, msp, H5S_ALL, H5P_DEFAULT, dst);
    H5Sclose(msp);
    H5Dclose(ds);
    return 0;
}

/**
 * Read a 1-D scan-axis dataset at 'path' into a newly malloc'd double[n].
 * Fills zeros if absent.
 */
static double *h5_read_scan_col(hid_t f, const char *path, size_t n)
{
    double *buf = (double *)calloc(n, sizeof(double));
    if (!buf) return NULL;
    if (!h5_exists(f, path)) return buf;   /* return all-zeros */

    hid_t ds  = H5Dopen2(f, path, H5P_DEFAULT);
    hid_t fsp = H5Dget_space(ds);
    /* Use the actual dataset length (may differ from n in corrupt files) */
    hsize_t actual = 0;
    H5Sget_simple_extent_dims(fsp, &actual, NULL);
    hsize_t cnt = (actual < (hsize_t)n) ? actual : (hsize_t)n;

    hid_t msp = H5Screate_simple(1, &cnt, NULL);
    hsize_t off = 0;
    H5Sselect_hyperslab(fsp, H5S_SELECT_SET, &off, NULL, &cnt, NULL);
    H5Dread(ds, H5T_NATIVE_DOUBLE, msp, fsp, H5P_DEFAULT, buf);
    H5Sclose(msp); H5Sclose(fsp); H5Dclose(ds);
    return buf;
}

/**
 * Determine n_points from the first available scan dataset.
 * Tries TAS_COL_S2 first, then falls through all columns.
 */
static size_t detect_n_points(hid_t f)
{
    /* Prefer the s2 dataset; any 1-D scan dataset works */
    int preferred[] = { TAS_COL_S2, TAS_COL_COUNTS, TAS_COL_MONITOR,
                        TAS_COL_EI, TAS_COL_S1, -1 };
    for (int i = 0; preferred[i] >= 0; i++) {
        int c = preferred[i];
        if (!h5_exists(f, COL_META[c].hdf5_path)) continue;
        hid_t ds  = H5Dopen2(f, COL_META[c].hdf5_path, H5P_DEFAULT);
        hid_t sp  = H5Dget_space(ds);
        hsize_t n = 0;
        H5Sget_simple_extent_dims(sp, &n, NULL);
        H5Sclose(sp); H5Dclose(ds);
        if (n > 0) return (size_t)n;
    }
    /* Last resort: walk all columns */
    for (int c = 0; c < TAS_N_COLS; c++) {
        if (!h5_exists(f, COL_META[c].hdf5_path)) continue;
        hid_t ds = H5Dopen2(f, COL_META[c].hdf5_path, H5P_DEFAULT);
        hid_t sp = H5Dget_space(ds);
        hsize_t n = 0;
        H5Sget_simple_extent_dims(sp, &n, NULL);
        H5Sclose(sp); H5Dclose(ds);
        if (n > 0) return (size_t)n;
    }
    return 0;
}

/* =========================================================================
 * Read user string array from /entry/metadata/user
 * ========================================================================= */

static void read_users(hid_t f, tas_scan_info_t *info)
{
    const char *path = "/entry/metadata/user";
    info->n_users = 0;
    for (int i = 0; i < TAS_MAX_USERS; i++) info->users[i] = NULL;

    if (!h5_exists(f, path)) {
        info->users[0] = strdup("");
        info->n_users  = 0;
        return;
    }

    hid_t ds  = H5Dopen2(f, path, H5P_DEFAULT);
    hid_t sp  = H5Dget_space(ds);
    hid_t typ = H5Dget_type(ds);

    int ndims = H5Sget_simple_extent_ndims(sp);

    if (ndims == 0) {
        /* Scalar — single user string */
        char *s = h5_read_str(f, path, "");
        info->users[0] = s;
        info->n_users  = 1;
    } else {
        /* Array of strings */
        hsize_t count = 0;
        H5Sget_simple_extent_dims(sp, &count, NULL);
        if (count > TAS_MAX_USERS) count = TAS_MAX_USERS;

        if (H5Tis_variable_str(typ)) {
            char **bufs = (char **)calloc((size_t)count, sizeof(char *));
            hid_t mem = H5Tcopy(H5T_C_S1);
            H5Tset_size(mem, H5T_VARIABLE);
            H5Tset_cset(mem, H5T_CSET_UTF8);
            H5Tset_strpad(mem, H5T_STR_NULLTERM);
            H5Dread(ds, mem, H5S_ALL, H5S_ALL, H5P_DEFAULT, bufs);
            for (hsize_t i = 0; i < count; i++)
                info->users[i] = strdup(bufs[i] ? bufs[i] : "");
            H5Dvlen_reclaim(mem, sp, H5P_DEFAULT, bufs);
            free(bufs);
            H5Tclose(mem);
        } else {
            /* Fixed-length string array */
            size_t sz = H5Tget_size(typ);
            char *flat = (char *)calloc((size_t)count * (sz + 1), 1);
            H5Dread(ds, typ, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat);
            for (hsize_t i = 0; i < count; i++)
                info->users[i] = strdup(flat + i * sz);
            free(flat);
        }
        info->n_users = (int)count;
    }

    H5Tclose(typ); H5Sclose(sp); H5Dclose(ds);
}

/* =========================================================================
 * tas_read_hdf — main entry point
 * ========================================================================= */

tas_scan_t *tas_read_hdf(const char *filepath, int flags)
{
    if (!filepath) {
        set_error("tas_read_hdf: filepath is NULL");
        return NULL;
    }

    /* Suppress HDF5 default error printing */
    H5E_auto2_t old_func; void *old_data;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_data);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    hid_t fid = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fid < 0) {
        H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
        set_error("Cannot open file: %s", filepath);
        return NULL;
    }

    /* ── Allocate scan ────────────────────────────────────────────────── */
    tas_scan_t *scan = (tas_scan_t *)calloc(1, sizeof(tas_scan_t));
    if (!scan) { H5Fclose(fid); set_error("Out of memory"); return NULL; }

    /* ── 1.  Determine n_points ────────────────────────────────────────── */
    scan->n_points = detect_n_points(fid);
    if (scan->n_points == 0) {
        set_error("No scan data found in %s", filepath);
        free(scan); H5Fclose(fid);
        H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
        return NULL;
    }
    size_t n = scan->n_points;

    /* ── 2.  Scalar channels ───────────────────────────────────────────── */
    for (int c = 0; c < TAS_N_COLS; c++) {
        scan->data[c] = h5_read_scan_col(fid, COL_META[c].hdf5_path, n);
        if (!scan->data[c]) {
            set_error("Out of memory allocating column %d (%s)",
                      c, COL_META[c].name);
            tas_scan_free(scan);
            H5Fclose(fid);
            H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
            return NULL;
        }
    }

    /* ── 3.  PSD ───────────────────────────────────────────────────────── */
    const char *psd_path = "/entry/instrument/det_group/psd/data";
    scan->psd = NULL;
    scan->info.psd_present = false;

    if ((flags & TAS_LOAD_PSD) && h5_exists(fid, psd_path)) {
        size_t total = n * (size_t)TAS_PSD_FRAME_SIZE;
        double *psd  = (double *)calloc(total, sizeof(double));
        if (!psd) {
            set_error("Out of memory allocating PSD (%zu doubles)", total);
            tas_scan_free(scan);
            H5Fclose(fid);
            H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
            return NULL;
        }

        hid_t ds  = H5Dopen2(fid, psd_path, H5P_DEFAULT);
        hid_t fsp = H5Dget_space(ds);
        /* Read [n][NY][NX] — flatten via memory space */
        hsize_t mdims[3] = {(hsize_t)n, TAS_PSD_NY, TAS_PSD_NX};
        hid_t msp = H5Screate_simple(3, mdims, NULL);
        H5Dread(ds, H5T_NATIVE_DOUBLE, msp, fsp, H5P_DEFAULT, psd);
        H5Sclose(msp); H5Sclose(fsp); H5Dclose(ds);

        /* Check if all-zero (treat as absent, same as Python) */
        bool nonzero = false;
        for (size_t k = 0; k < total && !nonzero; k++)
            if (psd[k] != 0.0) nonzero = true;

        if (nonzero) {
            scan->psd = psd;
            scan->info.psd_present = true;
        } else {
            free(psd);  /* all-zero → return NULL like Python */
        }
    }

    /* ── 4.  String metadata ───────────────────────────────────────────── */
    tas_scan_info_t *info = &scan->info;

    /* Map table-driven string fields to struct members */
    for (int i = 0; META_STR[i].key != NULL; i++) {
        char *val = h5_read_str(fid, META_STR[i].path, "");
#define ASSIGN_STR(field) \
        if (strcmp(META_STR[i].key, #field) == 0) { info->field = val; continue; }
        ASSIGN_STR(start_time)
        ASSIGN_STR(end_time)
        ASSIGN_STR(title)
        ASSIGN_STR(command)
        ASSIGN_STR(filename)
        ASSIGN_STR(scanning_axis)
        ASSIGN_STR(scan_no)
        ASSIGN_STR(facility)
        ASSIGN_STR(source)
        ASSIGN_STR(instrument_name)
        ASSIGN_STR(local_contact)
        ASSIGN_STR(experiment_id)
        ASSIGN_STR(proposal_no)
        ASSIGN_STR(sample_name)
        ASSIGN_STR(sample_type)
        ASSIGN_STR(mono_crystal)
        ASSIGN_STR(ana_crystal)
        ASSIGN_STR(sense)
#undef ASSIGN_STR
        free(val);   /* unrecognised key — shouldn't happen */
    }

    /* Ensure all string pointers are non-NULL (strdup("") as fallback) */
#define ENSURE_STR(field) if (!info->field) info->field = strdup("")
    ENSURE_STR(start_time);  ENSURE_STR(end_time);    ENSURE_STR(title);
    ENSURE_STR(command);     ENSURE_STR(filename);     ENSURE_STR(scanning_axis);
    ENSURE_STR(scan_no);     ENSURE_STR(facility);     ENSURE_STR(source);
    ENSURE_STR(instrument_name); ENSURE_STR(local_contact); ENSURE_STR(experiment_id);
    ENSURE_STR(proposal_no); ENSURE_STR(sample_name);  ENSURE_STR(sample_type);
    ENSURE_STR(mono_crystal); ENSURE_STR(ana_crystal); ENSURE_STR(sense);
#undef ENSURE_STR

    /* Fallback: scanning_axis from NXdata @axes attribute */
    if (info->scanning_axis[0] == '\0') {
        const char *nxdata = "/entry/data";
        if (h5_exists(fid, nxdata)) {
            hid_t grp = H5Gopen2(fid, nxdata, H5P_DEFAULT);
            if (H5Aexists(grp, "axes") > 0) {
                hid_t at = H5Aopen(grp, "axes", H5P_DEFAULT);
                hid_t at_type = H5Aget_type(at);
                if (H5Tis_variable_str(at_type)) {
                    char *buf = NULL;
                    H5Aread(at, at_type, &buf);
                    free(info->scanning_axis);
                    info->scanning_axis = strdup(buf ? buf : "");
                    H5free_memory(buf);
                } else {
                    size_t sz = H5Tget_size(at_type);
                    char *buf = (char *)calloc(sz + 1, 1);
                    H5Aread(at, at_type, buf);
                    free(info->scanning_axis);
                    info->scanning_axis = strdup(buf);
                    free(buf);
                }
                H5Tclose(at_type); H5Aclose(at);
            }
            H5Gclose(grp);
        }
    }

    /* Users */
    read_users(fid, info);

    /* ── 5.  Float metadata ────────────────────────────────────────────── */
    info->sample_mosaic        = h5_read_double(fid, "/entry/sample/sample_mosaic",          0.0);
    info->distance_vs_mono     = h5_read_double(fid, "/entry/instrument/distance_vs_mono",   0.0);
    info->distance_mono_sample = h5_read_double(fid, "/entry/instrument/distance_mono_sample",2.0);
    info->distance_sample_ana  = h5_read_double(fid, "/entry/instrument/distance_sample_ana", 1.5);
    info->distance_ana_det     = h5_read_double(fid, "/entry/instrument/distance_ana_det",    0.5);

    /* ── 6.  Array metadata ────────────────────────────────────────────── */
    memset(info->sample_v1, 0, sizeof(info->sample_v1));
    memset(info->sample_v2, 0, sizeof(info->sample_v2));
    memset(info->unit_cell, 0, sizeof(info->unit_cell));
    memset(info->ub_matrix, 0, sizeof(info->ub_matrix));

    h5_read_double_array(fid, "/entry/sample/sample_v1",  info->sample_v1, 3);
    h5_read_double_array(fid, "/entry/sample/sample_v2",  info->sample_v2, 3);
    h5_read_double_array(fid, "/entry/sample/unit_cell",  info->unit_cell, 6);
    h5_read_double_array(fid, "/entry/sample/ub_matrix",  (double *)info->ub_matrix, 9);

    /* ── 7.  Derived fields ────────────────────────────────────────────── */
    info->num_points = (int)n;

    H5Fclose(fid);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
    return scan;
}

/* =========================================================================
 * tas_scan_free
 * ========================================================================= */

void tas_scan_free(tas_scan_t *scan)
{
    if (!scan) return;

    for (int c = 0; c < TAS_N_COLS; c++) {
        free(scan->data[c]);
        scan->data[c] = NULL;
    }
    free(scan->psd);

    tas_scan_info_t *info = &scan->info;
#define FREE_STR(f) free(info->f); info->f = NULL
    FREE_STR(start_time);   FREE_STR(end_time);    FREE_STR(title);
    FREE_STR(command);      FREE_STR(filename);     FREE_STR(scanning_axis);
    FREE_STR(scan_no);      FREE_STR(facility);     FREE_STR(source);
    FREE_STR(instrument_name); FREE_STR(local_contact); FREE_STR(experiment_id);
    FREE_STR(proposal_no);  FREE_STR(sample_name);  FREE_STR(sample_type);
    FREE_STR(mono_crystal); FREE_STR(ana_crystal);  FREE_STR(sense);
#undef FREE_STR

    for (int i = 0; i < TAS_MAX_USERS; i++) {
        free(info->users[i]);
        info->users[i] = NULL;
    }

    free(scan);
}

/* =========================================================================
 * Column access helpers
 * ========================================================================= */

const double *tas_col_ptr(const tas_scan_t *scan, tas_col_id_t col)
{
    if (!scan || col < 0 || col >= TAS_N_COLS) return NULL;
    return scan->data[col];
}

double tas_col_value(const tas_scan_t *scan, tas_col_id_t col, size_t idx)
{
    if (!scan || col < 0 || col >= TAS_N_COLS || idx >= scan->n_points)
        return 0.0;
    return scan->data[col][idx];
}

const double *tas_psd_frame(const tas_scan_t *scan, size_t idx)
{
    if (!scan || !scan->psd || idx >= scan->n_points) return NULL;
    return scan->psd + idx * (size_t)TAS_PSD_FRAME_SIZE;
}

/* =========================================================================
 * tas_print_summary — mirrors Python __main__ print block
 * ========================================================================= */

void tas_print_summary(const tas_scan_t *scan)
{
    if (!scan) { printf("(null scan)\n"); return; }
    const tas_scan_info_t *info = &scan->info;
    size_t n = scan->n_points;

    printf("\n%s\n", "=======================================================");
    printf("Points   : %zu\n", n);
    printf("Facility : %s\n", info->facility);
    printf("Instr.   : %s\n", info->instrument_name);
    printf("Sample   : %s\n", info->sample_name);
    printf("Title    : %s\n", info->title);
    printf("Axis     : %s\n", info->scanning_axis);
    printf("Start    : %s\n", info->start_time);
    printf("End      : %s\n", info->end_time);

    printf("\n── Scalar channels ────────────────────────────────────\n");

    /* Report range of the scanning axis */
    const char *ax = info->scanning_axis;
    int ax_col = tas_col_from_name(ax);
    if (ax_col >= 0 && scan->data[ax_col]) {
        double mn = scan->data[ax_col][0], mx = scan->data[ax_col][0];
        for (size_t i = 1; i < n; i++) {
            if (scan->data[ax_col][i] < mn) mn = scan->data[ax_col][i];
            if (scan->data[ax_col][i] > mx) mx = scan->data[ax_col][i];
        }
        printf("  %s range    : %.4f → %.4f\n", ax, mn, mx);
    }

    /* Counts peak */
    const double *cnt = scan->data[TAS_COL_COUNTS];
    if (cnt) {
        double peak = cnt[0]; size_t peak_i = 0;
        for (size_t i = 1; i < n; i++) if (cnt[i] > peak) { peak = cnt[i]; peak_i = i; }
        double ax_at_peak = (ax_col >= 0 && scan->data[ax_col])
                          ? scan->data[ax_col][peak_i] : 0.0;
        printf("  counts peak  : %.1f  (at %s=%.4f)\n", peak, ax, ax_at_peak);
    }

    printf("\n── Metadata ───────────────────────────────────────────\n");
    printf("  facility          : %s\n", info->facility);
    printf("  source            : %s\n", info->source);
    printf("  instrument        : %s\n", info->instrument_name);
    printf("  experiment_id     : %s\n", info->experiment_id);
    printf("  proposal_no       : %s\n", info->proposal_no);
    printf("  local_contact     : %s\n", info->local_contact);
    printf("  sample_name       : %s\n", info->sample_name);
    printf("  sample_type       : %s\n", info->sample_type);
    printf("  sample_mosaic     : %.4f arcmin\n", info->sample_mosaic);
    printf("  mono_crystal      : %s\n", info->mono_crystal);
    printf("  ana_crystal       : %s\n", info->ana_crystal);
    printf("  sense             : %s\n", info->sense);
    printf("  distance_vs_mono  : %.3f m\n", info->distance_vs_mono);
    printf("  distance_m_s      : %.3f m\n", info->distance_mono_sample);
    printf("  distance_s_a      : %.3f m\n", info->distance_sample_ana);
    printf("  distance_a_d      : %.3f m\n", info->distance_ana_det);
    printf("  unit_cell         : [%.3f %.3f %.3f  %.1f %.1f %.1f]\n",
           info->unit_cell[0], info->unit_cell[1], info->unit_cell[2],
           info->unit_cell[3], info->unit_cell[4], info->unit_cell[5]);
    printf("  ub_matrix[0]      : [%.4f %.4f %.4f]\n",
           info->ub_matrix[0][0], info->ub_matrix[0][1], info->ub_matrix[0][2]);
    printf("  ub_matrix[1]      : [%.4f %.4f %.4f]\n",
           info->ub_matrix[1][0], info->ub_matrix[1][1], info->ub_matrix[1][2]);
    printf("  ub_matrix[2]      : [%.4f %.4f %.4f]\n",
           info->ub_matrix[2][0], info->ub_matrix[2][1], info->ub_matrix[2][2]);
    printf("  sample_v1         : [%.3f %.3f %.3f]\n",
           info->sample_v1[0], info->sample_v1[1], info->sample_v1[2]);
    printf("  sample_v2         : [%.3f %.3f %.3f]\n",
           info->sample_v2[0], info->sample_v2[1], info->sample_v2[2]);

    if (info->n_users > 0) {
        printf("  users             : ");
        for (int i = 0; i < info->n_users; i++)
            printf("%s%s", info->users[i], i < info->n_users - 1 ? ", " : "");
        printf("\n");
    }

    printf("\n── PSD ────────────────────────────────────────────────\n");
    if (scan->psd) {
        size_t total = n * (size_t)TAS_PSD_FRAME_SIZE;
        double sum = 0.0;
        for (size_t k = 0; k < total; k++) sum += scan->psd[k];
        double frame0_sum = 0.0;
        for (int k = 0; k < TAS_PSD_FRAME_SIZE; k++) frame0_sum += scan->psd[k];
        printf("  shape       : (%zu, %d, %d)\n", n, TAS_PSD_NY, TAS_PSD_NX);
        printf("  frame[0] sum: %.0f\n", frame0_sum);
        printf("  total sum   : %.0f\n", sum);
    } else {
        printf("  Not present or all-zero.\n");
    }
    printf("\n");
}

/* =========================================================================
 * tas_export_text — mirrors Python export_to_text()
 * ========================================================================= */

/* Column display order matching the Python function */
static const int TEXT_COL_PRIORITY[] = {
    /* Fixed start (Q/E) */
    TAS_COL_QH, TAS_COL_QK, TAS_COL_QL,
    TAS_COL_EN, TAS_COL_EI, TAS_COL_EF,
    /* Fixed middle (key motors) */
    TAS_COL_M1, TAS_COL_M2,
    TAS_COL_S1, TAS_COL_S2,
    TAS_COL_A1, TAS_COL_A2,
    /* Fixed end (detectors) */
    TAS_COL_COUNTS, TAS_COL_MONITOR,
    /* Remaining: source, mono, ana, slits, collimators, sample motors, env */
    TAS_COL_VS_LEFT, TAS_COL_VS_RIGHT,
    TAS_COL_MONOVF,  TAS_COL_MONOHF,  TAS_COL_MONOTILT, TAS_COL_MONOTRANS,
    TAS_COL_PS_LEFT, TAS_COL_PS_RIGHT, TAS_COL_PS_TOP,   TAS_COL_PS_BOTTOM,
    TAS_COL_PA_LEFT, TAS_COL_PA_RIGHT, TAS_COL_PA_TOP,   TAS_COL_PA_BOTTOM,
    TAS_COL_COL1,    TAS_COL_COL2,     TAS_COL_COL3,     TAS_COL_COL4,
    TAS_COL_COLL_ALPHA1, TAS_COL_COLL_ALPHA2,
    TAS_COL_COLL_ALPHA3, TAS_COL_COLL_ALPHA4,
    TAS_COL_SGU,  TAS_COL_SGL,  TAS_COL_STU,  TAS_COL_STL,
    TAS_COL_ANAVF, TAS_COL_ANAHF, TAS_COL_ANATILT, TAS_COL_ANATRANS,
    TAS_COL_TEMP1, TAS_COL_TEMP2, TAS_COL_TEMP3,   TAS_COL_TEMP4,
    TAS_COL_MFIELD, TAS_COL_EFIELD, TAS_COL_PRESSURE,
    TAS_COL_CRYO_HE, TAS_COL_CRYO_NEEDLE,
};
static const int N_TEXT_COLS = (int)(sizeof(TEXT_COL_PRIORITY) / sizeof(int));

/* Write the UB matrix on one line as Python does */
static void write_ub_oneline(FILE *out, const double ub[3][3])
{
    fprintf(out,
            "[[%.4f %.4f %.4f], [%.4f %.4f %.4f], [%.4f %.4f %.4f]]",
            ub[0][0], ub[0][1], ub[0][2],
            ub[1][0], ub[1][1], ub[1][2],
            ub[2][0], ub[2][1], ub[2][2]);
}

int tas_export_text(const tas_scan_t *scan, const char *output_filepath)
{
    if (!scan || !output_filepath) {
        set_error("tas_export_text: NULL argument");
        return TAS_ERR_ARG;
    }

    FILE *out = fopen(output_filepath, "w");
    if (!out) {
        set_error("Cannot open output file: %s (%s)",
                  output_filepath, strerror(errno));
        return TAS_ERR_GENERAL;
    }

    const tas_scan_info_t *info = &scan->info;
    size_t n = scan->n_points;

    /* ── 1. Metadata header ─────────────────────────────────────────────── */

    fprintf(out, "# TAS_NeXus_Version = %s\n",   "");  /* not stored in scan_info */
    fprintf(out, "# software_version = %s\n",    "");
    fprintf(out, "# filename = %s\n",             info->filename);
    fprintf(out, "# facility = %s\n",             info->facility);
    fprintf(out, "# source = %s\n",               info->source);
    fprintf(out, "# instrument = %s\n",           info->instrument_name);
    fprintf(out, "# experiment_id = %s\n",        info->experiment_id);
    fprintf(out, "# proposal = %s\n",             info->proposal_no);

    /* users list */
    fprintf(out, "# user(s) = [");
    for (int i = 0; i < info->n_users; i++) {
        fprintf(out, "%s%s", info->users[i],
                i < info->n_users - 1 ? ", " : "");
    }
    fprintf(out, "]\n");

    fprintf(out, "# local_contact = %s\n",        info->local_contact);
    fprintf(out, "# mono_crystal = %s\n",          info->mono_crystal);
    fprintf(out, "# ana_crystal = %s\n",           info->ana_crystal);
    fprintf(out, "# sense = %s\n",                 info->sense);
    fprintf(out, "# distance_vs_mono = %.4g\n",    info->distance_vs_mono);
    fprintf(out, "# distance_mono_sample = %.4g\n",info->distance_mono_sample);
    fprintf(out, "# distance_sample_ana = %.4g\n", info->distance_sample_ana);
    fprintf(out, "# distance_ana_det = %.4g\n",    info->distance_ana_det);
    fprintf(out, "# sample_name = %s\n",           info->sample_name);
    fprintf(out, "# sample_type = %s\n",           info->sample_type);
    fprintf(out, "# sample_mosaic = %.4g\n",        info->sample_mosaic);
    fprintf(out, "# sample_v1 = [%.4g %.4g %.4g]\n",
            info->sample_v1[0], info->sample_v1[1], info->sample_v1[2]);
    fprintf(out, "# sample_v2 = [%.4g %.4g %.4g]\n",
            info->sample_v2[0], info->sample_v2[1], info->sample_v2[2]);
    fprintf(out, "# unit_cell = [%.4g %.4g %.4g %.4g %.4g %.4g]\n",
            info->unit_cell[0], info->unit_cell[1], info->unit_cell[2],
            info->unit_cell[3], info->unit_cell[4], info->unit_cell[5]);
    fprintf(out, "# ub_matrix = ");
    write_ub_oneline(out, info->ub_matrix);
    fprintf(out, "\n");
    fprintf(out, "# scan_no = %s\n",          info->scan_no);
    fprintf(out, "# scantitle = %s\n",         info->title);
    fprintf(out, "# command = %s\n",           info->command);
    fprintf(out, "# start_time = %s\n",        info->start_time);
    fprintf(out, "# end_time = %s\n",          info->end_time);
    fprintf(out, "# scanning_axis = %s\n",     info->scanning_axis);
    fprintf(out, "# \n");

    /* ── 2. Build the ordered column list ──────────────────────────────── */
    /* The scanning axis goes first, then the priority list (deduped) */

    int ordered[TAS_N_COLS];
    bool used[TAS_N_COLS];
    int n_ordered = 0;
    memset(used, 0, sizeof(used));

    /* Scanning axis first */
    int ax_col = tas_col_from_name(info->scanning_axis);
    if (ax_col >= 0 && !used[ax_col]) {
        ordered[n_ordered++] = ax_col;
        used[ax_col] = true;
    }

    /* Priority list */
    for (int i = 0; i < N_TEXT_COLS; i++) {
        int c = TEXT_COL_PRIORITY[i];
        if (!used[c]) {
            ordered[n_ordered++] = c;
            used[c] = true;
        }
    }

    /* Any remaining columns not in priority list */
    for (int c = 0; c < TAS_N_COLS; c++) {
        if (!used[c]) {
            ordered[n_ordered++] = c;
            used[c] = true;
        }
    }

    /* ── 3. Column header ───────────────────────────────────────────────── */
    /* Python uses col_space=16, right-justified; 'counts' → 'detector' */
#define COL_W 16
    fprintf(out, "# ");
    fprintf(out, "%*s", COL_W, "Pt.");
    for (int i = 0; i < n_ordered; i++) {
        int c = ordered[i];
        const char *name = (c == TAS_COL_COUNTS) ? "detector" : COL_META[c].name;
        fprintf(out, "%*s", COL_W, name);
    }
    fprintf(out, "\n");

    /* ── 4. Data rows ────────────────────────────────────────────────────── */
    for (size_t row = 0; row < n; row++) {
        fprintf(out, "%*zu", COL_W, row + 1);
        for (int i = 0; i < n_ordered; i++) {
            int c = ordered[i];
            double v = scan->data[c] ? scan->data[c][row] : 0.0;
            fprintf(out, "%*.*f", COL_W, 4, v);
        }
        fprintf(out, "\n");
    }
#undef COL_W

    /* ── 5. PSD slices ──────────────────────────────────────────────────── */
    if (scan->psd && info->psd_present) {
        fprintf(out, "\n################PSD Data################\n");
        for (size_t fr = 0; fr < n; fr++) {
            fprintf(out, "#SLICE_%03zu\n", fr + 1);
            const double *frame = scan->psd + fr * (size_t)TAS_PSD_FRAME_SIZE;
            for (int row = 0; row < TAS_PSD_NY; row++) {
                for (int col = 0; col < TAS_PSD_NX; col++) {
                    fprintf(out, "%12.1f", frame[row * TAS_PSD_NX + col]);
                    if (col < TAS_PSD_NX - 1) fputc(' ', out);
                }
                fputc('\n', out);
            }
            fputc('\n', out);
        }
    }

    fclose(out);
    return TAS_OK;
}
