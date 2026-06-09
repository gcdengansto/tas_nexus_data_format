/**
 * @file tas_nexus_reader.h
 * @brief TAS (Triple-Axis Spectrometer) NeXus/HDF5 reader — C public API
 *
 * Reads files produced by tas_nexus_writer (C or Python) and exposes the
 * same three logical objects returned by the Python load_from_hdf():
 *
 *   1.  Scalar scan data   — one double array per motor/detector channel,
 *                            accessible by tas_col_id_t index or name.
 *   2.  PSD data           — optional (n_points × 128 × 128) double array.
 *   3.  Scan metadata      — all string, float, and array fields from the
 *                            NeXus file hierarchy.
 *
 * Additionally, an export function (tas_export_text) reproduces the text
 * dump produced by the Python export_to_text() function.
 *
 * ---------------------------------------------------------------------------
 * QUICK START
 * ---------------------------------------------------------------------------
 *
 *   tas_scan_t *scan = tas_read_hdf("my_scan.h5", TAS_LOAD_PSD);
 *   if (!scan) { fprintf(stderr, "%s\n", tas_last_error()); return 1; }
 *
 *   // Scalar channel access
 *   size_t n = scan->n_points;
 *   const double *s2     = tas_col_ptr(scan, TAS_COL_S2);
 *   const double *counts = tas_col_ptr(scan, TAS_COL_COUNTS);
 *
 *   // Metadata
 *   printf("Facility : %s\n", scan->info.facility);
 *   printf("Sample   : %s\n", scan->info.sample_name);
 *   printf("Points   : %zu\n", n);
 *
 *   // PSD (may be NULL if all-zero or load_psd was TAS_SKIP_PSD)
 *   if (scan->psd) {
 *       // frame i, pixel (row r, col c):
 *       double v = scan->psd[(size_t)i * TAS_PSD_FRAME_SIZE + r*TAS_PSD_NX + c];
 *   }
 *
 *   // Optional text export
 *   tas_export_text(scan, "output.txt");
 *
 *   tas_scan_free(scan);
 *
 * ---------------------------------------------------------------------------
 * BUILDING
 * ---------------------------------------------------------------------------
 *
 *   h5cc -O2 my_app.c tas_nexus_reader.c tas_nexus_writer.c -lm -o my_app
 *
 * The reader shares the column-metadata table in tas_nexus_writer.c and
 * therefore must be linked together with it.  The writer object itself is
 * not required at runtime if you only read files; link order does not matter.
 *
 * ---------------------------------------------------------------------------
 * MEMORY OWNERSHIP
 * ---------------------------------------------------------------------------
 *
 * All memory returned inside a tas_scan_t is owned by that object.
 * Call tas_scan_free() exactly once when done; after that, all pointers
 * into the struct are invalid.
 */

#ifndef TAS_NEXUS_READER_H
#define TAS_NEXUS_READER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "tas_nexus_writer.h"   /* tas_col_id_t, TAS_N_COLS, TAS_PSD_* */

#include <stddef.h>
#include <stdbool.h>

/* =========================================================================
 * Load flags (passed to tas_read_hdf)
 * ========================================================================= */

/** Load the full PSD array from disk (default / recommended). */
#define TAS_LOAD_PSD   1
/** Skip reading PSD data (faster for scalar-only analysis).   */
#define TAS_SKIP_PSD   0

/* =========================================================================
 * Scan metadata structure
 *
 * All char* fields point into memory owned by the parent tas_scan_t.
 * Do NOT free them individually.
 * ========================================================================= */

/** Maximum number of user strings stored in tas_scan_info_t.users[]. */
#define TAS_MAX_USERS 16

/**
 * @brief Scalar and array metadata read from the NeXus file.
 *
 * Mirrors the scan_info dict returned by the Python load_from_hdf().
 * Every string field is a heap-allocated C string owned by the parent
 * tas_scan_t; it is freed by tas_scan_free().
 */
typedef struct {
    /* ── Scan identification ────────────────────────────────────────── */
    char   *start_time;        /**< ISO-8601, e.g. "2024-01-15T08:30:00Z" */
    char   *end_time;
    char   *title;
    char   *command;           /**< Acquisition command string             */
    char   *filename;          /**< Original DAQ filename                  */
    char   *scanning_axis;     /**< e.g. "s2"                             */
    char   *scan_no;

    /* ── Facility / instrument ──────────────────────────────────────── */
    char   *facility;
    char   *source;
    char   *instrument_name;
    char   *local_contact;
    char   *experiment_id;
    char   *proposal_no;

    /* ── Users ──────────────────────────────────────────────────────── */
    char   *users[TAS_MAX_USERS]; /**< NULL-terminated list of user strings */
    int     n_users;

    /* ── Sample ─────────────────────────────────────────────────────── */
    char   *sample_name;
    char   *sample_type;
    double  sample_mosaic;      /**< minutes of arc                         */
    double  sample_v1[3];       /**< scattering-plane vector 1              */
    double  sample_v2[3];       /**< scattering-plane vector 2              */
    double  unit_cell[6];       /**< a, b, c, α, β, γ  (Å, °)             */
    double  ub_matrix[3][3];    /**< UB orientation matrix                  */

    /* ── Optics ─────────────────────────────────────────────────────── */
    char   *mono_crystal;       /**< e.g. "PG"                             */
    char   *ana_crystal;
    char   *sense;              /**< e.g. "+-+"                            */

    /* ── Distances (metres) ─────────────────────────────────────────── */
    double  distance_vs_mono;
    double  distance_mono_sample;
    double  distance_sample_ana;
    double  distance_ana_det;

    /* ── Derived flags ──────────────────────────────────────────────── */
    int     num_points;         /**< Total scan points read                 */
    bool    psd_present;        /**< true if PSD array contains non-zero data */
} tas_scan_info_t;

/* =========================================================================
 * Top-level scan structure
 * ========================================================================= */

/**
 * @brief Complete in-memory representation of one TAS NeXus scan.
 *
 * Obtain with tas_read_hdf(); release with tas_scan_free().
 */
typedef struct {
    /** Number of scan points. */
    size_t  n_points;

    /**
     * Scalar channel data: data[col][point_index].
     * Indexed by tas_col_id_t.  Channels not present in the file are
     * filled with zeros.  Shape: [TAS_N_COLS][n_points].
     * Access via tas_col_ptr() for convenience and bounds checking.
     */
    double *data[TAS_N_COLS];

    /**
     * PSD data: flat row-major array [n_points][TAS_PSD_NY][TAS_PSD_NX].
     * NULL if TAS_SKIP_PSD was passed or the dataset was all-zero.
     * Total elements: n_points * TAS_PSD_FRAME_SIZE.
     */
    double *psd;

    /** All metadata fields. */
    tas_scan_info_t info;
} tas_scan_t;

/* =========================================================================
 * Public API
 * ========================================================================= */

/**
 * @brief Read a TAS NeXus HDF5 file into memory.
 *
 * @param filepath   Path to the .tas.nxs.h5 file.
 * @param flags      TAS_LOAD_PSD or TAS_SKIP_PSD.
 *
 * @return Heap-allocated tas_scan_t on success; NULL on failure.
 *         Call tas_last_error() for a description of the failure.
 *         Free the result with tas_scan_free().
 */
tas_scan_t *tas_read_hdf(const char *filepath, int flags);

/**
 * @brief Free all memory associated with a scan returned by tas_read_hdf().
 *
 * After this call @p scan is invalid; do not dereference it.
 * Passing NULL is safe (no-op).
 */
void tas_scan_free(tas_scan_t *scan);

/**
 * @brief Return a pointer to the data array for column @p col.
 *
 * The returned pointer is valid for scan->n_points doubles.
 *
 * @return Non-NULL pointer (possibly all-zeros if the channel was absent),
 *         or NULL if @p col is out of range or @p scan is NULL.
 */
const double *tas_col_ptr(const tas_scan_t *scan, tas_col_id_t col);

/**
 * @brief Return the value of column @p col at scan point @p idx.
 *
 * @return The value, or 0.0 if arguments are out of range.
 */
double tas_col_value(const tas_scan_t *scan, tas_col_id_t col, size_t idx);

/**
 * @brief Return a pointer to one PSD frame (128×128 doubles, row-major).
 *
 * @param scan   Scan returned by tas_read_hdf().
 * @param idx    Scan-point index (0-based).
 *
 * @return Pointer to TAS_PSD_FRAME_SIZE doubles, or NULL if PSD is absent
 *         or @p idx is out of range.
 */
const double *tas_psd_frame(const tas_scan_t *scan, size_t idx);

/**
 * @brief Export the scan to a human-readable text file.
 *
 * Reproduces the format of the Python export_to_text() function:
 *   - "#  key = value" metadata header lines
 *   - A tab-separated data table with one row per scan point
 *   - Optional PSD slices appended at the end
 *
 * @param scan            Scan to export.
 * @param output_filepath Destination file path (created or overwritten).
 *
 * @return TAS_OK on success, negative TAS_ERR_* on failure.
 */
int tas_export_text(const tas_scan_t *scan, const char *output_filepath);

/**
 * @brief Print a summary of the scan to stdout (diagnostic helper).
 *
 * Reproduces the print block in the Python __main__ self-test.
 */
void tas_print_summary(const tas_scan_t *scan);

#ifdef __cplusplus
}
#endif

#endif /* TAS_NEXUS_READER_H */
