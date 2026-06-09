/**
 * @file tas_nexus_writer.h
 * @brief TAS (Triple-Axis Spectrometer) NeXus/HDF5 writer — C public API
 *
 * Produces NeXus-compliant HDF5 files that match the structure of the
 * reference Python implementation (tas_nexus_writer_v1.py).
 *
 * Two write modes are provided:
 *
 *  1.  Batch   — supply all scan points at once via tas_save_hdf().
 *  2.  Point-by-point — open a file with tas_pbp_open(), write one point
 *      per call with tas_pbp_write_point(), then close with tas_pbp_close().
 *
 * ---------------------------------------------------------------------------
 * QUICK START — BATCH
 * ---------------------------------------------------------------------------
 *
 *   tas_scan_data_t  scan  = {0};
 *   tas_metadata_t   meta  = {0};
 *   tas_scalar_col_t cols[2];
 *
 *   double s2[101], counts[101];  // fill arrays ...
 *   cols[0] = (tas_scalar_col_t){ TAS_COL_S2,     s2     };
 *   cols[1] = (tas_scalar_col_t){ TAS_COL_COUNTS,  counts };
 *
 *   scan.n_points  = 101;
 *   scan.columns   = cols;
 *   scan.n_columns = 2;
 *   scan.psd_data  = NULL;   // or (n_points × 128 × 128) row-major double*
 *
 *   meta.facility        = "ANSTO";
 *   meta.instrument_name = "TAIPAN";
 *
 *   int rc = tas_save_hdf("scan.h5", &scan, &meta, TAS_OVERWRITE);
 *
 * ---------------------------------------------------------------------------
 * QUICK START — POINT-BY-POINT
 * ---------------------------------------------------------------------------
 *
 *   tas_pbp_handle_t *h = tas_pbp_open("scan.h5", 101, &meta);
 *
 *   for (int i = 0; i < 101; i++) {
 *       tas_point_t pt = {0};
 *       pt.values[TAS_COL_S2]     = s2_pos[i];
 *       pt.values[TAS_COL_COUNTS] = cnt[i];
 *       pt.psd_frame = NULL;   // or double[128][128]
 *       tas_pbp_write_point(h, &pt);
 *   }
 *   tas_pbp_close(h);
 *
 * ---------------------------------------------------------------------------
 * BUILDING
 * ---------------------------------------------------------------------------
 *
 *   gcc -O2 -c tas_nexus_writer.c $(h5cc -show | grep -oP '\-[IL]\S+') \
 *       -o tas_nexus_writer.o
 *
 * Or simply let h5cc drive the compilation of your project:
 *
 *   h5cc -O2 my_app.c tas_nexus_writer.c -o my_app
 *
 * ---------------------------------------------------------------------------
 * THREAD SAFETY
 * ---------------------------------------------------------------------------
 *
 * The library itself does NOT use global state (other than passing the
 * TAS_TRACK_ORDER flag to HDF5).  Concurrent writes to DIFFERENT files with
 * DIFFERENT handles are safe provided HDF5 was built with thread safety.
 * Do NOT share a tas_pbp_handle_t across threads.
 */

#ifndef TAS_NEXUS_WRITER_H
#define TAS_NEXUS_WRITER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>   /* size_t  */
#include <stdint.h>   /* int64_t */

/* =========================================================================
 * Return codes
 * ========================================================================= */

/** Returned by all API functions on success. */
#define TAS_OK              0
/** Generic error — call tas_last_error() for a human-readable message.    */
#define TAS_ERR_GENERAL    -1
/** Bad argument (NULL pointer, out-of-range index, …).                    */
#define TAS_ERR_ARG        -2
/** HDF5 library call failed — underlying error printed to stderr by HDF5. */
#define TAS_ERR_HDF5       -3
/** File already exists and TAS_NO_OVERWRITE was requested.                */
#define TAS_ERR_EXISTS     -4
/** Scalar column array pointer is NULL for a column that was requested.   */
#define TAS_ERR_NULL_DATA  -5

/* =========================================================================
 * Scalar column identifiers
 *
 * Used as indices into tas_point_t.values[] and as the 'column' field of
 * tas_scalar_col_t.  Keep TAS_N_COLS last — it gives the array size.
 * ========================================================================= */

typedef enum {
    /* Sample motors */
    TAS_COL_S1 = 0, TAS_COL_S2, TAS_COL_SGU, TAS_COL_SGL,
    TAS_COL_STU, TAS_COL_STL,
    /* Virtual Q/E */
    TAS_COL_QH, TAS_COL_QK, TAS_COL_QL,
    TAS_COL_EI, TAS_COL_EF, TAS_COL_EN,
    /* Monochromator */
    TAS_COL_M1, TAS_COL_M2,
    TAS_COL_MONOVF, TAS_COL_MONOHF, TAS_COL_MONOTILT, TAS_COL_MONOTRANS,
    /* Analyser */
    TAS_COL_A1, TAS_COL_A2,
    TAS_COL_ANAVF, TAS_COL_ANAHF, TAS_COL_ANATILT, TAS_COL_ANATRANS,
    /* Source */
    TAS_COL_VS_LEFT, TAS_COL_VS_RIGHT,
    /* Pre-sample slits */
    TAS_COL_PS_LEFT, TAS_COL_PS_RIGHT, TAS_COL_PS_TOP, TAS_COL_PS_BOTTOM,
    /* Post-analyser slits */
    TAS_COL_PA_LEFT, TAS_COL_PA_RIGHT, TAS_COL_PA_TOP, TAS_COL_PA_BOTTOM,
    /* Collimators — physical motors */
    TAS_COL_COL1, TAS_COL_COL2, TAS_COL_COL3, TAS_COL_COL4,
    /* Collimators — angular values */
    TAS_COL_COLL_ALPHA1, TAS_COL_COLL_ALPHA2,
    TAS_COL_COLL_ALPHA3, TAS_COL_COLL_ALPHA4,
    /* Sample environment */
    TAS_COL_TEMP1, TAS_COL_TEMP2, TAS_COL_TEMP3, TAS_COL_TEMP4,
    TAS_COL_MFIELD, TAS_COL_EFIELD, TAS_COL_PRESSURE,
    TAS_COL_CRYO_HE, TAS_COL_CRYO_NEEDLE,
    /* Detectors */
    TAS_COL_COUNTS, TAS_COL_MONITOR,

    TAS_N_COLS   /**< Sentinel — number of scalar columns */
} tas_col_id_t;

/* =========================================================================
 * PSD geometry
 * ========================================================================= */

/** Width (fast axis) of the position-sensitive detector in pixels. */
#define TAS_PSD_NX  128
/** Height (slow axis) of the position-sensitive detector in pixels. */
#define TAS_PSD_NY  128
/** Total pixels per PSD frame (NX × NY). */
#define TAS_PSD_FRAME_SIZE  (TAS_PSD_NX * TAS_PSD_NY)

/* =========================================================================
 * Data structures
 * ========================================================================= */

/**
 * @brief One (column, data-pointer) pair used in the batch API.
 *
 * The data pointer must point to an array of at least n_points doubles,
 * stored contiguously in row-major order (C array).
 */
typedef struct {
    tas_col_id_t  column;   /**< Which scalar column this array provides. */
    const double *data;     /**< Array of length scan.n_points.           */
} tas_scalar_col_t;

/**
 * @brief Complete scan data for the batch API (tas_save_hdf).
 */
typedef struct {
    /** Number of scan points (rows). */
    size_t              n_points;

    /** Array of (column, pointer) pairs, length n_columns. */
    const tas_scalar_col_t *columns;
    size_t              n_columns;

    /**
     * PSD data: row-major array of shape [n_points][TAS_PSD_NY][TAS_PSD_NX].
     * Pass NULL to leave the PSD dataset filled with zeros.
     * Total element count: n_points * TAS_PSD_FRAME_SIZE.
     */
    const double       *psd_data;
} tas_scan_data_t;

/**
 * @brief Single-point data for the point-by-point API.
 */
typedef struct {
    /**
     * Motor/detector values indexed by tas_col_id_t.
     * Unset entries default to 0.0.
     */
    double        values[TAS_N_COLS];

    /**
     * PSD frame: row-major [TAS_PSD_NY][TAS_PSD_NX] array, or NULL.
     * Total element count: TAS_PSD_FRAME_SIZE.
     */
    const double *psd_frame;
} tas_point_t;

/**
 * @brief Experiment / instrument metadata.
 *
 * All fields are optional — NULL or 0.0 causes a sensible NeXus default to
 * be written.  String pointers are shallow-copied (kept only during the
 * call); they need not outlive the write call.
 *
 * users[] is a NULL-terminated array of C strings, e.g.:
 *     const char *u[] = { "Alice", "Bob", NULL };
 *     meta.users = u;
 */
typedef struct {
    /* Facility / instrument */
    const char  *facility;           /**< e.g. "ANSTO"              */
    const char  *source;             /**< e.g. "OPAL Reactor"       */
    const char  *instrument_name;    /**< e.g. "TAIPAN"             */
    const char  *software_version;   /**< e.g. "1.0.0"              */
    const char  *tas_nexus_version;  /**< e.g. "0.9.1"              */

    /* Proposal */
    const char  *experiment_id;
    const char  *proposal_no;
    const char  *local_contact;

    /* Scan */
    const char  *scan_no;
    const char  *title;
    const char  *command;
    const char  *filename;           /**< Original data-acquisition filename */
    const char  *scanning_axis;      /**< e.g. "s2"                 */
    const char  *start_time;         /**< ISO-8601 or NULL → now    */
    const char  *end_time;           /**< ISO-8601 or NULL → now    */

    /* Sample */
    const char  *sample_name;
    const char  *sample_type;        /**< e.g. "crystal"            */
    double       sample_mosaic;      /**< arcmin, 0 → written as 0  */
    double       sample_v1[3];       /**< scattering plane v1 vector */
    double       sample_v2[3];       /**< scattering plane v2 vector */
    double       unit_cell[6];       /**< a,b,c,α,β,γ (Å, °)       */
    double       ub_matrix[3][3];    /**< UB orientation matrix      */

    /* Optics */
    const char  *mono_crystal;       /**< e.g. "PG"                 */
    const char  *ana_crystal;        /**< e.g. "PG"                 */
    const char  *sense;              /**< e.g. "+-+"                */

    /* Distances (metres) */
    double       distance_vs_mono;
    double       distance_mono_sample;
    double       distance_sample_ana;
    double       distance_ana_det;

    /* Users — NULL-terminated array of C strings */
    const char * const *users;       /**< e.g. { "Alice", "Bob", NULL } */
} tas_metadata_t;

/**
 * @brief Opaque handle for point-by-point writing.
 *
 * Obtain with tas_pbp_open(), release with tas_pbp_close().
 * Do not access members directly.
 */
typedef struct tas_pbp_handle tas_pbp_handle_t;

/* =========================================================================
 * File-creation flags
 * ========================================================================= */

/** Silently replace an existing file (default-safe behaviour). */
#define TAS_OVERWRITE     1
/** Fail with TAS_ERR_EXISTS if the file already exists.        */
#define TAS_NO_OVERWRITE  0

/* =========================================================================
 * Public API — Batch
 * ========================================================================= */

/**
 * @brief Write an entire scan to a new NeXus HDF5 file in one call.
 *
 * @param filepath   Output file path (created or replaced).
 * @param scan       Scan data (motor arrays + optional PSD array).
 * @param meta       Metadata (may be NULL for all defaults).
 * @param flags      TAS_OVERWRITE or TAS_NO_OVERWRITE.
 *
 * @return TAS_OK on success, negative TAS_ERR_* code on failure.
 *         Call tas_last_error() for a description.
 */
int tas_save_hdf(const char         *filepath,
                 const tas_scan_data_t *scan,
                 const tas_metadata_t  *meta,
                 int                    flags);

/* =========================================================================
 * Public API — Point-by-point
 * ========================================================================= */

/**
 * @brief Create a new NeXus HDF5 file and return an open write handle.
 *
 * The file is created immediately with all datasets pre-allocated to
 * n_points rows, matching the Python reference implementation.
 *
 * @param filepath   Output file path (always created fresh — any existing
 *                   file at this path is replaced).
 * @param n_points   Total number of scan points that will be written.
 *                   Must be ≥ 1.
 * @param meta       Metadata written into static datasets at creation time.
 *                   May be NULL for all defaults.
 *
 * @return Non-NULL handle on success; NULL on failure (see tas_last_error()).
 */
tas_pbp_handle_t *tas_pbp_open(const char           *filepath,
                                size_t                n_points,
                                const tas_metadata_t *meta);

/**
 * @brief Append one scan point to a handle opened with tas_pbp_open().
 *
 * Points are written in the order this function is called (0-based index
 * incremented automatically).  Writing beyond n_points is an error.
 *
 * @param handle     Open write handle.
 * @param point      Values + optional PSD frame for this point.
 *
 * @return 1-based current point index on success (same as Python return),
 *         or negative TAS_ERR_* code on failure.
 */
int tas_pbp_write_point(tas_pbp_handle_t  *handle,
                        const tas_point_t *point);

/**
 * @brief Flush, stamp end_time, and close a point-by-point handle.
 *
 * After this call @p handle is invalid; do not use it again.
 *
 * @return TAS_OK on success, negative TAS_ERR_* code on failure.
 */
int tas_pbp_close(tas_pbp_handle_t *handle);

/* =========================================================================
 * Diagnostics
 * ========================================================================= */

/**
 * @brief Return a human-readable description of the last error.
 *
 * The returned pointer is valid until the next library call from the same
 * thread.  Thread-local on platforms that support it; otherwise a shared
 * static buffer (not safe for concurrent use from multiple threads).
 */
const char *tas_last_error(void);

/**
 * @brief Map a column ID to its canonical name string (e.g. "s2").
 *
 * @return The column name, or NULL if col is out of range.
 */
const char *tas_col_name(tas_col_id_t col);

/**
 * @brief Map a column name string to its ID.
 *
 * @return The column ID, or -1 if not found.
 */
int tas_col_from_name(const char *name);

#ifdef __cplusplus
}
#endif

#endif /* TAS_NEXUS_WRITER_H */
