/**
 * @file example.c
 * @brief TAS NeXus writer — C self-test matching the Python __main__ block.
 *
 * Compiles with:
 *   h5cc -O2 -Wall example.c tas_nexus_writer.c -lm -o example
 *
 * Produces two NeXus files:
 *   test_batch_c.tas.nxs.h5      — written in one shot (tas_save_hdf)
 *   test_pbp_c.tas.nxs.h5        — written point-by-point (tas_pbp_*)
 *
 * After running, verify with h5ls / HDFView / Python:
 *   python3 -c "
 *   import h5py, numpy as np
 *   a = h5py.File('test_batch_c.tas.nxs.h5','r')
 *   b = h5py.File('test_pbp_c.tas.nxs.h5','r')
 *   p = '/entry/sample/s2/value'
 *   print('Match:', np.allclose(a[p][()], b[p][()]))
 *   "
 */

#define _USE_MATH_DEFINES  /* MSVC */
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

#include "tas_nexus_writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N_POINTS  101
#define S2_START  (-65.0)
#define S2_END    (-61.0)

/* Simple PRNG (xoshiro128+ substitute for portability) */
static unsigned long rng_state = 42;
static double rng_normal(void)
{
    /* Box-Muller (cheap, good enough for test data) */
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    double u = ((rng_state >> 33) + 0.5) / (double)(1UL << 31);
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    double v = ((rng_state >> 33) + 0.5) / (double)(1UL << 31);
    return sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
}
static double rng_poisson_approx(double lam)
{
    /* Normal approximation to Poisson for lam >> 1 */
    double v = lam + sqrt(lam) * rng_normal();
    return v < 0.0 ? 0.0 : v;
}

/* -------------------------------------------------------------------------
 * Build per-column arrays for all N_POINTS scan points
 * ---------------------------------------------------------------------- */

typedef struct {
    double s2[N_POINTS];
    double ql[N_POINTS];
    double ei[N_POINTS];
    double ef[N_POINTS];
    double en[N_POINTS];
    double m1[N_POINTS], m2[N_POINTS];
    double monovf[N_POINTS], monohf[N_POINTS];
    double monotilt[N_POINTS], monotrans[N_POINTS];
    double a1[N_POINTS], a2[N_POINTS];
    double anavf[N_POINTS], anahf[N_POINTS];
    double anatilt[N_POINTS], anatrans[N_POINTS];
    double vs_left[N_POINTS], vs_right[N_POINTS];
    double ps_left[N_POINTS], ps_right[N_POINTS];
    double ps_top[N_POINTS],  ps_bottom[N_POINTS];
    double pa_left[N_POINTS], pa_right[N_POINTS];
    double pa_top[N_POINTS],  pa_bottom[N_POINTS];
    double col1[N_POINTS], col2[N_POINTS];
    double col3[N_POINTS], col4[N_POINTS];
    double temp1[N_POINTS], temp2[N_POINTS];
    double temp3[N_POINTS], temp4[N_POINTS];
    double cryo_he[N_POINTS], cryo_needle[N_POINTS];
    double counts[N_POINTS];
    double monitor[N_POINTS];
    /* PSD: flat [N_POINTS * 128 * 128] row-major */
    double *psd;
} scan_arrays_t;

static void fill_arrays(scan_arrays_t *a)
{
    for (int i = 0; i < N_POINTS; i++) {
        double t = (double)i / (double)(N_POINTS - 1);
        a->s2[i]  = S2_START + t * (S2_END - S2_START);
        a->ql[i]  = 1.4 + t * 0.2;
        a->ei[i]  = 14.6;
        a->ef[i]  = 14.6;
        a->en[i]  = 0.0;

        a->m1[i] = a->m2[i] = 45.0;
        a->monovf[i] = a->monohf[i] = 45.0;
        a->monotilt[i] = 45.0; a->monotrans[i] = 45.0;
        a->a1[i] = a->a2[i] = 60.0;
        a->anavf[i] = a->anahf[i] = 60.0;
        a->anatilt[i] = 60.0; a->anatrans[i] = 60.0;

        a->vs_left[i] = a->vs_right[i] = 10.0;
        a->ps_left[i] = a->ps_right[i] = 10.0;
        a->ps_top[i]  = a->ps_bottom[i]= 10.0;
        a->pa_left[i] = a->pa_right[i] = 10.0;
        a->pa_top[i]  = a->pa_bottom[i]= 10.0;
        a->col1[i] = a->col2[i] = 1.0;
        a->col3[i] = a->col4[i] = 1.0;
        a->temp1[i] = a->temp2[i] = 295.0;
        a->temp3[i] = a->temp4[i] = 295.0;
        a->cryo_he[i]     = 10.0;
        a->cryo_needle[i] = 7.0;
        a->monitor[i]     = 10000.0;

        double peak = 1000.0 * exp(-pow(a->s2[i] + 63.0, 2) / (2.0 * 0.25));
        a->counts[i] = peak + 20.0 * rng_normal();
        if (a->counts[i] < 0.0) a->counts[i] = 0.0;
    }

    /* PSD: Poisson noise around mean 100 */
    size_t total = (size_t)N_POINTS * TAS_PSD_FRAME_SIZE;
    a->psd = (double *)malloc(total * sizeof(double));
    for (size_t k = 0; k < total; k++)
        a->psd[k] = rng_poisson_approx(100.0);
}

/* -------------------------------------------------------------------------
 * Build metadata
 * ---------------------------------------------------------------------- */

static void fill_metadata(tas_metadata_t *m)
{
    static const char *users[] = {"Andrew Brown", "Alex Green", NULL};
    memset(m, 0, sizeof(*m));

    m->facility          = "ANSTO";
    m->source            = "OPAL Reactor";
    m->instrument_name   = "TAIPAN";
    m->experiment_id     = "exp1234";
    m->proposal_no       = "P21234";
    m->local_contact     = "Dr. Jane Doe";
    m->users             = users;
    m->mono_crystal      = "PG";
    m->ana_crystal       = "PG";
    m->sense             = "+-+";
    m->distance_vs_mono     = 3.5;
    m->distance_mono_sample = 2.0;
    m->distance_sample_ana  = 1.5;
    m->distance_ana_det     = 0.5;
    m->scan_no           = "1234567";
    m->title             = "s2 scan of Bragg (0 0 1.5)";
    m->command           = "scan s2 -65 -61 0.04 mon 10000";
    m->filename          = "TAIPAN_#1001235.tas.nxs.h5";
    m->sample_name       = "YBCO";
    m->sample_type       = "crystal";
    m->sample_mosaic     = 0.3;
    m->scanning_axis     = "s2";

    m->sample_v1[0] = 1.0;
    m->sample_v2[1] = 1.0;

    m->unit_cell[0] = 3.82;  m->unit_cell[1] = 3.82; m->unit_cell[2] = 11.68;
    m->unit_cell[3] = 90.0;  m->unit_cell[4] = 90.0; m->unit_cell[5] = 90.0;

    m->ub_matrix[0][0] = 0.1254; m->ub_matrix[0][1] = 0.0021;
    m->ub_matrix[1][0] = 0.0018; m->ub_matrix[1][1] = 0.1189;
    m->ub_matrix[2][2] = 0.0765;
}

/* =========================================================================
 * main
 * ========================================================================= */

int main(void)
{
    scan_arrays_t arr;
    fill_arrays(&arr);

    tas_metadata_t meta;
    fill_metadata(&meta);

    /* ── column table for the batch API ─────────────────────────────────── */
    tas_scalar_col_t cols[] = {
        {TAS_COL_S2,          arr.s2},
        {TAS_COL_QL,          arr.ql},
        {TAS_COL_EI,          arr.ei},
        {TAS_COL_EF,          arr.ef},
        {TAS_COL_EN,          arr.en},
        {TAS_COL_M1,          arr.m1},
        {TAS_COL_M2,          arr.m2},
        {TAS_COL_MONOVF,      arr.monovf},
        {TAS_COL_MONOHF,      arr.monohf},
        {TAS_COL_MONOTILT,    arr.monotilt},
        {TAS_COL_MONOTRANS,   arr.monotrans},
        {TAS_COL_A1,          arr.a1},
        {TAS_COL_A2,          arr.a2},
        {TAS_COL_ANAVF,       arr.anavf},
        {TAS_COL_ANAHF,       arr.anahf},
        {TAS_COL_ANATILT,     arr.anatilt},
        {TAS_COL_ANATRANS,    arr.anatrans},
        {TAS_COL_VS_LEFT,     arr.vs_left},
        {TAS_COL_VS_RIGHT,    arr.vs_right},
        {TAS_COL_PS_LEFT,     arr.ps_left},
        {TAS_COL_PS_RIGHT,    arr.ps_right},
        {TAS_COL_PS_TOP,      arr.ps_top},
        {TAS_COL_PS_BOTTOM,   arr.ps_bottom},
        {TAS_COL_PA_LEFT,     arr.pa_left},
        {TAS_COL_PA_RIGHT,    arr.pa_right},
        {TAS_COL_PA_TOP,      arr.pa_top},
        {TAS_COL_PA_BOTTOM,   arr.pa_bottom},
        {TAS_COL_COL1,        arr.col1},
        {TAS_COL_COL2,        arr.col2},
        {TAS_COL_COL3,        arr.col3},
        {TAS_COL_COL4,        arr.col4},
        {TAS_COL_TEMP1,       arr.temp1},
        {TAS_COL_TEMP2,       arr.temp2},
        {TAS_COL_TEMP3,       arr.temp3},
        {TAS_COL_TEMP4,       arr.temp4},
        {TAS_COL_CRYO_HE,     arr.cryo_he},
        {TAS_COL_CRYO_NEEDLE, arr.cryo_needle},
        {TAS_COL_COUNTS,      arr.counts},
        {TAS_COL_MONITOR,     arr.monitor},
    };
    size_t ncols = sizeof(cols) / sizeof(cols[0]);

    tas_scan_data_t scan = {
        .n_points  = N_POINTS,
        .columns   = cols,
        .n_columns = ncols,
        .psd_data  = arr.psd,
    };

    /* ── Test 1: batch write ─────────────────────────────────────────────── */
    const char *f1 = "test_batch_c.tas.nxs.h5";
    int rc = tas_save_hdf(f1, &scan, &meta, TAS_OVERWRITE);
    if (rc != TAS_OK) {
        fprintf(stderr, "ERROR (batch): %s\n", tas_last_error());
        free(arr.psd); return 1;
    }

    /* ── Test 2: point-by-point write ────────────────────────────────────── */
    const char *f2 = "test_pbp_c.tas.nxs.h5";
    tas_pbp_handle_t *h = tas_pbp_open(f2, N_POINTS, &meta);
    if (!h) {
        fprintf(stderr, "ERROR (pbp open): %s\n", tas_last_error());
        free(arr.psd); return 1;
    }

    for (int i = 0; i < N_POINTS; i++) {
        tas_point_t pt;
        memset(&pt, 0, sizeof(pt));

        pt.values[TAS_COL_S2]          = arr.s2[i];
        pt.values[TAS_COL_QL]          = arr.ql[i];
        pt.values[TAS_COL_EI]          = arr.ei[i];
        pt.values[TAS_COL_EF]          = arr.ef[i];
        pt.values[TAS_COL_EN]          = arr.en[i];
        pt.values[TAS_COL_M1]          = arr.m1[i];
        pt.values[TAS_COL_M2]          = arr.m2[i];
        pt.values[TAS_COL_MONOVF]      = arr.monovf[i];
        pt.values[TAS_COL_MONOHF]      = arr.monohf[i];
        pt.values[TAS_COL_MONOTILT]    = arr.monotilt[i];
        pt.values[TAS_COL_MONOTRANS]   = arr.monotrans[i];
        pt.values[TAS_COL_A1]          = arr.a1[i];
        pt.values[TAS_COL_A2]          = arr.a2[i];
        pt.values[TAS_COL_ANAVF]       = arr.anavf[i];
        pt.values[TAS_COL_ANAHF]       = arr.anahf[i];
        pt.values[TAS_COL_ANATILT]     = arr.anatilt[i];
        pt.values[TAS_COL_ANATRANS]    = arr.anatrans[i];
        pt.values[TAS_COL_VS_LEFT]     = arr.vs_left[i];
        pt.values[TAS_COL_VS_RIGHT]    = arr.vs_right[i];
        pt.values[TAS_COL_PS_LEFT]     = arr.ps_left[i];
        pt.values[TAS_COL_PS_RIGHT]    = arr.ps_right[i];
        pt.values[TAS_COL_PS_TOP]      = arr.ps_top[i];
        pt.values[TAS_COL_PS_BOTTOM]   = arr.ps_bottom[i];
        pt.values[TAS_COL_PA_LEFT]     = arr.pa_left[i];
        pt.values[TAS_COL_PA_RIGHT]    = arr.pa_right[i];
        pt.values[TAS_COL_PA_TOP]      = arr.pa_top[i];
        pt.values[TAS_COL_PA_BOTTOM]   = arr.pa_bottom[i];
        pt.values[TAS_COL_COL1]        = arr.col1[i];
        pt.values[TAS_COL_COL2]        = arr.col2[i];
        pt.values[TAS_COL_COL3]        = arr.col3[i];
        pt.values[TAS_COL_COL4]        = arr.col4[i];
        pt.values[TAS_COL_TEMP1]       = arr.temp1[i];
        pt.values[TAS_COL_TEMP2]       = arr.temp2[i];
        pt.values[TAS_COL_TEMP3]       = arr.temp3[i];
        pt.values[TAS_COL_TEMP4]       = arr.temp4[i];
        pt.values[TAS_COL_CRYO_HE]     = arr.cryo_he[i];
        pt.values[TAS_COL_CRYO_NEEDLE] = arr.cryo_needle[i];
        pt.values[TAS_COL_COUNTS]      = arr.counts[i];
        pt.values[TAS_COL_MONITOR]     = arr.monitor[i];

        /* PSD frame: pointer arithmetic into the flat array */
        pt.psd_frame = arr.psd + (size_t)i * TAS_PSD_FRAME_SIZE;

        rc = tas_pbp_write_point(h, &pt);
        if (rc < 0) {
            fprintf(stderr, "ERROR (pbp point %d): %s\n", i, tas_last_error());
            tas_pbp_close(h);
            free(arr.psd); return 1;
        }
    }

    rc = tas_pbp_close(h);
    if (rc != TAS_OK) {
        fprintf(stderr, "ERROR (pbp close): %s\n", tas_last_error());
        free(arr.psd); return 1;
    }

    free(arr.psd);

    printf("\nBoth files written successfully.\n");
    printf("  %s\n", f1);
    printf("  %s\n", f2);
    printf("\nVerify with:\n");
    printf("  h5ls -r %s\n", f1);
    printf("  python3 -c \"\n");
    printf("    import h5py, numpy as np\n");
    printf("    a=h5py.File('%s','r'); b=h5py.File('%s','r')\n", f1, f2);
    printf("    p='/entry/sample/s2/value'\n");
    printf("    print('s2 match:', np.allclose(a[p][()],b[p][()]))\n");
    printf("  \"\n");
    return 0;
}
