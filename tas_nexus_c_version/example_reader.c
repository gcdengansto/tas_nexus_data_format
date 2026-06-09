/**
 * @file example_reader.c
 * @brief TAS NeXus reader — C self-test mirroring tas_nexus_reader.py __main__
 *
 * Reads the file produced by the writer self-test (test_batch_c.tas.nxs.h5),
 * performs all round-trip value checks from the Python reference, prints a
 * summary, and writes an export text file.
 *
 * Build (from the project directory):
 *
 *   make example_reader
 *
 * Run:
 *
 *   ./example_reader [path/to/scan.h5]
 *   # default: test_batch_c.tas.nxs.h5
 */

#define _POSIX_C_SOURCE 200809L
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

#include "tas_nexus_reader.h"
#include "tas_nexus_writer.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Simple round-trip check helpers ────────────────────────────────────── */

static int checks_passed = 0;
static int checks_total  = 0;

static void check_dbl(const char *label, double got, double expected, double tol)
{
    checks_total++;
    bool ok = fabs(got - expected) <= tol;
    if (ok) checks_passed++;
    printf("  %s  %-22s  got=%-12.6g  expected=%-12.6g\n",
           ok ? "✓" : "✗", label, got, expected);
}

static void check_str(const char *label, const char *got, const char *expected)
{
    checks_total++;
    bool ok = (got && strcmp(got, expected) == 0);
    if (ok) checks_passed++;
    printf("  %s  %-22s  got=\"%s\"  expected=\"%s\"\n",
           ok ? "✓" : "✗", label, got ? got : "(null)", expected);
}

static void check_int(const char *label, int got, int expected)
{
    checks_total++;
    bool ok = (got == expected);
    if (ok) checks_passed++;
    printf("  %s  %-22s  got=%-6d  expected=%d\n",
           ok ? "✓" : "✗", label, got, expected);
}

/* =========================================================================
 * main
 * ========================================================================= */

int main(int argc, char *argv[])
{
    const char *test_file = (argc > 1) ? argv[1] : "test_batch_c.tas.nxs.h5";

    printf("\n%s\n", "=======================================================");
    printf("Reading : %s\n", test_file);
    printf("%s\n",   "=======================================================");

    /* ── Load ──────────────────────────────────────────────────────────── */
    tas_scan_t *scan = tas_read_hdf(test_file, TAS_LOAD_PSD);
    if (!scan) {
        fprintf(stderr, "\nERROR: %s\n", tas_last_error());
        fprintf(stderr, "Run the writer self-test first:\n"
                        "  make test\n");
        return 1;
    }

    /* ── Summary ───────────────────────────────────────────────────────── */
    tas_print_summary(scan);

    /* ── Round-trip value checks ────────────────────────────────────────── */
    printf("── Round-trip value checks ───────────────────────────────\n");

    size_t n = scan->n_points;
    const double *s2      = tas_col_ptr(scan, TAS_COL_S2);
    const double *ei_col  = tas_col_ptr(scan, TAS_COL_EI);
    const double *mon_col = tas_col_ptr(scan, TAS_COL_MONITOR);

    check_dbl("s2 first",       s2      ? s2[0]     : -999.0, -65.0, 1e-6);
    check_dbl("s2 last",        s2      ? s2[n-1]   : -999.0, -61.0, 1e-6);
    check_dbl("ei mean",
              (ei_col ? (({ double sum=0; for(size_t i=0;i<n;i++) sum+=ei_col[i]; sum; })/n) : -1.0),
              14.6, 1e-6);
    check_dbl("monitor[0]",     mon_col ? mon_col[0] : -999.0, 10000.0, 1e-6);
    check_int("n_points",       scan->info.num_points, 101);
    check_dbl("dist_vs_mono",   scan->info.distance_vs_mono,     3.5, 1e-9);
    check_dbl("dist_mono_samp", scan->info.distance_mono_sample, 2.0, 1e-9);
    check_dbl("sample_mosaic",  scan->info.sample_mosaic,        0.3, 1e-9);
    check_str("scanning_axis",  scan->info.scanning_axis, "s2");
    check_str("sense",          scan->info.sense,        "+-+");
    check_str("mono_crystal",   scan->info.mono_crystal, "PG");
    check_str("facility",       scan->info.facility,     "ANSTO");
    check_str("sample_name",    scan->info.sample_name,  "YBCO");

    /* unit_cell[0] == 3.82 */
    check_dbl("unit_cell[0]",   scan->info.unit_cell[0], 3.82, 1e-6);
    check_dbl("unit_cell[2]",   scan->info.unit_cell[2], 11.68, 1e-6);
    /* UB matrix corner */
    check_dbl("ub[0][0]",       scan->info.ub_matrix[0][0], 0.1254, 1e-6);
    check_dbl("ub[2][2]",       scan->info.ub_matrix[2][2], 0.0765, 1e-6);
    /* sample vectors */
    check_dbl("v1[0]",          scan->info.sample_v1[0], 1.0, 1e-9);
    check_dbl("v2[1]",          scan->info.sample_v2[1], 1.0, 1e-9);
    /* PSD */
    check_int("psd_present",    (int)scan->info.psd_present, 1);

    /* ── Column accessor checks ─────────────────────────────────────────── */
    check_dbl("tas_col_value s2[0]",
              tas_col_value(scan, TAS_COL_S2, 0), -65.0, 1e-6);
    check_dbl("tas_col_value ei[0]",
              tas_col_value(scan, TAS_COL_EI, 0), 14.6, 1e-6);

    /* PSD frame accessor */
    const double *frame0 = tas_psd_frame(scan, 0);
    if (frame0) {
        double fsum = 0.0;
        for (int k = 0; k < TAS_PSD_FRAME_SIZE; k++) fsum += frame0[k];
        /* Poisson mean=100, 128*128 pixels → expect ~1 638 400 ± noise */
        check_dbl("psd frame[0] sum > 0", fsum > 0 ? 1.0 : 0.0, 1.0, 1e-9);
    } else {
        check_int("psd_frame(0) non-null", 0, 1);
    }

    printf("\n%s\n\n",
           (checks_passed == checks_total)
           ? "All round-trip checks passed."
           : "*** SOME CHECKS FAILED ***");

    /* ── Export text ───────────────────────────────────────────────────── */
    const char *txt = "test_batch_c_exported.txt";
    int rc = tas_export_text(scan, txt);
    if (rc == TAS_OK)
        printf("Text export written: %s\n", txt);
    else
        fprintf(stderr, "Export failed: %s\n", tas_last_error());

    /* ── Scalar-only load check ─────────────────────────────────────────── */
    printf("\nTesting TAS_SKIP_PSD ...\n");
    tas_scan_t *scan2 = tas_read_hdf(test_file, TAS_SKIP_PSD);
    if (scan2) {
        printf("  n_points = %zu  psd = %s\n",
               scan2->n_points, scan2->psd ? "non-NULL (unexpected)" : "NULL (correct)");
        tas_scan_free(scan2);
    }

    tas_scan_free(scan);
    return (checks_passed == checks_total) ? 0 : 1;
}
