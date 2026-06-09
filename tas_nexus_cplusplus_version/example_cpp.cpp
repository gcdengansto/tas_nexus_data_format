/**
 * @file example_cpp.cpp
 * @brief TAS NeXus C++ API self-test
 *
 * Exercises BatchWriter, PointByPointWriter, and Reader against the files
 * produced by the C writer self-test, and also writes new .h5 files via the
 * C++ API and cross-checks them.
 *
 * Build:
 *   make example_cpp
 *
 * Run:
 *   ./example_cpp
 */

#include "tas_nexus.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

/* ── Simple PRNG (same seed as the C example) ──────────────────────────── */
static unsigned long rng_state = 42;
static double rng_normal() {
    constexpr double PI = 3.14159265358979323846;
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    double u = ((rng_state >> 33) + 0.5) / (double)(1UL << 31);
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    double v  = ((rng_state >> 33) + 0.5) / (double)(1UL << 31);
    return std::sqrt(-2.0 * std::log(u)) * std::cos(2.0 * PI * v);
}
static double poisson_approx(double lam) {
    double v = lam + std::sqrt(lam) * rng_normal();
    return v < 0.0 ? 0.0 : v;
}

/* ── Check helpers ──────────────────────────────────────────────────────── */
static int g_pass = 0, g_total = 0;

template<typename T>
static void check(const char *label, T got, T expected, double tol = 0.0) {
    ++g_total;
    bool ok;
    if (tol > 0.0)
        ok = std::abs(static_cast<double>(got) - static_cast<double>(expected)) <= tol;
    else
        ok = (got == expected);
    if (ok) ++g_pass;
    std::cout << (ok ? "  \u2713  " : "  \u2717  ")
              << std::left << std::setw(28) << label
              << "  got=" << got
              << "  expected=" << expected << "\n";
}

static void check_str(const char *label, const std::string &got, const std::string &expected) {
    ++g_total;
    bool ok = (got == expected);
    if (ok) ++g_pass;
    std::cout << (ok ? "  \u2713  " : "  \u2717  ")
              << std::left << std::setw(28) << label
              << "  got=\"" << got << "\""
              << "  expected=\"" << expected << "\"\n";
}

/* ── Build metadata ─────────────────────────────────────────────────────── */
static tas::Metadata build_meta() {
    tas::Metadata m;
    m.facility         = "ANSTO";
    m.source           = "OPAL Reactor";
    m.instrument_name  = "TAIPAN";
    m.experiment_id    = "exp1234";
    m.proposal_no      = "P21234";
    m.local_contact    = "Dr. Jane Doe";
    m.users            = {"Andrew Brown", "Alex Green"};
    m.mono_crystal     = "PG";
    m.ana_crystal      = "PG";
    m.sense            = "+-+";
    m.distance_vs_mono      = 3.5;
    m.distance_mono_sample  = 2.0;
    m.distance_sample_ana   = 1.5;
    m.distance_ana_det      = 0.5;
    m.scan_no          = "1234567";
    m.title            = "s2 scan of Bragg (0 0 1.5)";
    m.command          = "scan s2 -65 -61 0.04 mon 10000";
    m.filename         = "TAIPAN_#1001235.tas.nxs.h5";
    m.sample_name      = "YBCO";
    m.sample_type      = "crystal";
    m.sample_mosaic    = 0.3;
    m.scanning_axis    = "s2";
    m.sample_v1        = {1.0, 0.0, 0.0};
    m.sample_v2        = {0.0, 1.0, 0.0};
    m.unit_cell        = {3.82, 3.82, 11.68, 90.0, 90.0, 90.0};
    m.ub_matrix = {
        0.1254, 0.0021, 0.0,
        0.0018, 0.1189, 0.0,
        0.0,    0.0,    0.0765
    };
    return m;
}

/* ── Build scan arrays ──────────────────────────────────────────────────── */
static constexpr int N = 101;
static constexpr double S2_START = -65.0, S2_END = -61.0;

struct Arrays {
    std::vector<double> s2, ql, ei, ef, en;
    std::vector<double> m1, m2, monovf, monohf, monotilt, monotrans;
    std::vector<double> a1, a2, anavf, anahf, anatilt, anatrans;
    std::vector<double> vs_left, vs_right;
    std::vector<double> ps_left, ps_right, ps_top, ps_bottom;
    std::vector<double> pa_left, pa_right, pa_top, pa_bottom;
    std::vector<double> col1, col2, col3, col4;
    std::vector<double> temp1, temp2, temp3, temp4;
    std::vector<double> cryo_he, cryo_needle;
    std::vector<double> counts, monitor;
    std::vector<double> psd;  /* flat [N][128][128] */
};

static Arrays make_arrays() {
    Arrays a;
    auto fill = [](std::vector<double> &v, int n, double val) {
        v.assign(n, val);
    };
    a.s2.resize(N);    a.ql.resize(N);
    a.ei.resize(N);    a.ef.resize(N);    a.en.resize(N);
    a.counts.resize(N); a.monitor.resize(N);

    for (int i = 0; i < N; i++) {
        double t = (double)i / (N - 1);
        a.s2[i]  = S2_START + t * (S2_END - S2_START);
        a.ql[i]  = 1.4 + t * 0.2;
        a.ei[i]  = 14.6;  a.ef[i] = 14.6;  a.en[i] = 0.0;
        a.monitor[i] = 10000.0;
        double peak = 1000.0 * std::exp(-std::pow(a.s2[i] + 63.0, 2) / (2.0 * 0.25));
        a.counts[i]  = std::max(0.0, peak + 20.0 * rng_normal());
    }

    fill(a.m1, N, 45.0);     fill(a.m2, N, 45.0);
    fill(a.monovf, N, 45.0); fill(a.monohf, N, 45.0);
    fill(a.monotilt, N, 45.0); fill(a.monotrans, N, 45.0);
    fill(a.a1, N, 60.0);     fill(a.a2, N, 60.0);
    fill(a.anavf, N, 60.0);  fill(a.anahf, N, 60.0);
    fill(a.anatilt, N, 60.0); fill(a.anatrans, N, 60.0);
    fill(a.vs_left, N, 10.0); fill(a.vs_right, N, 10.0);
    fill(a.ps_left, N, 10.0); fill(a.ps_right, N, 10.0);
    fill(a.ps_top, N, 10.0);  fill(a.ps_bottom, N, 10.0);
    fill(a.pa_left, N, 10.0); fill(a.pa_right, N, 10.0);
    fill(a.pa_top, N, 10.0);  fill(a.pa_bottom, N, 10.0);
    fill(a.col1, N, 1.0); fill(a.col2, N, 1.0);
    fill(a.col3, N, 1.0); fill(a.col4, N, 1.0);
    fill(a.temp1, N, 295.0); fill(a.temp2, N, 295.0);
    fill(a.temp3, N, 295.0); fill(a.temp4, N, 295.0);
    fill(a.cryo_he, N, 10.0); fill(a.cryo_needle, N, 7.0);

    size_t total = (size_t)N * tas::PSD_FRAME_SIZE;
    a.psd.resize(total);
    for (auto &v : a.psd) v = poisson_approx(100.0);

    return a;
}

static void populate_writer(tas::BatchWriter &w, const Arrays &a) {
    using C = tas::ColId;
    w.setColumn(C::S2,          a.s2);
    w.setColumn(C::QL,          a.ql);
    w.setColumn(C::EI,          a.ei);
    w.setColumn(C::EF,          a.ef);
    w.setColumn(C::EN,          a.en);
    w.setColumn(C::M1,          a.m1);
    w.setColumn(C::M2,          a.m2);
    w.setColumn(C::MONOVF,      a.monovf);
    w.setColumn(C::MONOHF,      a.monohf);
    w.setColumn(C::MONOTILT,    a.monotilt);
    w.setColumn(C::MONOTRANS,   a.monotrans);
    w.setColumn(C::A1,          a.a1);
    w.setColumn(C::A2,          a.a2);
    w.setColumn(C::ANAVF,       a.anavf);
    w.setColumn(C::ANAHF,       a.anahf);
    w.setColumn(C::ANATILT,     a.anatilt);
    w.setColumn(C::ANATRANS,    a.anatrans);
    w.setColumn(C::VS_LEFT,     a.vs_left);
    w.setColumn(C::VS_RIGHT,    a.vs_right);
    w.setColumn(C::PS_LEFT,     a.ps_left);
    w.setColumn(C::PS_RIGHT,    a.ps_right);
    w.setColumn(C::PS_TOP,      a.ps_top);
    w.setColumn(C::PS_BOTTOM,   a.ps_bottom);
    w.setColumn(C::PA_LEFT,     a.pa_left);
    w.setColumn(C::PA_RIGHT,    a.pa_right);
    w.setColumn(C::PA_TOP,      a.pa_top);
    w.setColumn(C::PA_BOTTOM,   a.pa_bottom);
    w.setColumn(C::COL1,        a.col1);
    w.setColumn(C::COL2,        a.col2);
    w.setColumn(C::COL3,        a.col3);
    w.setColumn(C::COL4,        a.col4);
    w.setColumn(C::TEMP1,       a.temp1);
    w.setColumn(C::TEMP2,       a.temp2);
    w.setColumn(C::TEMP3,       a.temp3);
    w.setColumn(C::TEMP4,       a.temp4);
    w.setColumn(C::CRYO_HE,     a.cryo_he);
    w.setColumn(C::CRYO_NEEDLE, a.cryo_needle);
    w.setColumn(C::COUNTS,      a.counts);
    w.setColumn(C::MONITOR,     a.monitor);
    w.setPsd(a.psd);
}

/* =========================================================================
 * Tests
 * ========================================================================= */

static void test_col_id_helpers() {
    std::cout << "\n── ColId helpers ────────────────────────────────────────\n";
    check_str("colName(S2)",    tas::colName(tas::ColId::S2),     "s2");
    check_str("colName(EI)",    tas::colName(tas::ColId::EI),     "ei");
    check_str("colName(COUNTS)",tas::colName(tas::ColId::COUNTS), "counts");

    auto id = tas::colId("s2");
    check("colId(s2) has value", id.has_value(), true);
    check("colId(s2)==S2",
          static_cast<int>(id.value()),
          static_cast<int>(tas::ColId::S2));

    check("colId(unknown) empty", tas::colId("__bad__").has_value(), false);
}

static void test_batch_write(const Arrays &a, const tas::Metadata &meta) {
    std::cout << "\n── BatchWriter ──────────────────────────────────────────\n";
    const char *path = "test_batch_cpp.tas.nxs.h5";

    tas::BatchWriter w(path, meta);
    populate_writer(w, a);
    w.write();   // throws on error

    std::cout << "  Written: " << path << "\n";
    check("nPoints", w.nPoints(), (size_t)N);
}

static void test_pbp_write(const Arrays &a, const tas::Metadata &meta) {
    std::cout << "\n── PointByPointWriter ───────────────────────────────────\n";
    const char *path = "test_pbp_cpp.tas.nxs.h5";

    {
        tas::PointByPointWriter pbp(path, N, meta);
        using C = tas::ColId;

        for (int i = 0; i < N; i++) {
            tas::ScanPoint pt;
            pt[C::S2]          = a.s2[i];
            pt[C::QL]          = a.ql[i];
            pt[C::EI]          = a.ei[i];
            pt[C::EF]          = a.ef[i];
            pt[C::EN]          = a.en[i];
            pt[C::M1]          = a.m1[i];
            pt[C::M2]          = a.m2[i];
            pt[C::MONOVF]      = a.monovf[i];
            pt[C::MONOHF]      = a.monohf[i];
            pt[C::MONOTILT]    = a.monotilt[i];
            pt[C::MONOTRANS]   = a.monotrans[i];
            pt[C::A1]          = a.a1[i];
            pt[C::A2]          = a.a2[i];
            pt[C::ANAVF]       = a.anavf[i];
            pt[C::ANAHF]       = a.anahf[i];
            pt[C::ANATILT]     = a.anatilt[i];
            pt[C::ANATRANS]    = a.anatrans[i];
            pt[C::VS_LEFT]     = a.vs_left[i];
            pt[C::VS_RIGHT]    = a.vs_right[i];
            pt[C::PS_LEFT]     = a.ps_left[i];
            pt[C::PS_RIGHT]    = a.ps_right[i];
            pt[C::PS_TOP]      = a.ps_top[i];
            pt[C::PS_BOTTOM]   = a.ps_bottom[i];
            pt[C::PA_LEFT]     = a.pa_left[i];
            pt[C::PA_RIGHT]    = a.pa_right[i];
            pt[C::PA_TOP]      = a.pa_top[i];
            pt[C::PA_BOTTOM]   = a.pa_bottom[i];
            pt[C::COL1]        = a.col1[i];
            pt[C::COL2]        = a.col2[i];
            pt[C::COL3]        = a.col3[i];
            pt[C::COL4]        = a.col4[i];
            pt[C::TEMP1]       = a.temp1[i];
            pt[C::TEMP2]       = a.temp2[i];
            pt[C::TEMP3]       = a.temp3[i];
            pt[C::TEMP4]       = a.temp4[i];
            pt[C::CRYO_HE]     = a.cryo_he[i];
            pt[C::CRYO_NEEDLE] = a.cryo_needle[i];
            pt[C::COUNTS]      = a.counts[i];
            pt[C::MONITOR]     = a.monitor[i];

            /* PSD frame */
            std::vector<double> frame(
                a.psd.begin() + (size_t)i * tas::PSD_FRAME_SIZE,
                a.psd.begin() + (size_t)(i + 1) * tas::PSD_FRAME_SIZE);
            pt.setPsd(frame);

            pbp.writePoint(pt);
        }
        check("written()", pbp.written(), (size_t)N);
        check("isOpen()",  pbp.isOpen(),  true);
    } // ← RAII close here

    std::cout << "  Written: " << path << " (RAII close)\n";
}

static void test_reader(const Arrays &a) {
    std::cout << "\n── Reader — load batch file ─────────────────────────────\n";

    auto scan = tas::Reader::load("test_batch_cpp.tas.nxs.h5");

    /* Summary */
    scan.printSummary();

    /* Round-trip checks */
    std::cout << "── Round-trip value checks ──────────────────────────────\n";

    const auto &s2  = scan.column(tas::ColId::S2);
    const auto &ei  = scan.column(tas::ColId::EI);
    const auto &mon = scan.column(tas::ColId::MONITOR);

    check("nPoints",      scan.nPoints(), (size_t)N);
    check("s2[0]",        s2[0],   S2_START,  1e-9);
    check("s2[N-1]",      s2[N-1], S2_END,    1e-9);
    double ei_mean = std::accumulate(ei.begin(), ei.end(), 0.0) / ei.size();
    check("ei mean",      ei_mean, 14.6, 1e-9);
    check("monitor[0]",   mon[0],  10000.0,   1e-9);

    const auto &inf = scan.info();
    check_str("facility",      inf.facility,          "ANSTO");
    check_str("instrument",    inf.instrument_name,    "TAIPAN");
    check_str("sample_name",   inf.sample_name,        "YBCO");
    check_str("scanning_axis", inf.scanning_axis,       "s2");
    check_str("sense",         inf.sense,              "+-+");
    check_str("mono_crystal",  inf.mono_crystal,        "PG");
    check("sample_mosaic",   inf.sample_mosaic,         0.3,  1e-9);
    check("dist_vs_mono",    inf.distance_vs_mono,      3.5,  1e-9);
    check("dist_mono_samp",  inf.distance_mono_sample,  2.0,  1e-9);
    check("unit_cell[0]",    inf.unit_cell[0],           3.82, 1e-6);
    check("unit_cell[2]",    inf.unit_cell[2],          11.68, 1e-6);
    check("ub[0][0]",        inf.ub_matrix[0],          0.1254, 1e-6);
    check("ub[8]",           inf.ub_matrix[8],          0.0765, 1e-6);
    check("v1[0]",           inf.sample_v1[0],           1.0,  1e-9);
    check("v2[1]",           inf.sample_v2[1],           1.0,  1e-9);
    check("n_users",         (int)inf.users.size(),       2);
    check_str("user[0]",     inf.users[0], "Andrew Brown");

    /* PSD checks */
    check("hasPsd()",    scan.hasPsd(),  true);
    auto f0 = scan.psdFrame(0);
    double fsum = 0.0;
    for (int r = 0; r < tas::PSD_NY; r++)
        for (int c = 0; c < tas::PSD_NX; c++)
            fsum += f0(r, c);
    check("psd frame[0] sum > 0", fsum > 0.0, true);

    /* value() accessor */
    check("value(S2,0)",    scan.value(tas::ColId::S2, 0), S2_START, 1e-9);
    check("value(EI,0)",    scan.value(tas::ColId::EI, 0), 14.6,     1e-9);

    /* Data identity: C++ batch == C++ pbp */
    std::cout << "\n── Cross-file identity check (batch vs pbp) ─────────────\n";
    auto scan_pbp = tas::Reader::load("test_pbp_cpp.tas.nxs.h5");
    const auto &s2b = scan_pbp.column(tas::ColId::S2);
    bool s2_match = std::equal(s2.begin(), s2.end(), s2b.begin(),
                               [](double x, double y){ return std::abs(x-y) < 1e-12; });
    check("s2 batch==pbp",  s2_match, true);

    const auto &cnt_b = scan.column(tas::ColId::COUNTS);
    const auto &cnt_p = scan_pbp.column(tas::ColId::COUNTS);
    bool cnt_match = std::equal(cnt_b.begin(), cnt_b.end(), cnt_p.begin(),
                                [](double x, double y){ return std::abs(x-y) < 1e-12; });
    check("counts batch==pbp", cnt_match, true);

    /* Text export */
    scan.exportText("test_batch_cpp_exported.txt");
    std::cout << "  Text export written: test_batch_cpp_exported.txt\n";
}

static void test_error_handling() {
    std::cout << "\n── Error handling ───────────────────────────────────────\n";

    /* Reading a non-existent file must throw */
    bool threw = false;
    try { tas::Reader::load("__does_not_exist__.h5"); }
    catch (const tas::Error &e) { threw = true; (void)e; }
    check("load missing file throws", threw, true);

    /* BatchWriter with mismatched column sizes must throw */
    threw = false;
    try {
        tas::BatchWriter w("dummy.h5", tas::Metadata{});
        w.setColumn(tas::ColId::S2,     std::vector<double>(10, 0.0));
        w.setColumn(tas::ColId::COUNTS, std::vector<double>(5,  0.0)); // wrong size
    }
    catch (const tas::Error &e) { threw = true; (void)e; }
    check("mismatched column sizes throws", threw, true);

    /* ScanPoint bad PSD size must throw */
    threw = false;
    try {
        tas::ScanPoint pt;
        pt.setPsd(std::vector<double>(100, 0.0)); // too short
    }
    catch (const tas::Error &e) { threw = true; (void)e; }
    check("bad PSD frame size throws", threw, true);

    /* TAS_SKIP_PSD: psd should be absent */
    auto scan = tas::Reader::load("test_batch_cpp.tas.nxs.h5", /*loadPsd=*/false);
    check("SKIP_PSD → hasPsd()=false", scan.hasPsd(), false);
}

/* =========================================================================
 * main
 * ========================================================================= */

int main() {
    std::cout << "=======================================================\n";
    std::cout << "TAS NeXus C++ API self-test\n";
    std::cout << "=======================================================\n";

    try {
        auto meta = build_meta();
        auto arrs = make_arrays();

        test_col_id_helpers();
        test_batch_write(arrs, meta);
        test_pbp_write(arrs, meta);
        test_reader(arrs);
        test_error_handling();

    } catch (const tas::Error &e) {
        std::cerr << "\nFATAL tas::Error: " << e.what()
                  << "  (code=" << e.code() << ")\n";
        return 1;
    } catch (const std::exception &e) {
        std::cerr << "\nFATAL std::exception: " << e.what() << "\n";
        return 1;
    }

    std::cout << "\n=======================================================\n";
    if (g_pass == g_total)
        std::cout << "All " << g_total << " checks passed.\n";
    else
        std::cout << g_pass << "/" << g_total << " checks passed — SOME FAILED.\n";
    std::cout << "=======================================================\n";
    return (g_pass == g_total) ? 0 : 1;
}
