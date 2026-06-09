/**
 * @file tas_nexus.hpp
 * @brief TAS NeXus HDF5 writer and reader — C++17 public API
 *
 * A modern C++ wrapper around the C API (tas_nexus_writer / tas_nexus_reader).
 * All HDF5 I/O is performed by the proven C layer; this header provides
 * idiomatic C++ types, RAII resource management, and std::optional / exception
 * error handling.
 *
 * Requires C++17.  Link together with tas_nexus_writer.c and tas_nexus_reader.c
 * (and libhdf5).
 *
 * ---------------------------------------------------------------------------
 * NAMESPACE  tas
 * ---------------------------------------------------------------------------
 *
 *  Writer classes
 *  ──────────────
 *   tas::Metadata          — all experiment / instrument metadata fields
 *   tas::ScanPoint         — values for one scan point (motor + detector)
 *   tas::BatchWriter       — write an entire scan in one call
 *   tas::PointByPointWriter— RAII writer for live point-by-point acquisition
 *
 *  Reader classes
 *  ──────────────
 *   tas::ScanData          — complete in-memory scan (scalars + PSD + metadata)
 *   tas::Reader            — factory: Reader::load() → ScanData
 *
 *  Shared types
 *  ────────────
 *   tas::ColId             — scoped enum wrapping tas_col_id_t
 *   tas::PsdFrame          — lightweight non-owning view of one PSD frame
 *   tas::Error             — exception thrown on I/O or argument errors
 *
 * ---------------------------------------------------------------------------
 * QUICK START — WRITING (BATCH)
 * ---------------------------------------------------------------------------
 *
 *   #include "tas_nexus.hpp"
 *
 *   tas::Metadata meta;
 *   meta.facility        = "ANSTO";
 *   meta.instrument_name = "TAIPAN";
 *   meta.sample_name     = "YBCO";
 *
 *   std::vector<double> s2(101), counts(101);
 *   // … fill arrays …
 *
 *   tas::BatchWriter writer("scan.h5", meta);
 *   writer.setColumn(tas::ColId::S2,     s2);
 *   writer.setColumn(tas::ColId::COUNTS, counts);
 *   writer.write();               // throws tas::Error on failure
 *
 * ---------------------------------------------------------------------------
 * QUICK START — WRITING (POINT-BY-POINT)
 * ---------------------------------------------------------------------------
 *
 *   tas::PointByPointWriter pbp("scan.h5", 101, meta);
 *
 *   for (int i = 0; i < 101; i++) {
 *       tas::ScanPoint pt;
 *       pt[tas::ColId::S2]     = s2_pos[i];
 *       pt[tas::ColId::COUNTS] = counts[i];
 *       pbp.writePoint(pt);       // throws tas::Error on failure
 *   }
 *   // file is closed and end_time stamped when pbp goes out of scope
 *
 * ---------------------------------------------------------------------------
 * QUICK START — READING
 * ---------------------------------------------------------------------------
 *
 *   auto scan = tas::Reader::load("scan.h5");   // loads PSD by default
 *
 *   std::cout << "Points : " << scan.nPoints()         << "\n";
 *   std::cout << "Facility: " << scan.info().facility  << "\n";
 *
 *   const auto& s2 = scan.column(tas::ColId::S2);   // std::vector<double>&
 *   double peak = *std::max_element(s2.begin(), s2.end());
 *
 *   if (scan.hasPsd()) {
 *       tas::PsdFrame f = scan.psdFrame(0);   // non-owning view, 128×128
 *       double pixel = f(row, col);
 *   }
 *
 *   scan.exportText("output.txt");
 *   scan.printSummary();
 *
 * ---------------------------------------------------------------------------
 * BUILDING
 * ---------------------------------------------------------------------------
 *
 *   h5c++ -O2 -std=c++17 my_app.cpp tas_nexus_writer.c tas_nexus_reader.c \
 *         -lstdc++ -lm -o my_app
 *
 * Or via CMake (see CMakeLists.txt).
 *
 * ---------------------------------------------------------------------------
 * ERROR HANDLING
 * ---------------------------------------------------------------------------
 *
 * All public methods throw tas::Error (derived from std::runtime_error) on
 * failure.  No error codes or output-parameter checks needed.
 *
 * ---------------------------------------------------------------------------
 * THREAD SAFETY
 * ---------------------------------------------------------------------------
 *
 * ScanData is read-only after construction — safe to share across threads.
 * PointByPointWriter is NOT thread-safe; use one instance per thread.
 */

#pragma once
#ifndef TAS_NEXUS_HPP
#define TAS_NEXUS_HPP

#include "tas_nexus_writer.h"   /* C enum, constants, C API */
#include "tas_nexus_reader.h"   /* C reader API              */

#include <algorithm>
#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace tas {

/* =========================================================================
 * Error type
 * ========================================================================= */

/**
 * @brief Exception thrown by all tas:: functions on failure.
 *
 * Carries the same message that would be returned by tas_last_error() in
 * the C API, plus an optional integer error code (TAS_ERR_*).
 */
class Error : public std::runtime_error {
public:
    explicit Error(const std::string &msg, int code = TAS_ERR_GENERAL)
        : std::runtime_error(msg), code_(code) {}

    /** TAS_ERR_* code from tas_nexus_writer.h */
    int code() const noexcept { return code_; }

private:
    int code_;
};

/* =========================================================================
 * ColId — scoped enum wrapping tas_col_id_t
 * ========================================================================= */

/**
 * @brief Strongly-typed column identifier.
 *
 * Maps 1:1 onto tas_col_id_t.  Use tas::colId() / tas::colName() to
 * convert between this enum and C identifiers or string names.
 */
enum class ColId : int {
    S1 = TAS_COL_S1, S2 = TAS_COL_S2, SGU = TAS_COL_SGU, SGL = TAS_COL_SGL,
    STU = TAS_COL_STU, STL = TAS_COL_STL,
    QH = TAS_COL_QH, QK = TAS_COL_QK, QL = TAS_COL_QL,
    EI = TAS_COL_EI, EF = TAS_COL_EF, EN = TAS_COL_EN,
    M1 = TAS_COL_M1, M2 = TAS_COL_M2,
    MONOVF = TAS_COL_MONOVF, MONOHF = TAS_COL_MONOHF,
    MONOTILT = TAS_COL_MONOTILT, MONOTRANS = TAS_COL_MONOTRANS,
    A1 = TAS_COL_A1, A2 = TAS_COL_A2,
    ANAVF = TAS_COL_ANAVF, ANAHF = TAS_COL_ANAHF,
    ANATILT = TAS_COL_ANATILT, ANATRANS = TAS_COL_ANATRANS,
    VS_LEFT = TAS_COL_VS_LEFT, VS_RIGHT = TAS_COL_VS_RIGHT,
    PS_LEFT = TAS_COL_PS_LEFT, PS_RIGHT = TAS_COL_PS_RIGHT,
    PS_TOP  = TAS_COL_PS_TOP,  PS_BOTTOM = TAS_COL_PS_BOTTOM,
    PA_LEFT = TAS_COL_PA_LEFT, PA_RIGHT = TAS_COL_PA_RIGHT,
    PA_TOP  = TAS_COL_PA_TOP,  PA_BOTTOM = TAS_COL_PA_BOTTOM,
    COL1 = TAS_COL_COL1, COL2 = TAS_COL_COL2,
    COL3 = TAS_COL_COL3, COL4 = TAS_COL_COL4,
    COLL_ALPHA1 = TAS_COL_COLL_ALPHA1, COLL_ALPHA2 = TAS_COL_COLL_ALPHA2,
    COLL_ALPHA3 = TAS_COL_COLL_ALPHA3, COLL_ALPHA4 = TAS_COL_COLL_ALPHA4,
    TEMP1 = TAS_COL_TEMP1, TEMP2 = TAS_COL_TEMP2,
    TEMP3 = TAS_COL_TEMP3, TEMP4 = TAS_COL_TEMP4,
    MFIELD = TAS_COL_MFIELD, EFIELD = TAS_COL_EFIELD,
    PRESSURE = TAS_COL_PRESSURE,
    CRYO_HE = TAS_COL_CRYO_HE, CRYO_NEEDLE = TAS_COL_CRYO_NEEDLE,
    COUNTS = TAS_COL_COUNTS, MONITOR = TAS_COL_MONITOR,
    N_COLS = TAS_N_COLS   /**< Sentinel */
};

/** Convert ColId to C tas_col_id_t */
inline tas_col_id_t toCId(ColId c) noexcept {
    return static_cast<tas_col_id_t>(static_cast<int>(c));
}

/** Convert tas_col_id_t to ColId */
inline ColId fromCId(tas_col_id_t c) noexcept {
    return static_cast<ColId>(static_cast<int>(c));
}

/** Return the canonical name string for a column (e.g. "s2"). */
inline std::string colName(ColId c) {
    const char *n = tas_col_name(toCId(c));
    return n ? std::string(n) : std::string();
}

/** Look up a ColId by name.  Returns std::nullopt if not found. */
inline std::optional<ColId> colId(const std::string &name) {
    int rc = tas_col_from_name(name.c_str());
    if (rc < 0) return std::nullopt;
    return fromCId(static_cast<tas_col_id_t>(rc));
}

/* =========================================================================
 * PSD constants (mirrors C macros)
 * ========================================================================= */

static constexpr int PSD_NX         = TAS_PSD_NX;          ///< 128
static constexpr int PSD_NY         = TAS_PSD_NY;          ///< 128
static constexpr int PSD_FRAME_SIZE = TAS_PSD_FRAME_SIZE;  ///< 128*128

/* =========================================================================
 * PsdFrame — non-owning 2-D view of one PSD frame
 * ========================================================================= */

/**
 * @brief Lightweight, non-owning 2-D view of one 128×128 PSD frame.
 *
 * The pointed-to data is owned by the ScanData that produced this view.
 * Do not use the view after the parent ScanData has been destroyed.
 */
class PsdFrame {
public:
    /** Construct a view over @p data (must point to PSD_NY * PSD_NX doubles). */
    explicit PsdFrame(const double *data) noexcept : data_(data) {}

    /** Access pixel at (row, col), 0-based. */
    double operator()(int row, int col) const noexcept {
        return data_[row * PSD_NX + col];
    }

    /** Raw pointer to the PSD_FRAME_SIZE doubles, row-major. */
    const double *data() const noexcept { return data_; }

private:
    const double *data_;
};

/* =========================================================================
 * Metadata
 * ========================================================================= */

/**
 * @brief All experiment / instrument metadata for one TAS scan.
 *
 * Every field has a sensible default (empty string / 0.0 / zero array).
 * Unset fields are written with the same defaults as the C API.
 */
struct Metadata {
    /* Facility / instrument */
    std::string facility;
    std::string source;
    std::string instrument_name;
    std::string software_version;
    std::string tas_nexus_version;

    /* Proposal */
    std::string experiment_id;
    std::string proposal_no;
    std::string local_contact;

    /* Scan */
    std::string scan_no;
    std::string title;
    std::string command;
    std::string filename;
    std::string scanning_axis;
    std::string start_time;   ///< ISO-8601; empty → current time
    std::string end_time;     ///< ISO-8601; empty → current time at close

    /* Sample */
    std::string sample_name;
    std::string sample_type {"crystal"};
    double      sample_mosaic {0.0};             ///< minutes of arc
    std::array<double, 3> sample_v1 {1., 0., 0.};
    std::array<double, 3> sample_v2 {0., 1., 0.};
    std::array<double, 6> unit_cell {3.82, 3.82, 11.68, 90., 90., 90.};
    std::array<double, 9> ub_matrix {};          ///< row-major 3×3; zeros → default

    /* Optics */
    std::string mono_crystal {"PG"};
    std::string ana_crystal  {"PG"};
    std::string sense        {"+-+"};

    /* Distances (metres) */
    double distance_vs_mono      {0.0};
    double distance_mono_sample  {2.0};
    double distance_sample_ana   {1.5};
    double distance_ana_det      {0.5};

    /* Users */
    std::vector<std::string> users;

    // ── Internal conversion helpers ─────────────────────────────────────

    /**
     * @brief Build a C tas_metadata_t that borrows from this struct's strings.
     *
     * The returned struct is only valid while *this and the c_users vector
     * remain alive.  Use immediately; do not store.
     */
    tas_metadata_t toCMeta(std::vector<const char *> &c_users) const {
        tas_metadata_t m{};

        auto s = [](const std::string &str) -> const char * {
            return str.empty() ? nullptr : str.c_str();
        };

        m.facility          = s(facility);
        m.source            = s(source);
        m.instrument_name   = s(instrument_name);
        m.software_version  = s(software_version);
        m.tas_nexus_version = s(tas_nexus_version);
        m.experiment_id     = s(experiment_id);
        m.proposal_no       = s(proposal_no);
        m.local_contact     = s(local_contact);
        m.scan_no           = s(scan_no);
        m.title             = s(title);
        m.command           = s(command);
        m.filename          = s(filename);
        m.scanning_axis     = s(scanning_axis);
        m.start_time        = s(start_time);
        m.end_time          = s(end_time);
        m.sample_name       = s(sample_name);
        m.sample_type       = s(sample_type);
        m.sample_mosaic     = sample_mosaic;
        m.mono_crystal      = s(mono_crystal);
        m.ana_crystal       = s(ana_crystal);
        m.sense             = s(sense);
        m.distance_vs_mono      = distance_vs_mono;
        m.distance_mono_sample  = distance_mono_sample;
        m.distance_sample_ana   = distance_sample_ana;
        m.distance_ana_det      = distance_ana_det;

        for (int i = 0; i < 3; i++) {
            m.sample_v1[i] = sample_v1[i];
            m.sample_v2[i] = sample_v2[i];
        }
        for (int i = 0; i < 6; i++) m.unit_cell[i] = unit_cell[i];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                m.ub_matrix[r][c] = ub_matrix[r * 3 + c];

        c_users.clear();
        for (const auto &u : users) c_users.push_back(u.c_str());
        c_users.push_back(nullptr);
        m.users = c_users.data();

        return m;
    }

    /** Build a Metadata from a C tas_scan_info_t (used by Reader). */
    static Metadata fromCInfo(const tas_scan_info_t &info) {
        Metadata m;
        auto str = [](const char *s) -> std::string {
            return s ? std::string(s) : std::string();
        };
        m.facility          = str(info.facility);
        m.source            = str(info.source);
        m.instrument_name   = str(info.instrument_name);
        m.experiment_id     = str(info.experiment_id);
        m.proposal_no       = str(info.proposal_no);
        m.local_contact     = str(info.local_contact);
        m.scan_no           = str(info.scan_no);
        m.title             = str(info.title);
        m.command           = str(info.command);
        m.filename          = str(info.filename);
        m.scanning_axis     = str(info.scanning_axis);
        m.start_time        = str(info.start_time);
        m.end_time          = str(info.end_time);
        m.sample_name       = str(info.sample_name);
        m.sample_type       = str(info.sample_type);
        m.sample_mosaic     = info.sample_mosaic;
        m.mono_crystal      = str(info.mono_crystal);
        m.ana_crystal       = str(info.ana_crystal);
        m.sense             = str(info.sense);
        m.distance_vs_mono      = info.distance_vs_mono;
        m.distance_mono_sample  = info.distance_mono_sample;
        m.distance_sample_ana   = info.distance_sample_ana;
        m.distance_ana_det      = info.distance_ana_det;
        for (int i = 0; i < 3; i++) {
            m.sample_v1[i] = info.sample_v1[i];
            m.sample_v2[i] = info.sample_v2[i];
        }
        for (int i = 0; i < 6; i++) m.unit_cell[i] = info.unit_cell[i];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                m.ub_matrix[r * 3 + c] = info.ub_matrix[r][c];
        for (int i = 0; i < info.n_users; i++)
            if (info.users[i]) m.users.push_back(info.users[i]);
        return m;
    }
};

/* =========================================================================
 * ScanPoint — one scan point for point-by-point writing
 * ========================================================================= */

/**
 * @brief Values for a single scan point (all motors + detectors).
 *
 * Unset channels default to 0.0.  The optional psd field holds one full
 * PSD_NY × PSD_NX frame (16 384 doubles).
 */
class ScanPoint {
public:
    ScanPoint() { values_.fill(0.0); }

    /** Set / get motor or detector value by ColId. */
    double &operator[](ColId c) { return values_[static_cast<int>(c)]; }
    double  operator[](ColId c) const { return values_[static_cast<int>(c)]; }

    /**
     * @brief Attach a PSD frame.
     *
     * @p frame must contain exactly PSD_FRAME_SIZE doubles (row-major).
     * The data is copied into the ScanPoint's internal storage.
     */
    void setPsd(const std::vector<double> &frame) {
        if (static_cast<int>(frame.size()) != PSD_FRAME_SIZE)
            throw Error("ScanPoint::setPsd: frame must have PSD_FRAME_SIZE elements",
                        TAS_ERR_ARG);
        psd_ = frame;
    }

    /** @overload */
    void setPsd(std::vector<double> &&frame) {
        if (static_cast<int>(frame.size()) != PSD_FRAME_SIZE)
            throw Error("ScanPoint::setPsd: frame must have PSD_FRAME_SIZE elements",
                        TAS_ERR_ARG);
        psd_ = std::move(frame);
    }

    /** Clear any attached PSD frame. */
    void clearPsd() noexcept { psd_.clear(); }

    /** True if a PSD frame has been attached. */
    bool hasPsd() const noexcept { return !psd_.empty(); }

    /** Raw pointer to PSD data (nullptr if no frame attached). */
    const double *psdData() const noexcept {
        return psd_.empty() ? nullptr : psd_.data();
    }

    /** Raw pointer to all TAS_N_COLS values (needed by BatchWriter). */
    const double *rawValues() const noexcept { return values_.data(); }

private:
    std::array<double, TAS_N_COLS> values_;
    std::vector<double> psd_;
};

/* =========================================================================
 * BatchWriter
 * ========================================================================= */

/**
 * @brief Write an entire scan to a NeXus HDF5 file in one call.
 *
 * Usage:
 *   tas::BatchWriter w("out.h5", meta);
 *   w.setColumn(tas::ColId::S2,     s2_vec);
 *   w.setColumn(tas::ColId::COUNTS, cnt_vec);
 *   w.setPsd(psd_flat);   // optional
 *   w.write();            // throws on error
 *
 * All column vectors must have the same length (n_points).  n_points is
 * inferred from the first column set; subsequent columns are checked.
 */
class BatchWriter {
public:
    /**
     * @param filepath   Output file path.
     * @param meta       Experiment metadata.
     * @param overwrite  If false, throws Error if the file already exists.
     */
    BatchWriter(std::string filepath, Metadata meta, bool overwrite = true)
        : filepath_(std::move(filepath))
        , meta_(std::move(meta))
        , overwrite_(overwrite)
    {}

    /**
     * @brief Set the data array for a scalar column.
     *
     * @p data is copied into internal storage.  Repeated calls for the same
     * column overwrite the previous data.
     *
     * @throws Error if data.size() is inconsistent with previously set columns.
     */
    BatchWriter &setColumn(ColId col, const std::vector<double> &data) {
        checkSize(data.size());
        columns_[static_cast<int>(col)] = data;
        return *this;
    }

    /** @overload — move from a temporary vector. */
    BatchWriter &setColumn(ColId col, std::vector<double> &&data) {
        checkSize(data.size());
        columns_[static_cast<int>(col)] = std::move(data);
        return *this;
    }

    /**
     * @brief Set the PSD dataset from a flat row-major array.
     *
     * @p psd must have exactly n_points * PSD_FRAME_SIZE elements.
     * n_points must already be determined (call setColumn first).
     *
     * @throws Error on size mismatch.
     */
    BatchWriter &setPsd(const std::vector<double> &psd) {
        if (n_points_ == 0)
            throw Error("BatchWriter::setPsd: set at least one column before setPsd",
                        TAS_ERR_ARG);
        if (psd.size() != n_points_ * static_cast<size_t>(PSD_FRAME_SIZE))
            throw Error("BatchWriter::setPsd: psd.size() must equal n_points * PSD_FRAME_SIZE",
                        TAS_ERR_ARG);
        psd_ = psd;
        return *this;
    }

    /** @overload */
    BatchWriter &setPsd(std::vector<double> &&psd) {
        if (n_points_ == 0)
            throw Error("BatchWriter::setPsd: set at least one column before setPsd",
                        TAS_ERR_ARG);
        if (psd.size() != n_points_ * static_cast<size_t>(PSD_FRAME_SIZE))
            throw Error("BatchWriter::setPsd: psd.size() must equal n_points * PSD_FRAME_SIZE",
                        TAS_ERR_ARG);
        psd_ = std::move(psd);
        return *this;
    }

    /**
     * @brief Build and write the HDF5 file.
     *
     * @throws Error on any I/O or argument failure.
     */
    void write() {
        if (n_points_ == 0)
            throw Error("BatchWriter::write: no column data has been set", TAS_ERR_ARG);

        /* Build C column list */
        std::vector<tas_scalar_col_t> c_cols;
        for (int i = 0; i < TAS_N_COLS; i++) {
            if (columns_[i].empty()) continue;
            tas_scalar_col_t col{};
            col.column = static_cast<tas_col_id_t>(i);
            col.data   = columns_[i].data();
            c_cols.push_back(col);
        }

        tas_scan_data_t scan{};
        scan.n_points  = n_points_;
        scan.columns   = c_cols.data();
        scan.n_columns = c_cols.size();
        scan.psd_data  = psd_.empty() ? nullptr : psd_.data();

        std::vector<const char *> c_users;
        tas_metadata_t c_meta = meta_.toCMeta(c_users);

        int flags = overwrite_ ? TAS_OVERWRITE : TAS_NO_OVERWRITE;
        int rc    = tas_save_hdf(filepath_.c_str(), &scan, &c_meta, flags);
        if (rc != TAS_OK)
            throw Error(std::string(tas_last_error()), rc);
    }

    /** Number of scan points determined so far (0 if no columns set). */
    size_t nPoints() const noexcept { return n_points_; }

private:
    void checkSize(size_t sz) {
        if (sz == 0)
            throw Error("BatchWriter: column data must not be empty", TAS_ERR_ARG);
        if (n_points_ == 0) {
            n_points_ = sz;
        } else if (sz != n_points_) {
            throw Error("BatchWriter: all columns must have the same length", TAS_ERR_ARG);
        }
    }

    std::string  filepath_;
    Metadata     meta_;
    bool         overwrite_;
    size_t       n_points_ {0};

    std::array<std::vector<double>, TAS_N_COLS> columns_;
    std::vector<double> psd_;
};

/* =========================================================================
 * PointByPointWriter
 * ========================================================================= */

/**
 * @brief RAII writer for live point-by-point acquisition.
 *
 * The HDF5 file is created on construction and closed (with end_time
 * stamped) when the object is destroyed or close() is explicitly called.
 *
 * Usage:
 *   {
 *       tas::PointByPointWriter pbp("out.h5", 101, meta);
 *       for (auto &pt : points) pbp.writePoint(pt);
 *   }  // ← file closed here
 */
class PointByPointWriter {
public:
    /**
     * @param filepath   Output file path (always created fresh).
     * @param n_points   Total number of points that will be written.
     * @param meta       Experiment metadata.
     *
     * @throws Error if the file cannot be created.
     */
    PointByPointWriter(std::string filepath,
                       size_t      n_points,
                       Metadata    meta)
        : filepath_(std::move(filepath))
        , meta_(std::move(meta))
        , n_points_(n_points)
    {
        std::vector<const char *> c_users;
        tas_metadata_t c_meta = meta_.toCMeta(c_users);
        handle_ = tas_pbp_open(filepath_.c_str(), n_points_, &c_meta);
        if (!handle_)
            throw Error(std::string(tas_last_error()), TAS_ERR_HDF5);
    }

    /** Destructor — closes the file gracefully (swallows errors). */
    ~PointByPointWriter() noexcept {
        if (handle_) {
            tas_pbp_close(handle_);
            handle_ = nullptr;
        }
    }

    /* Non-copyable, movable */
    PointByPointWriter(const PointByPointWriter &)            = delete;
    PointByPointWriter &operator=(const PointByPointWriter &) = delete;

    PointByPointWriter(PointByPointWriter &&o) noexcept
        : filepath_(std::move(o.filepath_))
        , meta_(std::move(o.meta_))
        , n_points_(o.n_points_)
        , written_(o.written_)
        , handle_(o.handle_)
    {
        o.handle_ = nullptr;
    }

    PointByPointWriter &operator=(PointByPointWriter &&o) noexcept {
        if (this != &o) {
            if (handle_) tas_pbp_close(handle_);
            filepath_ = std::move(o.filepath_);
            meta_     = std::move(o.meta_);
            n_points_ = o.n_points_;
            written_  = o.written_;
            handle_   = o.handle_;
            o.handle_ = nullptr;
        }
        return *this;
    }

    /**
     * @brief Write one scan point to the file.
     *
     * @return 1-based index of the point just written.
     * @throws Error on failure or if already closed.
     */
    int writePoint(const ScanPoint &pt) {
        if (!handle_)
            throw Error("PointByPointWriter: file is already closed", TAS_ERR_ARG);

        tas_point_t c_pt{};
        std::copy(pt.rawValues(), pt.rawValues() + TAS_N_COLS, c_pt.values);
        c_pt.psd_frame = pt.psdData();

        int rc = tas_pbp_write_point(handle_, &c_pt);
        if (rc < 0)
            throw Error(std::string(tas_last_error()), rc);
        ++written_;
        return rc;
    }

    /**
     * @brief Explicitly flush, stamp end_time, and close the file.
     *
     * After this call the object is in a "closed" state.  Calling
     * writePoint() afterwards throws.  The destructor becomes a no-op.
     *
     * @throws Error on failure.
     */
    void close() {
        if (!handle_) return;
        int rc = tas_pbp_close(handle_);
        handle_ = nullptr;
        if (rc != TAS_OK)
            throw Error(std::string(tas_last_error()), rc);
    }

    /** True if the file is still open. */
    bool isOpen()   const noexcept { return handle_ != nullptr; }
    /** Number of points written so far. */
    size_t written() const noexcept { return written_; }
    /** Total pre-allocated points. */
    size_t nPoints() const noexcept { return n_points_; }

private:
    std::string          filepath_;
    Metadata             meta_;
    size_t               n_points_  {0};
    size_t               written_   {0};
    tas_pbp_handle_t    *handle_    {nullptr};
};

/* =========================================================================
 * ScanData — complete in-memory scan returned by Reader
 * ========================================================================= */

/**
 * @brief Complete in-memory representation of one TAS NeXus scan.
 *
 * Produced by Reader::load().  All data is owned by this object.
 * ScanData is move-only (it owns large data vectors).
 */
class ScanData {
public:
    /* --- Construction (internal — use Reader::load) -------------------- */

    /** Build from a C tas_scan_t pointer (takes ownership; frees on construction). */
    explicit ScanData(tas_scan_t *raw) {
        if (!raw) throw Error("ScanData: null C scan pointer", TAS_ERR_ARG);

        n_points_ = raw->n_points;
        info_     = Metadata::fromCInfo(raw->info);

        /* Copy scalar channels */
        for (int c = 0; c < TAS_N_COLS; c++) {
            if (raw->data[c]) {
                cols_[c].assign(raw->data[c], raw->data[c] + n_points_);
            } else {
                cols_[c].assign(n_points_, 0.0);
            }
        }

        /* Copy PSD */
        if (raw->psd) {
            size_t total = n_points_ * static_cast<size_t>(PSD_FRAME_SIZE);
            psd_.assign(raw->psd, raw->psd + total);
            psd_present_ = raw->info.psd_present;
        }

        tas_scan_free(raw);   /* release C memory */
    }

    /* Move-only */
    ScanData(const ScanData &)            = delete;
    ScanData &operator=(const ScanData &) = delete;
    ScanData(ScanData &&)                 = default;
    ScanData &operator=(ScanData &&)      = default;

    /* --- Accessors ------------------------------------------------------ */

    /** Number of scan points. */
    size_t nPoints() const noexcept { return n_points_; }

    /** All metadata. */
    const Metadata &info() const noexcept { return info_; }

    /**
     * @brief Return the full data vector for a scalar column.
     *
     * Always returns a vector of length nPoints() (all-zeros if the
     * channel was absent in the file).
     *
     * @throws Error if col is out of range.
     */
    const std::vector<double> &column(ColId col) const {
        int idx = static_cast<int>(col);
        if (idx < 0 || idx >= TAS_N_COLS)
            throw Error("ScanData::column: col out of range", TAS_ERR_ARG);
        return cols_[idx];
    }

    /** Return the value for column @p col at scan point @p idx. */
    double value(ColId col, size_t idx) const {
        const auto &v = column(col);
        if (idx >= n_points_)
            throw Error("ScanData::value: idx out of range", TAS_ERR_ARG);
        return v[idx];
    }

    /** True if a non-zero PSD array was loaded. */
    bool hasPsd() const noexcept { return psd_present_; }

    /** Raw flat PSD array [n_points][PSD_NY][PSD_NX].  Empty if no PSD. */
    const std::vector<double> &psdRaw() const noexcept { return psd_; }

    /**
     * @brief Return a non-owning view of PSD frame @p idx.
     *
     * @throws Error if PSD is absent or idx is out of range.
     */
    PsdFrame psdFrame(size_t idx) const {
        if (!psd_present_ || psd_.empty())
            throw Error("ScanData::psdFrame: no PSD data loaded", TAS_ERR_ARG);
        if (idx >= n_points_)
            throw Error("ScanData::psdFrame: idx out of range", TAS_ERR_ARG);
        return PsdFrame(psd_.data() + idx * static_cast<size_t>(PSD_FRAME_SIZE));
    }

    /**
     * @brief Export the scan to a formatted text file.
     *
     * Calls the C tas_export_text() underneath; reproduces the Python format.
     *
     * @throws Error on failure.
     */
    void exportText(const std::string &filepath) const {
        /* Rebuild a temporary C tas_scan_t for the export function */
        tas_scan_t tmp{};
        tmp.n_points = n_points_;

        std::vector<std::vector<double>> col_bufs(TAS_N_COLS);
        for (int c = 0; c < TAS_N_COLS; c++) {
            col_bufs[c] = cols_[c];
            tmp.data[c] = col_bufs[c].data();
        }
        std::vector<double> psd_buf = psd_;
        tmp.psd = psd_buf.empty() ? nullptr : psd_buf.data();

        /* Build C info */
        tas_scan_info_t &ci = tmp.info;
        auto fill = [](char *&dest, const std::string &src) {
            dest = const_cast<char *>(src.c_str());
        };
        fill(ci.start_time,       info_.start_time);
        fill(ci.end_time,         info_.end_time);
        fill(ci.title,            info_.title);
        fill(ci.command,          info_.command);
        fill(ci.filename,         info_.filename);
        fill(ci.scanning_axis,    info_.scanning_axis);
        fill(ci.scan_no,          info_.scan_no);
        fill(ci.facility,         info_.facility);
        fill(ci.source,           info_.source);
        fill(ci.instrument_name,  info_.instrument_name);
        fill(ci.local_contact,    info_.local_contact);
        fill(ci.experiment_id,    info_.experiment_id);
        fill(ci.proposal_no,      info_.proposal_no);
        fill(ci.sample_name,      info_.sample_name);
        fill(ci.sample_type,      info_.sample_type);
        fill(ci.mono_crystal,     info_.mono_crystal);
        fill(ci.ana_crystal,      info_.ana_crystal);
        fill(ci.sense,            info_.sense);

        ci.sample_mosaic        = info_.sample_mosaic;
        ci.distance_vs_mono     = info_.distance_vs_mono;
        ci.distance_mono_sample = info_.distance_mono_sample;
        ci.distance_sample_ana  = info_.distance_sample_ana;
        ci.distance_ana_det     = info_.distance_ana_det;
        for (int i = 0; i < 3; i++) {
            ci.sample_v1[i] = info_.sample_v1[i];
            ci.sample_v2[i] = info_.sample_v2[i];
        }
        for (int i = 0; i < 6; i++) ci.unit_cell[i] = info_.unit_cell[i];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                ci.ub_matrix[r][c] = info_.ub_matrix[r * 3 + c];

        std::vector<char *> u_ptrs;
        for (const auto &u : info_.users)
            u_ptrs.push_back(const_cast<char *>(u.c_str()));
        ci.n_users = static_cast<int>(u_ptrs.size());
        for (int i = 0; i < ci.n_users && i < TAS_MAX_USERS; i++)
            ci.users[i] = u_ptrs[i];

        ci.num_points   = static_cast<int>(n_points_);
        ci.psd_present  = psd_present_;

        int rc = tas_export_text(&tmp, filepath.c_str());
        if (rc != TAS_OK)
            throw Error(std::string(tas_last_error()), rc);
    }

    /**
     * @brief Print a diagnostic summary to stdout.
     *
     * Rebuilds a temporary C tas_scan_t and calls tas_print_summary().
     */
    void printSummary() const {
        /* Build a minimal tas_scan_t pointing into our vectors */
        tas_scan_t tmp{};
        tmp.n_points = n_points_;
        for (int c = 0; c < TAS_N_COLS; c++)
            tmp.data[c] = const_cast<double *>(cols_[c].data());
        tmp.psd = psd_.empty() ? nullptr : const_cast<double *>(psd_.data());

        /* Populate info — borrow string pointers */
        tas_scan_info_t &ci = tmp.info;
        auto fill = [](char *&dest, const std::string &src) {
            dest = const_cast<char *>(src.c_str());
        };
        fill(ci.start_time,      info_.start_time);
        fill(ci.end_time,        info_.end_time);
        fill(ci.title,           info_.title);
        fill(ci.command,         info_.command);
        fill(ci.filename,        info_.filename);
        fill(ci.scanning_axis,   info_.scanning_axis);
        fill(ci.scan_no,         info_.scan_no);
        fill(ci.facility,        info_.facility);
        fill(ci.source,          info_.source);
        fill(ci.instrument_name, info_.instrument_name);
        fill(ci.local_contact,   info_.local_contact);
        fill(ci.experiment_id,   info_.experiment_id);
        fill(ci.proposal_no,     info_.proposal_no);
        fill(ci.sample_name,     info_.sample_name);
        fill(ci.sample_type,     info_.sample_type);
        fill(ci.mono_crystal,    info_.mono_crystal);
        fill(ci.ana_crystal,     info_.ana_crystal);
        fill(ci.sense,           info_.sense);

        ci.sample_mosaic        = info_.sample_mosaic;
        ci.distance_vs_mono     = info_.distance_vs_mono;
        ci.distance_mono_sample = info_.distance_mono_sample;
        ci.distance_sample_ana  = info_.distance_sample_ana;
        ci.distance_ana_det     = info_.distance_ana_det;
        for (int i = 0; i < 3; i++) {
            ci.sample_v1[i] = info_.sample_v1[i];
            ci.sample_v2[i] = info_.sample_v2[i];
        }
        for (int i = 0; i < 6; i++) ci.unit_cell[i] = info_.unit_cell[i];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                ci.ub_matrix[r][c] = info_.ub_matrix[r * 3 + c];

        std::vector<char *> u_ptrs;
        for (const auto &u : info_.users)
            u_ptrs.push_back(const_cast<char *>(u.c_str()));
        ci.n_users = static_cast<int>(u_ptrs.size());
        for (int i = 0; i < ci.n_users && i < TAS_MAX_USERS; i++)
            ci.users[i] = u_ptrs[i];

        ci.num_points  = static_cast<int>(n_points_);
        ci.psd_present = psd_present_;

        tas_print_summary(&tmp);
    }

private:
    size_t   n_points_    {0};
    bool     psd_present_ {false};
    Metadata info_;
    std::array<std::vector<double>, TAS_N_COLS> cols_;
    std::vector<double> psd_;
};

/* =========================================================================
 * Reader
 * ========================================================================= */

/**
 * @brief Factory class for loading TAS NeXus HDF5 files.
 *
 * All methods are static; Reader has no state.
 */
class Reader {
public:
    /**
     * @brief Load a TAS NeXus HDF5 file into memory.
     *
     * @param filepath  Path to the .tas.nxs.h5 file.
     * @param loadPsd   If true (default), load the PSD array.
     *
     * @return ScanData owning all loaded data.
     * @throws Error on I/O failure or missing data.
     */
    static ScanData load(const std::string &filepath, bool loadPsd = true) {
        int flags = loadPsd ? TAS_LOAD_PSD : TAS_SKIP_PSD;
        tas_scan_t *raw = tas_read_hdf(filepath.c_str(), flags);
        if (!raw)
            throw Error(std::string(tas_last_error()), TAS_ERR_HDF5);
        return ScanData(raw);   /* constructor calls tas_scan_free(raw) */
    }

    Reader() = delete;  /* pure static utility class */
};

} /* namespace tas */

#endif /* TAS_NEXUS_HPP */
