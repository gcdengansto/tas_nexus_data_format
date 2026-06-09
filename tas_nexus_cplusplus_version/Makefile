# Makefile — TAS NeXus HDF5 writer + reader  (C and C++)
# ─────────────────────────────────────────────────────────────────
# Targets
#   make              → all: C static+shared libs + all example binaries
#   make static       → C static library only
#   make shared       → C shared library only
#   make example      → C writer self-test binary
#   make example_reader → C reader self-test binary
#   make example_cpp  → C++ self-test binary
#   make test         → writer C test (produces .h5 files)
#   make test_reader  → reader C test (needs .h5 from above)
#   make test_cpp     → C++ test (writer+reader in one)
#   make test_all     → all three tests in sequence
#   make clean        → remove build artefacts
#   make install      → install headers + libs (PREFIX=/usr/local)
#
# Requirements:  libhdf5-dev  (apt-get install libhdf5-dev)
# ─────────────────────────────────────────────────────────────────

CC      ?= gcc
CXX     ?= g++
CFLAGS  := -O2 -Wall -Wextra -std=c11   -fPIC -D_POSIX_C_SOURCE=200809L
CXXFLAGS:= -O2 -Wall -Wextra -std=c++17 -fPIC -D_POSIX_C_SOURCE=200809L

# ── HDF5 flags ───────────────────────────────────────────────────
H5_SHOW   := $(shell h5cc -show 2>/dev/null | sed 's/^[^ ]* //')
H5_CFLAGS := $(filter -I%,$(H5_SHOW))
H5_LDADD  := $(filter-out -I% -D%,$(H5_SHOW))

ifeq ($(strip $(H5_LDADD)),)
  H5_CFLAGS := $(shell pkg-config --cflags hdf5 2>/dev/null)
  H5_LDADD  := $(shell pkg-config --libs   hdf5 2>/dev/null)
endif

# ── Source files ─────────────────────────────────────────────────
LIB_NAME := tas_nexus
C_SRCS   := tas_nexus_writer.c tas_nexus_reader.c
C_OBJS   := $(C_SRCS:.c=.o)
C_HDRS   := tas_nexus_writer.h tas_nexus_reader.h
CXX_HDR  := tas_nexus.hpp

STATIC_LIB := lib$(LIB_NAME).a
SHARED_LIB := lib$(LIB_NAME).so

WRITER_BIN  := example_tas
READER_BIN  := example_reader
CPP_BIN     := example_cpp

# ─────────────────────────────────────────────────────────────────
.PHONY: all static shared example example_reader example_cpp \
        test test_reader test_cpp test_all clean install

all: static shared example example_reader example_cpp

# ── C object files ────────────────────────────────────────────────
%.o: %.c $(C_HDRS)
	$(CC) $(CFLAGS) $(H5_CFLAGS) -c $< -o $@

# ── Static library ────────────────────────────────────────────────
static: $(STATIC_LIB)

$(STATIC_LIB): $(C_OBJS)
	ar rcs $@ $^
	@echo "Built: $@"

# ── Shared library ────────────────────────────────────────────────
shared: $(SHARED_LIB)

$(SHARED_LIB): $(C_OBJS)
	$(CC) -shared -o $@ $^ $(H5_LDADD) -lm
	@echo "Built: $@"

# ── C writer example ──────────────────────────────────────────────
example: $(WRITER_BIN)

$(WRITER_BIN): example.c $(STATIC_LIB) $(C_HDRS)
	$(CC) $(CFLAGS) $(H5_CFLAGS) $< $(STATIC_LIB) $(H5_LDADD) -lm -o $@
	@echo "Built: $@"

# ── C reader example ──────────────────────────────────────────────
example_reader: $(READER_BIN)

$(READER_BIN): example_reader.c $(STATIC_LIB) $(C_HDRS)
	$(CC) $(CFLAGS) $(H5_CFLAGS) $< $(STATIC_LIB) $(H5_LDADD) -lm -o $@
	@echo "Built: $@"

# ── C++ example ───────────────────────────────────────────────────
# The C++ binary links the C++ translation unit together with the
# C library objects so the linker sees COL_META and tas_last_error.
example_cpp: $(CPP_BIN)

$(CPP_BIN): example_cpp.cpp $(STATIC_LIB) $(C_HDRS) $(CXX_HDR)
	$(CXX) $(CXXFLAGS) $(H5_CFLAGS) $< $(STATIC_LIB) $(H5_LDADD) -lm -o $@
	@echo "Built: $@"

# ── Tests ─────────────────────────────────────────────────────────
test: $(WRITER_BIN)
	@echo "─── C writer self-test ──────────────────────────────────"
	./$(WRITER_BIN)
	@echo "─── Done ────────────────────────────────────────────────"

test_reader: $(READER_BIN)
	@echo "─── C reader self-test ──────────────────────────────────"
	./$(READER_BIN)
	@echo "─── Done ────────────────────────────────────────────────"

test_cpp: $(CPP_BIN)
	@echo "─── C++ self-test ───────────────────────────────────────"
	./$(CPP_BIN)
	@echo "─── Done ────────────────────────────────────────────────"

test_all: test test_reader test_cpp

# ── Install ──────────────────────────────────────────────────────
PREFIX ?= /usr/local

install: $(STATIC_LIB) $(SHARED_LIB)
	install -d $(DESTDIR)$(PREFIX)/include
	install -d $(DESTDIR)$(PREFIX)/lib
	install -m 644 $(C_HDRS) $(CXX_HDR) $(DESTDIR)$(PREFIX)/include/
	install -m 644 $(STATIC_LIB)         $(DESTDIR)$(PREFIX)/lib/
	install -m 755 $(SHARED_LIB)         $(DESTDIR)$(PREFIX)/lib/
	ldconfig $(DESTDIR)$(PREFIX)/lib || true

# ── Clean ────────────────────────────────────────────────────────
clean:
	$(RM) $(C_OBJS) $(STATIC_LIB) $(SHARED_LIB) \
	      $(WRITER_BIN) $(READER_BIN) $(CPP_BIN) \
	      test_batch_c.tas.nxs.h5   test_pbp_c.tas.nxs.h5 \
	      test_batch_cpp.tas.nxs.h5 test_pbp_cpp.tas.nxs.h5 \
	      test_batch_c_exported.txt test_batch_cpp_exported.txt
