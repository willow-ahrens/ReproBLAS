# Build mode e.g. debug, profile, release.  Build specific mode flags
# can be entered in $(MK)/build-$(BUILD_MODE).mk file e.g. for debug
# following seems to be a reasonable contents
#CFLAGS   += -ggdb
#CXXFLAGS += -ggdb
#CPPFLAGS += -DDEBUG
#LDFLAGS  += -ggdb
# If you don't plan on having different build modes then just comment
# below or set it to empty.
ifeq ($(BUILD_MODE),)
  BUILD_MODE := release
endif

# To have some per directory setting automatically propagated to all
# subdirs then uncomment below.  That way you can have all project
# compiled with "global" settings and easily switch flags for some
# subtree just by setting per directory settings at the top dir of the
# subtree.  You may of course overwrite inherited values and you can
# turn inheritance in some part by just clearing INHERIT_DIR_VARS_$(d)
# This is a global inheritance flag - you might want to turn it on only
# in some directory (just set INHERIT_DIR_VARS_$(d) there).
INHERIT_DIR_VARS := INCLUDES CPPFLAGS CFLAGS LDFLAGS MPICFLAGS MPILDFLAGS

# Default optimization flags.
ifeq ($(OPTFLAGS),)
  OPTFLAGS := -O3
endif

# Here's a good place to translate some of these settings into
# compilation flags/variables.  As an example a preprocessor macro for
# target endianess
ifeq ($(ENDIAN),big)
  CPPFLAGS += -DBIG_ENDIAN
else
  CPPFLAGS += -DLITTLE_ENDIAN
endif

# Detect C compiler in the following order if CC hasn't been set
ifeq ($(CC),)
  ifeq ("$(shell which gcc >/dev/null; echo $$?)", "0")
    CC := gcc
  else ifeq ("$(shell which icc >/dev/null; echo $$?)", "0")
    CC := icc
  else ifeq ("$(shell which pgcc >/dev/null; echo $$?)", "0")
    CC := pgcc
  else ifeq ("$(shell which craycc >/dev/null; echo $$?)", "0")
    CC := craycc
  else ifeq ("$(shell which clang >/dev/null; echo $$?)", "0")
    CC := clang
  endif
endif


#TODO do some better MPI flag detection.
# Detect MPI C compiler in the following order
ifeq ($(MPICC),)
  ifeq ("$(shell which mpicc >/dev/null; echo $$?)", "0")
    MPICC := mpicc
  else
    MPICC :=
  endif
endif

# Detect MPI C compiler flags in the following order if MPICFLAGS hasn't been set
ifeq ($(MPICFLAGS),)
  ifeq ($(MPICC), mpicc)
    MPICFLAGS := $(shell $(MPICC) --showme:compile)
  endif
endif

# Detect MPI C linker flags in the following order if MPILDFLAGS hasn't been set
ifeq ($(MPILDFLAGS),)
  ifeq ($(MPICC), mpicc)
    MPILDFLAGS := $(shell $(MPICC) --showme:link)
  endif
endif

# Detect python in the following order if PYTHON hasn't been set
ifeq ($(PYTHON),)
  ifeq ("$(shell which python3 >/dev/null; echo $$?)", "0")
    PYTHON := python3
  else ifeq ("$(shell which python >/dev/null; echo $$?)", "0")
    PYTHON := python
  endif
endif

# Use vectorization flags to determine compiler flags. MTARGET_ARCH is used to determine the target architectures vectorization settings.
ifeq ($(strip $(MMX))),true)
  CFLAGS += -mmmx
endif
ifeq ($(strip $(MMX)),false)
  CFLAGS += -mno-mmx
endif
ifeq ($(strip $(SSE)),true)
  CFLAGS += -msse
endif
ifeq ($(strip $(SSE)),false)
  CFLAGS += -mno-sse
endif
ifeq ($(strip $(SSE2)),true)
  CFLAGS += -msse2
endif
ifeq ($(strip $(SSE2)),false)
  CFLAGS += -mno-sse2
endif
ifeq ($(strip $(SSE3)),true)
  CFLAGS += -msse3
endif
ifeq ($(strip $(SSE3)),false)
  CFLAGS += -mno-sse3
endif
ifeq ($(strip $(SSE4_1)),true)
  CFLAGS += -msse4.1
endif
ifeq ($(strip $(SSE4_1)),false)
  CFLAGS += -mno-sse4.1
endif
ifeq ($(strip $(SSE4_2)),true)
  CFLAGS += -msse4.2
endif
ifeq ($(strip $(SSE4_2)),false)
  CFLAGS += -mno-sse4.2
endif
ifeq ($(strip $(AVX)),true)
  CFLAGS += -mavx
endif
ifeq ($(strip $(AVX)),false)
  CFLAGS += -mno-avx
endif
ifeq ($(strip $(AVX2)),true)
  CFLAGS += -mavx
endif
ifeq ($(strip $(AVX2)),false)
  CFLAGS += -mno-avx2
endif
ifeq ($(MMX)$(SSE)$(SSE1)$(SSE2)$(SSE3)$(SSE4_1)$(SSE4_2)$(AVX)$(AVX2),)
  ifeq ($(MTARGET_ARCH),)
    CFLAGS += -march=native
  else
    CFLAGS += -march=$(strip $(MTARGET_ARCH))
  endif
endif

# Determine BLAS link
ifeq ($(strip $(BLAS)),REF)
  LDFLAGS += -lblas
  CPPFLAGS += -DBLAS=1
else ifeq ($(strip $(BLAS)),MKL)
  LDFLAGS += -lmkl_sequential
  CPPFLAGS += -DCBLAS=1
else ifeq ($(strip $(BLAS)),ATLAS)
  LDFLAGS += -latlas
  CPPFLAGS += -DCBLAS=1
else ifeq ($(strip $(BLAS)),ACCELERATE)
  LDFLAGS += -framework Accelerate
  CPPFLAGS += -DCBLAS=1
endif

CALL_PYTHON = PYTHONPATH=$(TOP) $(PYTHON)

# Create cog compiler
COG = $(CALL_PYTHON) -m scripts.cogapp $(COGFLAGS)

# Create cog compiler
PSEUDOCOG = $(TOP)/cog/pseudocog.sh
