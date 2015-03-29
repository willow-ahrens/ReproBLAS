# Here's a good place to translate some of these settings into
# compilation flags/variables.  As an example a preprocessor macro for
# target endianess
ifeq ($(ENDIAN),big)
  CPPFLAGS += -DBIG_ENDIAN
else
  CPPFLAGS += -DLITTLE_ENDIAN
endif

# Detect C compiler in the following order if CC hasn't been set
ifeq ($(CC),cc)
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

# Detect MPI C compiler in the following order if MPICC hasn't been set
ifeq ($(MPICC),)
  ifeq ("$(shell which mpicc >/dev/null; echo $$?)", "0")
    MPICC := mpicc
  else
    MPICC := $(CC)
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
LDFLAGS += $(BLAS)
ifeq ($(strip $(BLAS)),-lblas)
  CPPFLAGS += -DBLAS
else ifeq ($(strip $(BLAS)),-lmkl_sequential)
  CPPFLAGS += -DCBLAS
else ifeq ($(strip $(BLAS)),-latlas)
  CPPFLAGS += -DCBLAS
else ifeq ($(strip $(BLAS)),-framework Accelerate)
  CPPFLAGS += -DCBLAS
endif

# Create cog compiler
COG = PYTHONPATH=$(TOP)/cog $(PYTHON) -m cogapp $(COGFLAGS) -D args=$(ARGS) -D params=$(PARAMS)

# Create cog compiler
PSEUDOCOG = $(TOP)/cog/pseudocog.sh
