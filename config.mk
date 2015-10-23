# select compiler (comment all for auto)
CC = cc
#CC = gcc
#CC = icc
#CC = pgcc
#CC = craycc
#CC = clang

# add CFLAGS
CFLAGS += -Wall

# add CPPFLAGS
CPPFLAGS +=

# add LDFLAGS
LDFLAGS +=

# select whether or not to build mpi (comment not to build mpi. Corresponding
# parts of ReproBLAS won't build)
#BUILD_MPI = true

# select MPI compiler flags (comment all for auto)
#MPICFLAGS = $(shell mpicc --showme compile)
#MPICFLAGS = $(shell mpicc -compile_info)
#MPICFLAGS =

# select MPI linker flags (comment all for auto)
#MPILDFLAGS = $(shell mpicc --showme link)
#MPICFLAGS = $(shell mpicc -link_info)
#MPICFLAGS =

# select python (comment all for auto)
#PYTHON = python
#PYTHON = python3

# select vectorization (comment all for auto)
#MMX := true
#MMX := false
#SSE := true
#SSE := false
#SSE2 := true
#SSE2 := false
#SSE3 := true
#SSE3 := false
#SSE4_1 := true
#SSE4_1 := false
#SSE4_2 := true
#SSE4_2 := false
#AVX := true
#AVX := false
#AVX2 := true
#AVX2 := false

# select optimization flags (comment for auto)
OPTFLAGS := -O3 -funroll-loops

# select endianness (by default we are running on the same architecture we are
# building - if you're cross compiling then you should set this manually)
ENDIAN := $(shell perl -le 'print unpack(N,pack(L,0x01020304)) == 0x01020304 ? big : little')

# select whether or not to use of BLAS library (if BLAS is not defined or set to value other that "true" an external BLAS library will not be used and corresponding parts of ReproBLAS won't build)
BLAS := true

# select BLAS library from preconfigured or custom options (if BLAS=true)

# reference BLAS
#LDFLAGS += -lblas
#CPPFLAGS += -DBLAS=1

# Intel MKL Sequential BLAS
LDFLAGS += ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_core.a ${MKLROOT}/lib/libmkl_sequential.a -lpthread -lm
CPPFLAGS += -DCBLAS=1

# Atlas BLAS
#LDFLAGS += -latlas
#CPPFLAGS += -DCBLAS=1

# Mac OS Accelerate BLAS
#LDFLAGS += -framework Accelerate
#CPPFLAGS += -DCBLAS=1

# custom BLAS LDFLAGS
#LDFLAGS += -lblas
# custom BLAS interface (reference fortran or cblas)
#CPPFLAGS += -DBLAS=1
#CPPFLAGS += -DCBLAS=1

# select build mode (comment all for auto)
#BUILD_MODE := release
#BUILD_MODE := debug
#BUILD_MODE := profile

#HOST_ARCH := profile

# select arguments file (contains all of the values for tuning parameters in the
# library. One can give this parameter on the command line as well. If you have
# run the autotuner, you can point this at the output file to use the tuned
# arguments)
ARGS = $(TOP)/src/tuned_args.json

# select parameters file (contains all tuning parameters in the library. This 
# parameter likely does not need modification)
PARAMS = $(TOP)/src/params.json

# select verbosity (if VERBOSE is not defined or set to value other than "true"
# you'll see just indication of what is being compiled or run (without details
# about options))
#VERBOSE := true
