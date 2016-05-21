# select binary installation directory
BIN_DIR := /usr/local/bin

# select library installation directory
LIB_DIR := /usr/local/lib

# select header installation directory
INC_DIR := /usr/local/include

# select documentation installation directory
DOC_DIR := /usr/local/share/doc/reproBLAS2

# select compiler (comment all for auto)
CC = cc
#CC = gcc
#CC = icc
#CC = pgcc
#CC = craycc
#CC = clang

# add CFLAGS
CFLAGS += -Wall -std=c99

# add CPPFLAGS
CPPFLAGS +=

# add LDFLAGS
LDFLAGS +=

# select whether or not to build mpi (if BUILD_MPI is not defined or set to value other that "true" MPI will not be used and corresponding parts of ReproBLAS won't build)
BUILD_MPI = false

# select MPI compiler flags (comment all for auto)
#MPICFLAGS = $(shell mpicc --showme:compile)
#MPICFLAGS = $(shell mpicc --showme compile)
#MPICFLAGS = $(shell mpicc -compile_info)
#MPICFLAGS =

# select MPI linker flags (comment all for auto)
#MPILDFLAGS = $(shell mpicc --showme:link)
#MPILDFLAGS = $(shell mpicc --showme link)
#MPICFLAGS = $(shell mpicc -link_info)
#MPICFLAGS =

# select python (comment all for auto)
#PYTHON = python
#PYTHON = python3

# optionally disable vectorization (comment all for best available)
#SSE2 := false
#AVX := false

# select optimization flags (comment for auto)
OPTFLAGS := -O3

# select endianness (by default we are running on the same architecture we are
# building - if you're cross compiling then you should set this manually)
ENDIAN := $(shell perl -le 'print unpack(N,pack(L,0x01020304)) == 0x01020304 ? big : little')

# select whether or not to use of BLAS library (if BUILD_BLAS is not defined or set to value other that "true" an external BLAS library will not be used and corresponding parts of ReproBLAS won't build)
#BUILD_BLAS := true

# select BLAS library from preconfigured or custom options (if BLAS=true)

# reference BLAS
#LDFLAGS += -lblas
#CPPFLAGS += -DBLAS=1

# Intel MKL BLAS
#LDFLAGS += -mkl
#CPPFLAGS += -DCBLAS=1

# Intel MKL Sequential BLAS
#LDFLAGS += -mkl=sequential
#LDFLAGS += ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_core.a ${MKLROOT}/lib/libmkl_sequential.a -lpthread -lm
#CPPFLAGS += -DCBLAS=1

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
