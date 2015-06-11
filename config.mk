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
#OPTFLAGS := -O2

# select endianness (by default we are running on the same architecture we are
# building - if you're cross compiling then you should set this manually)
ENDIAN := $(shell perl -le 'print unpack(N,pack(L,0x01020304)) == 0x01020304 ? big : little')

# select BLAS library from preconfigured or CUSTOM options. (comment all to not
# use external BLAS library. Corresponding parts of ReproBLAS won't build)
#BLAS := REF
#BLAS := ATLAS
#BLAS := MKL
#BLAS := ACCELERATE
#BLAS := CUSTOM

# select CUSTOM BLAS LDFLAGS (if BLAS == CUSTOM)
#LDFLAGS += -lblas
# select CUSTOM BLAS interface (reference fortran or cblas) (if BLAS == CUSTOM)
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
ARGS = $(TOP)/src/default_args.json

# select parameters file (contains all tuning parameters in the library. This 
# parameter likely does not need modification)
PARAMS = $(TOP)/src/params.json

# select verbosity (if VERBOSE is not defined or set to value other than "true"
# you'll see just indication of what is being compiled or run (without details
# about options))
#VERBOSE := true

# select coloring (uncomment if you don't like coloring of the output)
#COLOR_TTY := false
