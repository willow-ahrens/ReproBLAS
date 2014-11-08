###############################################
# C COMPILER                                  #
###############################################
# C_COMPILER = clang
C_COMPILER    = gcc
# C_COMPILER    = icc
# C_COMPILER    = CC

###############################################
# MPI C COMPILER                              #
###############################################
# MPIC_COMPILER = mpicc

###############################################
# COMPILING FLAGS                             #
###############################################

CFLAGS        = -O2 -Werror-implicit-function-declaration -march=core-avx-i
#CFLAGS        = -pg -Werror-implicit-function-declaration -march=core-avx-i
#CFLAGS        = -g -Werror-implicit-function-declaration -march=core-avx-i

# for ICC
# CFLAGS        = -O2 -fp-model strict -Werror-implicit-function-declaration

###############################################
# FLAG FOR OPENMP SUPPORT	              #
###############################################
OMP_FLAG     = -fopenmp

# for ICC
# OMP_FLAG     = -openmp

###############################################
# BLAS LIBRARY LOADER FLAG                    #
###############################################

# REFERENCE LIBRARBY
# LDFLAGS      += -lblas

# ATLAS
# LDFLGAS      += -latlas

# INTEL BLAS
# LDFLAGS      += -mkl_sequential

# MAC OS ACCELERATE FRAMEWORK
LDFLAGS      += -framework Accelerate

LDFLAGS += -lm
