TARGETS := libtestxblas.a
SUBDIRS :=

LIBTESTXBLAS := $(OBJPATH)/libtestxblas.a

libtestxblas.a_DEPS = testgen_BLAS_ddot.o testgen_BLAS_zdot.o testgen_BLAS_sdot.o testgen_BLAS_cdot.o
