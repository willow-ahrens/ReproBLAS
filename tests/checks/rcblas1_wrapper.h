#ifndef RCBLAS1_WRAPPER_H
#define RCBLAS1_WRAPPER_H

#include <reproBLAS.h>
#include <binnedBLAS.h>
#include <binned.h>
#include "../../config.h"

#define wrap_RCSUM  0
#define wrap_RSCASUM 1
#define wrap_RSCNRM2 2
#define wrap_RCDOTU  3
#define wrap_RCDOTC  4
static const int wrap_rcblas1_n_names = 5;
static const char* wrap_rcblas1_names[] = {"rcsum",
                                           "rscasum",
                                           "rscnrm2",
                                           "rcdotu",
                                           "rcdotc"};
static const char* wrap_rcblas1_descs[] = {"rcsum",
                                           "rscasum",
                                           "rscnrm2",
                                           "rcdotu",
                                           "rcdotc"};

typedef float complex (*wrap_rcblas1)(int, float complex*, int, float complex*, int);
typedef void (*wrap_ciblas1)(int, float complex*, int, float complex*, int, float_complex_binned*);


float complex wrap_rcsum(int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  float complex ret;
  reproBLAS_rcsum_sub(N, x, incx, &ret);
  return ret;
}

void wrap_cbcsum(int N, float complex *x, int incx, float complex *y, int incy, float_complex_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_cbcsum(SIDEFAULTFOLD, N, x, incx, z);
}

float complex wrap_rscasum(int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  return (float complex)reproBLAS_rscasum(N, x, incx);
}

void wrap_sbcasum(int N, float complex *x, int incx, float complex *y, int incy, float_complex_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_smcasum(SIDEFAULTFOLD, N, x, incx, z, 2, z + 2 * SIDEFAULTFOLD, 2);
}

float complex wrap_rcdotu(int N, float complex *x, int incx, float complex *y, int incy) {
  float complex ret;
  reproBLAS_rcdotu_sub(N, x, incx, y, incy, &ret);
  return ret;
}

void wrap_cbcdotu(int N, float complex *x, int incx, float complex *y, int incy, float_complex_binned *z) {
  binnedBLAS_cbcdotu(SIDEFAULTFOLD, N, x, incx, y, incy, z);
}

float complex wrap_rcdotc(int N, float complex *x, int incx, float complex *y, int incy) {
  float complex ret;
  reproBLAS_rcdotc_sub(N, x, incx, y, incy, &ret);
  return ret;
}

void wrap_cbcdotc(int N, float complex *x, int incx, float complex *y, int incy, float_complex_binned *z) {
  binnedBLAS_cbcdotc(SIDEFAULTFOLD, N, x, incx, y, incy, z);
}

float complex wrap_rscnrm2(int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  return reproBLAS_rscnrm2(N, x, incx);
}

void wrap_sicnrm(int N, float complex *x, int incx, float complex *y, int incy, float_complex_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_smcssq(SIDEFAULTFOLD, N, x, incx, 0.0, z, 2, z + 2 * SIDEFAULTFOLD, 2);
}

wrap_rcblas1 wrap_rcblas1_func(int func) {
  switch(func){
    case wrap_RCSUM:
      return wrap_rcsum;
    case wrap_RSCASUM:
      return wrap_rscasum;
    case wrap_RSCNRM2:
      return wrap_rscnrm2;
    case wrap_RCDOTU:
      return wrap_rcdotu;
    case wrap_RCDOTC:
      return wrap_rcdotc;
  }
  return NULL;
}

wrap_ciblas1 wrap_ciblas1_func(int func) {
  switch(func){
    case wrap_RCSUM:
      return wrap_cbcsum;
    case wrap_RSCASUM:
      return wrap_sbcasum;
    case wrap_RSCNRM2:
      return wrap_sicnrm;
    case wrap_RCDOTU:
      return wrap_cbcdotu;
    case wrap_RCDOTC:
      return wrap_cbcdotc;
  }
  return NULL;
}

#endif
