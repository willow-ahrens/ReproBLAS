#ifndef RZBLAS1_WRAPPER_H
#define RZBLAS1_WRAPPER_H

#include <reproBLAS.h>
#include <binnedBLAS.h>
#include <binned.h>
#include "../../config.h"

#define wrap_RZSUM  0
#define wrap_RDZASUM 1
#define wrap_RDZNRM2 2
#define wrap_RZDOTU  3
#define wrap_RZDOTC  4
static const int wrap_rzblas1_n_names = 5;
static const char* wrap_rzblas1_names[] = {"rzsum",
                                           "rdzasum",
                                           "rdznrm2",
                                           "rzdotu",
                                           "rzdotc"};
static const char* wrap_rzblas1_descs[] = {"rzsum",
                                           "rdzasum",
                                           "rdznrm2",
                                           "rzdotu",
                                           "rzdotc"};

typedef double complex (*wrap_rzblas1)(int, double complex*, int, double complex*, int);
typedef void (*wrap_ziblas1)(int, double complex*, int, double complex*, int, double_complex_binned*);

double complex wrap_rzsum(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  double complex ret;
  reproBLAS_rzsum_sub(N, x, incx, &ret);
  return ret;
}

void wrap_zbzsum(int N, double complex *x, int incx, double complex *y, int incy, double_complex_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_zbzsum(DIDEFAULTFOLD, N, x, incx, z);
}

double complex wrap_rdzasum(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  return (double complex)reproBLAS_rdzasum(N, x, incx);
}

void wrap_dbzasum(int N, double complex *x, int incx, double complex *y, int incy, double_complex_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_dmzasum(DIDEFAULTFOLD, N, x, incx, z, 2, z + 2 * DIDEFAULTFOLD, 2);
}

double complex wrap_rzdotu(int N, double complex *x, int incx, double complex *y, int incy) {
  double complex ret;
  reproBLAS_rzdotu_sub(N, x, incx, y, incy, &ret);
  return ret;
}

void wrap_zbzdotu(int N, double complex *x, int incx, double complex *y, int incy, double_complex_binned *z) {
  binnedBLAS_zbzdotu(DIDEFAULTFOLD, N, x, incx, y, incy, z);
}

double complex wrap_rzdotc(int N, double complex *x, int incx, double complex *y, int incy) {
  double complex ret;
  reproBLAS_rzdotc_sub(N, x, incx, y, incy, &ret);
  return ret;
}

void wrap_zbzdotc(int N, double complex *x, int incx, double complex *y, int incy, double_complex_binned *z) {
  binnedBLAS_zbzdotc(DIDEFAULTFOLD, N, x, incx, y, incy, z);
}

double complex wrap_rdznrm2(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  return reproBLAS_rdznrm2(N, x, incx);
}

void wrap_diznrm(int N, double complex *x, int incx, double complex *y, int incy, double_complex_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_dmzssq(DIDEFAULTFOLD, N, x, incx, 0.0, z, 2, z + 2 * DIDEFAULTFOLD, 2);
}

wrap_rzblas1 wrap_rzblas1_func(int func) {
  switch(func){
    case wrap_RZSUM:
      return wrap_rzsum;
    case wrap_RDZASUM:
      return wrap_rdzasum;
    case wrap_RDZNRM2:
      return wrap_rdznrm2;
    case wrap_RZDOTU:
      return wrap_rzdotu;
    case wrap_RZDOTC:
      return wrap_rzdotc;
  }
  return NULL;
}

wrap_ziblas1 wrap_ziblas1_func(int func) {
  switch(func){
    case wrap_RZSUM:
      return wrap_zbzsum;
    case wrap_RDZASUM:
      return wrap_dbzasum;
    case wrap_RDZNRM2:
      return wrap_diznrm;
    case wrap_RZDOTU:
      return wrap_zbzdotu;
    case wrap_RZDOTC:
      return wrap_zbzdotc;
  }
  return NULL;
}

#endif
