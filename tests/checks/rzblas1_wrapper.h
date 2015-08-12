#ifndef RZBLAS1_WRAPPER_H
#define RZBLAS1_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <idxd.h>
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
typedef void (*wrap_ziblas1)(int, double complex*, int, double complex*, int, double_complex_indexed*);

double complex wrap_rzsum(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  double complex ret;
  rzsum_sub(N, x, incx, &ret);
  return ret;
}

void wrap_zizsum(int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  idxdBLAS_zizsum(DIDEFAULTFOLD, N, x, incx, z);
}

double complex wrap_rdzasum(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  return (double complex)rdzasum(N, x, incx);
}

void wrap_dizasum(int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  idxdBLAS_dmzasum(DIDEFAULTFOLD, N, x, incx, z, 2, z + 2 * DIDEFAULTFOLD, 2);
}

double complex wrap_rzdotu(int N, double complex *x, int incx, double complex *y, int incy) {
  double complex ret;
  rzdotu_sub(N, x, incx, y, incy, &ret);
  return ret;
}

void wrap_zizdotu(int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  idxdBLAS_zizdotu(DIDEFAULTFOLD, N, x, incx, y, incy, z);
}

double complex wrap_rzdotc(int N, double complex *x, int incx, double complex *y, int incy) {
  double complex ret;
  rzdotc_sub(N, x, incx, y, incy, &ret);
  return ret;
}

void wrap_zizdotc(int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  idxdBLAS_zizdotc(DIDEFAULTFOLD, N, x, incx, y, incy, z);
}

double complex wrap_rdznrm2(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  return rdznrm2(N, x, incx);
}

void wrap_diznrm(int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  idxdBLAS_dmzssq(DIDEFAULTFOLD, N, x, incx, 0.0, z, 2, z + 2 * DIDEFAULTFOLD, 2);
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
      return wrap_zizsum;
    case wrap_RDZASUM:
      return wrap_dizasum;
    case wrap_RDZNRM2:
      return wrap_diznrm;
    case wrap_RZDOTU:
      return wrap_zizdotu;
    case wrap_RZDOTC:
      return wrap_zizdotc;
  }
  return NULL;
}

#endif
