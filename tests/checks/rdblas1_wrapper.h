#ifndef RDBLAS1_WRAPPER_H
#define RDBLAS1_WRAPPER_H

#include <reproBLAS.h>
#include <binnedBLAS.h>
#include <binned.h>
#include "../../config.h"

#define wrap_RDSUM  0
#define wrap_RDASUM 1
#define wrap_RDNRM2 2
#define wrap_RDDOT  3

typedef double (*wrap_rdblas1)(int, double*, int, double*, int);
typedef void (*wrap_diblas1)(int, double*, int, double*, int, double_binned*);
static const int wrap_rdblas1_n_names = 4;
static const char* wrap_rdblas1_names[] = {"rdsum",
                                           "rdasum",
                                           "rdnrm2",
                                           "rddot"};
static const char* wrap_rdblas1_descs[] = {"rdsum",
                                           "rdasum",
                                           "rdnrm2",
                                           "rddot"};

double wrap_rdsum(int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  return reproBLAS_rdsum(N, x, incx);
}

void wrap_dbdsum(int N, double *x, int incx, double *y, int incy, double_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_dbdsum(DIDEFAULTFOLD, N, x, incx, z);
}

double wrap_rdasum(int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  return reproBLAS_rdasum(N, x, incx);
}

void wrap_dbdasum(int N, double *x, int incx, double *y, int incy, double_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_dbdasum(DIDEFAULTFOLD, N, x, incx, z);
}

double wrap_rddot(int N, double *x, int incx, double *y, int incy) {
  return reproBLAS_rddot(N, x, incx, y, incy);
}

void wrap_dbddot(int N, double *x, int incx, double *y, int incy, double_binned *z) {
  binnedBLAS_dbddot(DIDEFAULTFOLD, N, x, incx, y, incy, z);
}

double wrap_rdnrm2(int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  return reproBLAS_rdnrm2(N, x, incx);
}

void wrap_didnrm(int N, double *x, int incx, double *y, int incy, double_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_dbdssq(DIDEFAULTFOLD, N, x, incx, 0.0, z);
}

wrap_rdblas1 wrap_rdblas1_func(int func) {
  switch(func){
    case wrap_RDSUM:
      return wrap_rdsum;
    case wrap_RDASUM:
      return wrap_rdasum;
    case wrap_RDNRM2:
      return wrap_rdnrm2;
    case wrap_RDDOT:
      return wrap_rddot;
  }
  return NULL;
}

wrap_diblas1 wrap_diblas1_func(int func) {
  switch(func){
    case wrap_RDSUM:
      return wrap_dbdsum;
    case wrap_RDASUM:
      return wrap_dbdasum;
    case wrap_RDNRM2:
      return wrap_didnrm;
    case wrap_RDDOT:
      return wrap_dbddot;
  }
  return NULL;
}

#endif
