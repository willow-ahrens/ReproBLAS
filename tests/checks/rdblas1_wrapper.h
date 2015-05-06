#ifndef RDBLAS1_WRAPPER_H
#define RDBLAS1_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>

#define wrap_RDSUM  0
#define wrap_RDASUM 1
#define wrap_RDNRM2 2
#define wrap_RDDOT  3

typedef double (*wrap_rdblas1)(int, double*, int, double*, int);
typedef void (*wrap_diblas1)(int, double*, int, double*, int, double_indexed*);
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
  return rdsum(N, x, incx);
}

void wrap_didsum(int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  didsum(DEFAULT_FOLD, N, x, incx, z);
}

double wrap_rdasum(int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  return rdasum(N, x, incx);
}

void wrap_didasum(int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  didasum(DEFAULT_FOLD, N, x, incx, z);
}

double wrap_rddot(int N, double *x, int incx, double *y, int incy) {
  return rddot(N, x, incx, y, incy);
}

void wrap_diddot(int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  diddot(DEFAULT_FOLD, N, x, incx, y, incy, z);
}

double wrap_rdnrm2(int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  return rdnrm2(N, x, incx);
}

void wrap_didnrm(int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  didnrm(DEFAULT_FOLD, N, x, incx, z);
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
      return wrap_didsum;
    case wrap_RDASUM:
      return wrap_didasum;
    case wrap_RDNRM2:
      return wrap_didnrm;
    case wrap_RDDOT:
      return wrap_diddot;
  }
  return NULL;
}

#endif
