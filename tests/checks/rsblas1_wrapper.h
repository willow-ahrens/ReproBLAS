#ifndef RSBLAS1_WRAPPER_H
#define RSBLAS1_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>
#include "../../config.h"

#define wrap_RSSUM  0
#define wrap_RSASUM 1
#define wrap_RSNRM2 2
#define wrap_RSDOT  3
static const int wrap_rsblas1_n_names = 4;
static const char* wrap_rsblas1_names[] = {"rssum",
                                           "rsasum",
                                           "rsnrm2",
                                           "rsdot"};
static const char* wrap_rsblas1_descs[] = {"rssum",
                                           "rdasum",
                                           "rdnrm2",
                                           "rddot"};

typedef float (*wrap_rsblas1)(int, float*, int, float*, int);
typedef void (*wrap_siblas1)(int, float*, int, float*, int, float_indexed*);

float wrap_rssum(int N, float *x, int incx, float *y, int incy) {
  (void)y;
  (void)incy;
  return rssum(N, x, incx);
}

void wrap_sissum(int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  (void)y;
  (void)incy;
  sissum(DEFAULT_FOLD, N, x, incx, z);
}

float wrap_rsasum(int N, float *x, int incx, float *y, int incy) {
  (void)y;
  (void)incy;
  return rsasum(N, x, incx);
}

void wrap_sisasum(int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  (void)y;
  (void)incy;
  sisasum(DEFAULT_FOLD, N, x, incx, z);
}

float wrap_rsdot(int N, float *x, int incx, float *y, int incy) {
  return rsdot(N, x, incx, y, incy);
}

void wrap_sisdot(int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  sisdot(DEFAULT_FOLD, N, x, incx, y, incy, z);
}

float wrap_rsnrm2(int N, float *x, int incx, float *y, int incy) {
  (void)y;
  (void)incy;
  return rsnrm2(N, x, incx);
}

void wrap_sisnrm(int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  (void)y;
  (void)incy;
  sisssq(DEFAULT_FOLD, N, x, incx, 0.0, z);
}

wrap_rsblas1 wrap_rsblas1_func(int func) {
  switch(func){
    case wrap_RSSUM:
      return wrap_rssum;
    case wrap_RSASUM:
      return wrap_rsasum;
    case wrap_RSNRM2:
      return wrap_rsnrm2;
    case wrap_RSDOT:
      return wrap_rsdot;
  }
  return NULL;
}

wrap_siblas1 wrap_siblas1_func(int func) {
  switch(func){
    case wrap_RSSUM:
      return wrap_sissum;
    case wrap_RSASUM:
      return wrap_sisasum;
    case wrap_RSNRM2:
      return wrap_sisnrm;
    case wrap_RSDOT:
      return wrap_sisdot;
  }
  return NULL;
}

#endif
