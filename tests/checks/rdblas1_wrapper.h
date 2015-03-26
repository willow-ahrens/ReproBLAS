#ifndef _RDBLAS1_WRAPPER_H
#define _RDBLAS1_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>

#define wrap_RDSUM  0
#define wrap_RDASUM 1
#define wrap_RDNRM2 2
#define wrap_RDDOT  3

typedef double (*wrap_rdblas1)(int, double*, int, double*, int);
typedef Idouble (*wrap_Idblas1)(int, double*, int, double*, int);
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
  return rdsum(N, x, incx);
}

Idouble wrap_dsumI(int N, double *x, int incx, double *y, int incy) {
  return dsumI(N, x, incx);
}

double wrap_rdasum(int N, double *x, int incx, double *y, int incy) {
  return rdasum(N, x, incx);
}

Idouble wrap_dasumI(int N, double *x, int incx, double *y, int incy) {
  return dasumI(N, x, incx);
}

double wrap_rddot(int N, double *x, int incx, double *y, int incy) {
  return rddot(N, x, incx, y, incy);
}

Idouble wrap_ddotI(int N, double *x, int incx, double *y, int incy) {
  return ddotI(N, x, incx, y, incy);
}

double wrap_rdnrm2(int N, double *x, int incx, double *y, int incy) {
  return rdnrm2(N, x, incx);
}

Idouble wrap_dnrm2I(int N, double *x, int incx, double *y, int incy) {
  Idouble nrm2;
  dISetZero(nrm2);
  dnrm2I(N, x, incx, &nrm2);
  return nrm2;
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

wrap_Idblas1 wrap_Idblas1_func(int func) {
  switch(func){
    case wrap_RDSUM:
      return wrap_dsumI;
    case wrap_RDASUM:
      return wrap_dasumI;
    case wrap_RDNRM2:
      return wrap_dnrm2I;
    case wrap_RDDOT:
      return wrap_ddotI;
  }
  return NULL;
}

#endif
