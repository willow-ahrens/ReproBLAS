#ifndef _RDBLAS1_WRAPPER_H
#define _RDBLAS1_WRAPPER_H

#include <rblas.h>
#include <IndexedFP.h>

#define verify_RDSUM  0
#define verify_RDASUM 1
#define verify_RDNRM2 2
#define verify_RDDOT  3

typedef double (*wrap_rdblas1)(int, double*, int, double*, int);
typedef Idouble (*wrap_Idblas1)(int, double*, int, double*, int);

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

const char* wrap_rdblas1_name(int func) {
  switch(func){
    case verify_RDSUM:
      return "rdsum";
    case verify_RDASUM:
      return "rdasum";
    case verify_RDNRM2:
      return "rdnrm2";
    case verify_RDDOT:
      return "rddot";
  }
  return "";
}

wrap_rdblas1 wrap_rdblas1_func(int func) {
  switch(func){
    case verify_RDSUM:
      return wrap_rdsum;
    case verify_RDASUM:
      return wrap_rdasum;
    case verify_RDNRM2:
      return wrap_rdnrm2;
    case verify_RDDOT:
      return wrap_rddot;
  }
  return NULL;
}

wrap_Idblas1 wrap_Idblas1_func(int func) {
  switch(func){
    case verify_RDSUM:
      return wrap_dsumI;
    case verify_RDASUM:
      return wrap_dasumI;
    case verify_RDNRM2:
      return wrap_dnrm2I;
    case verify_RDDOT:
      return wrap_ddotI;
  }
  return NULL;
}

#endif
