#ifndef _RSBLAS1_WRAPPER_H
#define _RSBLAS1_WRAPPER_H

#include <rblas.h>
#include <IndexedFP.h>

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
typedef Ifloat (*wrap_Isblas1)(int, float*, int, float*, int);

float wrap_rssum(int N, float *x, int incx, float *y, int incy) {
  return rssum(N, x, incx);
}

Ifloat wrap_ssumI(int N, float *x, int incx, float *y, int incy) {
  return ssumI(N, x, incx);
}

float wrap_rsasum(int N, float *x, int incx, float *y, int incy) {
  return rsasum(N, x, incx);
}

Ifloat wrap_sasumI(int N, float *x, int incx, float *y, int incy) {
  return sasumI(N, x, incx);
}

float wrap_rsdot(int N, float *x, int incx, float *y, int incy) {
  return rsdot(N, x, incx, y, incy);
}

Ifloat wrap_sdotI(int N, float *x, int incx, float *y, int incy) {
  return sdotI(N, x, incx, y, incy);
}

float wrap_rsnrm2(int N, float *x, int incx, float *y, int incy) {
  return rsnrm2(N, x, incx);
}

Ifloat wrap_snrm2I(int N, float *x, int incx, float *y, int incy) {
  Ifloat nrm2;
  sISetZero(nrm2);
  snrm2I(N, x, incx, &nrm2);
  return nrm2;
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

wrap_Isblas1 wrap_Isblas1_func(int func) {
  switch(func){
    case wrap_RSSUM:
      return wrap_ssumI;
    case wrap_RSASUM:
      return wrap_sasumI;
    case wrap_RSNRM2:
      return wrap_snrm2I;
    case wrap_RSDOT:
      return wrap_sdotI;
  }
  return NULL;
}

#endif
