#ifndef _RSBLAS1_WRAPPER_H
#define _RSBLAS1_WRAPPER_H

#include <rblas.h>
#include <IndexedFP.h>

#define verify_RSSUM  0
#define verify_RSASUM 1
#define verify_RSNRM2 2
#define verify_RSDOT  3

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

const char* wrap_rsblas1_name(int func) {
  switch(func){
    case verify_RSSUM:
      return "rssum";
    case verify_RSASUM:
      return "rsasum";
    case verify_RSNRM2:
      return "rsnrm2";
    case verify_RSDOT:
      return "rsdot";
  }
  return "";
}

wrap_rsblas1 wrap_rsblas1_func(int func) {
  switch(func){
    case verify_RSSUM:
      return wrap_rssum;
    case verify_RSASUM:
      return wrap_rsasum;
    case verify_RSNRM2:
      return wrap_rsnrm2;
    case verify_RSDOT:
      return wrap_rsdot;
  }
  return NULL;
}

wrap_Isblas1 wrap_Isblas1_func(int func) {
  switch(func){
    case verify_RSSUM:
      return wrap_ssumI;
    case verify_RSASUM:
      return wrap_sasumI;
    case verify_RSNRM2:
      return wrap_snrm2I;
    case verify_RSDOT:
      return wrap_sdotI;
  }
  return NULL;
}

#endif
