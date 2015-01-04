#ifndef _RCBLAS1_WRAPPER_H
#define _RCBLAS1_WRAPPER_H

#include <rblas.h>
#include <IndexedFP.h>

#define verify_RCSUM  0
#define verify_RSCASUM 1
#define verify_RSCNRM2 2
#define verify_RCDOTU  3
#define verify_RCDOTC  4
static int verify_rcblas1_n_names = 5;
static const char* verify_rcblas1_n_names[] = {"rcsum",
                                               "rscasum",
                                               "rscnrm2",
                                               "rcdotu",
                                               "rcdotc"}
static const char* verify_rcblas1_n_descs[] = {"rcsum",
                                               "rscasum",
                                               "rscnrm2",
                                               "rcdotu",
                                               "rcdotc"}

typedef float complex (*wrap_rcblas1)(int, float complex*, int, float complex*, int);
typedef I_float_Complex (*wrap_Icblas1)(int, float complex*, int, float complex*, int);


float complex wrap_rcsum(int N, float complex *x, int incx, float complex *y, int incy) {
  return rcsum(N, x, incx);
}

I_float_Complex wrap_csumI(int N, float complex *x, int incx, float complex *y, int incy) {
  return csumI(N, x, incx);
}

float complex wrap_rscasum(int N, float complex *x, int incx, float complex *y, int incy) {
  return (float complex)rscasum(N, x, incx);
}

I_float_Complex wrap_scasumI(int N, float complex *x, int incx, float complex *y, int incy) {
  I_float_Complex casum;
  Ifloat asum = scasumI(N, x, incx);
  csISet(casum, asum);
  return casum;
}

float complex wrap_rcdotu(int N, float complex *x, int incx, float complex *y, int incy) {
  return rcdotu(N, x, incx, y, incy);
}

I_float_Complex wrap_cdotuI(int N, float complex *x, int incx, float complex *y, int incy) {
  return cdotuI(N, x, incx, y, incy);
}

float complex wrap_rcdotc(int N, float complex *x, int incx, float complex *y, int incy) {
  return rcdotc(N, x, incx, y, incy);
}

I_float_Complex wrap_cdotcI(int N, float complex *x, int incx, float complex *y, int incy) {
  return cdotcI(N, x, incx, y, incy);
}

float complex wrap_rscnrm2(int N, float complex *x, int incx, float complex *y, int incy) {
  return rscnrm2(N, x, incx);
}

I_float_Complex wrap_scnrm2I(int N, float complex *x, int incx, float complex *y, int incy) {
  I_float_Complex cnrm2;
  Ifloat nrm2;
  sISetZero(nrm2);
  scnrm2I(N, x, incx, &nrm2);
  csISet(cnrm2, nrm2);
  return cnrm2;
}

wrap_rcblas1 wrap_rcblas1_func(int func) {
  switch(func){
    case verify_RCSUM:
      return wrap_rcsum;
    case verify_RSCASUM:
      return wrap_rscasum;
    case verify_RSCNRM2:
      return wrap_rscnrm2;
    case verify_RCDOTU:
      return wrap_rcdotu;
    case verify_RCDOTC:
      return wrap_rcdotc;
  }
  return NULL;
}

wrap_Icblas1 wrap_Icblas1_func(int func) {
  switch(func){
    case verify_RCSUM:
      return wrap_csumI;
    case verify_RSCASUM:
      return wrap_scasumI;
    case verify_RSCNRM2:
      return wrap_scnrm2I;
    case verify_RCDOTU:
      return wrap_cdotuI;
    case verify_RCDOTC:
      return wrap_cdotcI;
  }
  return NULL;
}

#endif
