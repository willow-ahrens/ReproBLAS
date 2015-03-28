#ifndef _RZBLAS1_WRAPPER_H
#define _RZBLAS1_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>

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
typedef I_double_Complex (*wrap_Izblas1)(int, double complex*, int, double complex*, int);

double complex wrap_rzsum(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  return rzsum(N, x, incx);
}

I_double_Complex wrap_zsumI(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  return zsumI(N, x, incx);
}

double complex wrap_rdzasum(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  return (double complex)rdzasum(N, x, incx);
}

I_double_Complex wrap_dzasumI(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  I_double_Complex casum;
  Idouble asum = dzasumI(N, x, incx);
  zdISet(casum, asum);
  return casum;
}

double complex wrap_rzdotu(int N, double complex *x, int incx, double complex *y, int incy) {
  return rzdotu(N, x, incx, y, incy);
}

I_double_Complex wrap_zdotuI(int N, double complex *x, int incx, double complex *y, int incy) {
  return zdotuI(N, x, incx, y, incy);
}

double complex wrap_rzdotc(int N, double complex *x, int incx, double complex *y, int incy) {
  return rzdotc(N, x, incx, y, incy);
}

I_double_Complex wrap_zdotcI(int N, double complex *x, int incx, double complex *y, int incy) {
  return zdotcI(N, x, incx, y, incy);
}

double complex wrap_rdznrm2(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  return rdznrm2(N, x, incx);
}

I_double_Complex wrap_dznrm2I(int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  I_double_Complex cnrm2;
  Idouble nrm2;
  dISetZero(nrm2);
  dznrm2I(N, x, incx, &nrm2);
  zdISet(cnrm2, nrm2);
  return cnrm2;
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

wrap_Izblas1 wrap_Izblas1_func(int func) {
  switch(func){
    case wrap_RZSUM:
      return wrap_zsumI;
    case wrap_RDZASUM:
      return wrap_dzasumI;
    case wrap_RDZNRM2:
      return wrap_dznrm2I;
    case wrap_RZDOTU:
      return wrap_zdotuI;
    case wrap_RZDOTC:
      return wrap_zdotcI;
  }
  return NULL;
}

#endif
