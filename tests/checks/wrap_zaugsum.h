#ifndef ZAUGSUM_WRAPPER_H
#define ZAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>

typedef enum wrap_zaugsum_func {
  wrap_zaugsum_RZSUM = 0,
  wrap_zaugsum_RDZASUM,
  wrap_zaugsum_RDZNRM2,
  wrap_zaugsum_RZDOTU,
  wrap_zaugsum_RZDOTC,
  wrap_zaugsum_ZIZIADD,
  wrap_zaugsum_ZIZADD,
  wrap_zaugsum_ZIZDEPOSIT
} wrap_zaugsum_func_t;

typedef double complex (*wrap_zaugsum)(int, int, double complex*, int, double complex*, int);
typedef void (*wrap_ziaugsum)(int, int, double complex*, int, double complex*, int, double_complex_indexed*);
static const int wrap_zaugsum_func_n_names = 8;
static const char* wrap_zaugsum_func_names[] = {"rzsum",
                                                "rdzasum",
                                                "rdznrm2",
                                                "rzdotu",
                                                "rzdotc",
                                                "ziziadd",
                                                "zizadd",
                                                "zizdeposit"};
static const char* wrap_zaugsum_func_descs[] = {"rzsum",
                                                "rdzasum",
                                                "rdznrm2",
                                                "rzdotu",
                                                "rzdotc",
                                                "ziziadd",
                                                "zizadd",
                                                "zizdeposit"};

double complex wrap_rzsum(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    double complex res;
    rzsum_sub(N, x, incx, &res);
    return res;
  }else{
    double_complex_indexed *ires = zialloc(fold);
    zisetzero(fold, ires);
    zizsum(fold, N, x, incx, ires);
    double complex res;
    zziconv_sub(fold, ires, &res);
    free(ires);
    return res;
  }
}

void wrap_zizsum(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  zizsum(fold, N, x, incx, z);
}

double complex wrap_rdzasum(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rdzasum(N, x, incx);
  }else{
    double_indexed *ires = dialloc(fold);
    disetzero(fold, ires);
    dizasum(fold, N, x, incx, ires);
    double complex res = ddiconv(fold, ires);
    free(ires);
    return res;
  }
}

void wrap_dizasum(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  dmzasum(fold, N, x, incx, z, 2, z + 2 * fold, 2);
}

double complex wrap_rdznrm2(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rdznrm2(N, x, incx);
  }else{
    double_indexed *ires = dialloc(fold);
    disetzero(fold, ires);
    double complex scale = diznrm(fold, N, x, incx, ires);
    double complex res = ddiconv(fold, ires);
    free(ires);
    return scale * sqrt(res);
  }
}

void wrap_diznrm(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  dmznrm(fold, N, x, incx, z, 2, z + 2 * fold, 2);
}

double complex wrap_rzdotu(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  if(fold == DEFAULT_FOLD){
    double complex res;
    rzdotu_sub(N, x, incx, y, incy, &res);
    return res;
  }else{
    double_complex_indexed *ires = zialloc(fold);
    zisetzero(fold, ires);
    zizdotu(fold, N, x, incx, y, incy, ires);
    double complex res;
    zziconv_sub(fold, ires, &res);
    free(ires);
    return res;
  }
}

void wrap_zizdotu(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  zizdotu(fold, N, x, incx, y, incy, z);
}

double complex wrap_rzdotc(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  if(fold == DEFAULT_FOLD){
    double complex res;
    rzdotc_sub(N, x, incx, y, incy, &res);
    return res;
  }else{
    double_complex_indexed *ires = zialloc(fold);
    zisetzero(fold, ires);
    zizdotc(fold, N, x, incx, y, incy, ires);
    double complex res;
    zziconv_sub(fold, ires, &res);
    free(ires);
    return res;
  }
}

void wrap_zizdotc(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  zizdotc(fold, N, x, incx, y, incy, z);
}

double complex wrap_rziziadd(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  double_complex_indexed *ires = zialloc(fold);
  double_complex_indexed *itmp = zialloc(fold);
  zisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    zizconv(fold, x + i * incx, itmp);
    ziziadd(fold, itmp, ires);
  }
  double complex res;
  zziconv_sub(fold, ires, &res);
  free(ires);
  free(itmp);
  return res;
}

void wrap_ziziadd(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  double_complex_indexed *itmp = zialloc(fold);
  int i;
  for(i = 0; i < N; i++){
    zizconv(fold, x + i * incx, itmp);
    ziziadd(fold, itmp, z);
  }
  free(itmp);
}

double complex wrap_rzizadd(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  double_complex_indexed *ires = zialloc(fold);
  zisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    zizadd(fold, x + i * incx, ires);
  }
  double complex res;
  zziconv_sub(fold, ires, &res);
  free(ires);
  return res;
}

void wrap_zizadd(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  int i;
  for(i = 0; i < N; i++){
    zizadd(fold, x + i * incx, z);
  }
}

double complex wrap_rzizdeposit(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  double_complex_indexed *ires = zialloc(fold);
  zisetzero(fold, ires);
  double complex amax;
  zamax_sub(N, x, incx, &amax);
  zizupdate(fold, &amax, ires);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= dicapacity()){
      zirenorm(fold, ires);
      j = 0;
    }
    zizdeposit(fold, x + i * incx, ires);
    j++;
  }
  zirenorm(fold, ires);
  double complex res;
  zziconv_sub(fold, ires, &res);
  free(ires);
  return res;
}

void wrap_zizdeposit(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  double complex amax;
  zamax_sub(N, x, incx, &amax);
  zizupdate(fold, &amax, z);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= dicapacity()){
      zirenorm(fold, z);
      j = 0;
    }
    zizdeposit(fold, x + i * incx, z);
    j++;
  }
  zirenorm(fold, z);
}

wrap_zaugsum wrap_zaugsum_func(wrap_zaugsum_func_t func) {
  switch(func){
    case wrap_zaugsum_RZSUM:
      return wrap_rzsum;
    case wrap_zaugsum_RDZASUM:
      return wrap_rdzasum;
    case wrap_zaugsum_RDZNRM2:
      return wrap_rdznrm2;
    case wrap_zaugsum_RZDOTU:
      return wrap_rzdotu;
    case wrap_zaugsum_RZDOTC:
      return wrap_rzdotc;
    case wrap_zaugsum_ZIZIADD:
      return wrap_rziziadd;
    case wrap_zaugsum_ZIZADD:
      return wrap_rzizadd;
    case wrap_zaugsum_ZIZDEPOSIT:
      return wrap_rzizdeposit;
  }
  return NULL;
}

wrap_ziaugsum wrap_ziaugsum_func(wrap_zaugsum_func_t func) {
  switch(func){
    case wrap_zaugsum_RZSUM:
      return wrap_zizsum;
    case wrap_zaugsum_RDZASUM:
      return wrap_dizasum;
    case wrap_zaugsum_RDZNRM2:
      return wrap_diznrm;
    case wrap_zaugsum_RZDOTU:
      return wrap_zizdotu;
    case wrap_zaugsum_RZDOTC:
      return wrap_zizdotc;
    case wrap_zaugsum_ZIZIADD:
      return wrap_ziziadd;
    case wrap_zaugsum_ZIZADD:
      return wrap_zizadd;
    case wrap_zaugsum_ZIZDEPOSIT:
      return wrap_zizdeposit;
  }
  return NULL;
}

#endif
