#ifndef RDBLAS1_WRAPPER_H
#define RDBLAS1_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>

typedef enum wrap_daugsum_func {
  wrap_daugsum_RDSUM = 0,
  wrap_daugsum_RDASUM,
  wrap_daugsum_RDNRM2,
  wrap_daugsum_RDDOT,
  wrap_daugsum_DIDIADD,
  wrap_daugsum_DIDADD,
  wrap_daugsum_DIDDEPOSIT
} wrap_daugsum_func_t;

typedef double (*wrap_daugsum)(int, int, double*, int, double*, int);
typedef void (*wrap_diaugsum)(int, int, double*, int, double*, int, double_indexed*);
static const int wrap_daugsum_func_n_names = 7;
static const char* wrap_daugsum_func_names[] = {"rdsum",
                                                "rdasum",
                                                "rdnrm2",
                                                "rddot",
                                                "didiadd",
                                                "didadd",
                                                "diddeposit"};
static const char* wrap_daugsum_func_descs[] = {"rdsum",
                                                "rdasum",
                                                "rdnrm2",
                                                "rddot",
                                                "didiadd",
                                                "didadd",
                                                "diddeposit"};

double wrap_rdsum(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rdsum(N, x, incx);
  }else{
    double_indexed *ires = dialloc(fold);
    disetzero(fold, ires);
    didsum(fold, N, x, incx, ires);
    double res = ddiconv(fold, ires);
    free(ires);
    return res;
  }
}

void wrap_didsum(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  didsum(fold, N, x, incx, z);
}

double wrap_rdasum(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rdasum(N, x, incx);
  }else{
    double_indexed *ires = dialloc(fold);
    disetzero(fold, ires);
    didasum(fold, N, x, incx, ires);
    double res = ddiconv(fold, ires);
    free(ires);
    return res;
  }
}

void wrap_didasum(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  didasum(fold, N, x, incx, z);
}

double wrap_rdnrm2(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rdnrm2(N, x, incx);
  }else{
    double_indexed *ires = dialloc(fold);
    disetzero(fold, ires);
    double scale = didnrm(fold, N, x, incx, ires);
    double res = ddiconv(fold, ires);
    free(ires);
    return scale * sqrt(res);
  }
}

void wrap_didnrm(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  didnrm(fold, N, x, incx, z);
}

double wrap_rddot(int fold, int N, double *x, int incx, double *y, int incy) {
  if(fold == DEFAULT_FOLD){
    return rddot(N, x, incx, y, incy);
  }else{
    double_indexed *ires = dialloc(fold);
    disetzero(fold, ires);
    diddot(fold, N, x, incx, y, incy, ires);
    double res = ddiconv(fold, ires);
    free(ires);
    return res;
  }
}

void wrap_diddot(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  diddot(fold, N, x, incx, y, incy, z);
}

double wrap_rdidiadd(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  double_indexed *ires = dialloc(fold);
  double_indexed *itmp = dialloc(fold);
  disetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    didconv(fold, x[i * incx], itmp);
    didiadd(fold, itmp, ires);
  }
  double res = ddiconv(fold, ires);
  free(ires);
  free(itmp);
  return res;
}

void wrap_didiadd(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  double_indexed *itmp = dialloc(fold);
  int i;
  for(i = 0; i < N; i++){
    didconv(fold, x[i * incx], itmp);
    didiadd(fold, itmp, z);
  }
  free(itmp);
}

double wrap_rdidadd(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  double_indexed *ires = dialloc(fold);
  disetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    didadd(fold, x[i * incx], ires);
  }
  double res = ddiconv(fold, ires);
  free(ires);
  return res;
}

void wrap_didadd(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  int i;
  for(i = 0; i < N; i++){
    didadd(fold, x[i * incx], z);
  }
}

double wrap_rdiddeposit(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  double_indexed *ires = dialloc(fold);
  disetzero(fold, ires);
  double amax = damax(N, x, incx);
  didupdate(fold, amax, ires);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= dicapacity()){
      direnorm(fold, ires);
      j = 0;
    }
    diddeposit(fold, x[i * incx], ires);
    j++;
  }
  direnorm(fold, ires);
  double res = ddiconv(fold, ires);
  free(ires);
  return res;
}

void wrap_diddeposit(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  double amax = damax(N, x, incx);
  didupdate(fold, amax, z);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= dicapacity()){
      direnorm(fold, z);
      j = 0;
    }
    diddeposit(fold, x[i * incx], z);
    j++;
  }
  direnorm(fold, z);
}

wrap_daugsum wrap_daugsum_func(wrap_daugsum_func_t func) {
  switch(func){
    case wrap_daugsum_RDSUM:
      return wrap_rdsum;
    case wrap_daugsum_RDASUM:
      return wrap_rdasum;
    case wrap_daugsum_RDNRM2:
      return wrap_rdnrm2;
    case wrap_daugsum_RDDOT:
      return wrap_rddot;
    case wrap_daugsum_DIDIADD:
      return wrap_rdidiadd;
    case wrap_daugsum_DIDADD:
      return wrap_rdidadd;
    case wrap_daugsum_DIDDEPOSIT:
      return wrap_rdiddeposit;
  }
  return NULL;
}

wrap_diaugsum wrap_diaugsum_func(wrap_daugsum_func_t func) {
  switch(func){
    case wrap_daugsum_RDSUM:
      return wrap_didsum;
    case wrap_daugsum_RDASUM:
      return wrap_didasum;
    case wrap_daugsum_RDNRM2:
      return wrap_didnrm;
    case wrap_daugsum_RDDOT:
      return wrap_diddot;
    case wrap_daugsum_DIDIADD:
      return wrap_didiadd;
    case wrap_daugsum_DIDADD:
      return wrap_didadd;
    case wrap_daugsum_DIDDEPOSIT:
      return wrap_diddeposit;
  }
  return NULL;
}

#endif
