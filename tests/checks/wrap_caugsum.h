#ifndef CAUGSUM_WRAPPER_H
#define CAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>

#include "../../src/Common/Common.h"

typedef enum wrap_caugsum_func {
  wrap_caugsum_RCSUM = 0,
  wrap_caugsum_RSCASUM,
  wrap_caugsum_RSCNRM2,
  wrap_caugsum_RCDOTU,
  wrap_caugsum_RCDOTC,
  wrap_caugsum_CICIADD,
  wrap_caugsum_CICADD,
  wrap_caugsum_CICDEPOSIT
} wrap_caugsum_func_t;

typedef float complex (*wrap_caugsum)(int, int, float complex*, int, float complex*, int);
typedef void (*wrap_ciaugsum)(int, int, float complex*, int, float complex*, int, float_complex_indexed*);
static const int wrap_caugsum_func_n_names = 8;
static const char* wrap_caugsum_func_names[] = {"rcsum",
                                                "rscasum",
                                                "rscnrm2",
                                                "rcdotu",
                                                "rcdotc",
                                                "ciciadd",
                                                "cicadd",
                                                "cicdeposit"};
static const char* wrap_caugsum_func_descs[] = {"rcsum",
                                                "rscasum",
                                                "rscnrm2",
                                                "rcdotu",
                                                "rcdotc",
                                                "ciciadd",
                                                "cicadd",
                                                "cicdeposit"};

float complex wrap_rcsum(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    float complex res;
    rcsum_sub(N, x, incx, &res);
    return res;
  }else{
    float_complex_indexed *ires = cialloc(fold);
    cisetzero(fold, ires);
    cicsum(fold, N, x, incx, ires);
    float complex res;
    cciconv_sub(fold, ires, &res);
    free(ires);
    return res;
  }
}

void wrap_cicsum(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  cicsum(fold, N, x, incx, c);
}

float complex wrap_rdcasum(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rscasum(N, x, incx);
  }else{
    float_indexed *ires = sialloc(fold);
    sisetzero(fold, ires);
    sicasum(fold, N, x, incx, ires);
    float complex res = ssiconv(fold, ires);
    free(ires);
    return res;
  }
}

void wrap_sicasum(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  smcasum(fold, N, x, incx, c, 2, c + 2 * fold, 2);
}

float complex wrap_rdcnrm2(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rscnrm2(N, x, incx);
  }else{
    float_indexed *ires = sialloc(fold);
    sisetzero(fold, ires);
    float complex scale = sicnrm(fold, N, x, incx, ires);
    float complex res = ssiconv(fold, ires);
    free(ires);
    return scale * sqrt(res);
  }
}

void wrap_sicnrm(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  smcnrm(fold, N, x, incx, c, 2, c + 2 * fold, 2);
}

float complex wrap_rcdotu(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  if(fold == DEFAULT_FOLD){
    float complex res;
    rcdotu_sub(N, x, incx, y, incy, &res);
    return res;
  }else{
    float_complex_indexed *ires = cialloc(fold);
    cisetzero(fold, ires);
    cicdotu(fold, N, x, incx, y, incy, ires);
    float complex res;
    cciconv_sub(fold, ires, &res);
    free(ires);
    return res;
  }
}

void wrap_cicdotu(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  cicdotu(fold, N, x, incx, y, incy, c);
}

float complex wrap_rcdotc(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  if(fold == DEFAULT_FOLD){
    float complex res;
    rcdotc_sub(N, x, incx, y, incy, &res);
    return res;
  }else{
    float_complex_indexed *ires = cialloc(fold);
    cisetzero(fold, ires);
    cicdotc(fold, N, x, incx, y, incy, ires);
    float complex res;
    cciconv_sub(fold, ires, &res);
    free(ires);
    return res;
  }
}

void wrap_cicdotc(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  cicdotc(fold, N, x, incx, y, incy, c);
}

float complex wrap_rciciadd(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  float_complex_indexed *ires = cialloc(fold);
  float_complex_indexed *itmp = cialloc(fold);
  cisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    cicconv(fold, x + i * incx, itmp);
    ciciadd(fold, itmp, ires);
  }
  float complex res;
  cciconv_sub(fold, ires, &res);
  free(ires);
  free(itmp);
  return res;
}

void wrap_ciciadd(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  float_complex_indexed *itmp = cialloc(fold);
  int i;
  for(i = 0; i < N; i++){
    cicconv(fold, x + i * incx, itmp);
    ciciadd(fold, itmp, c);
  }
  free(itmp);
}

float complex wrap_rcicadd(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  float_complex_indexed *ires = cialloc(fold);
  cisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    cicadd(fold, x + i * incx, ires);
  }
  float complex res;
  cciconv_sub(fold, ires, &res);
  free(ires);
  return res;
}

void wrap_cicadd(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  int i;
  for(i = 0; i < N; i++){
    cicadd(fold, x + i * incx, c);
  }
}

float complex wrap_rcicdeposit(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  float_complex_indexed *ires = cialloc(fold);
  cisetzero(fold, ires);
  float complex amax;
  camax_sub(N, x, incx, &amax);
  cicupdate(fold, &amax, ires);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= sicapacity()){
      cirenorm(fold, ires);
      j = 0;
    }
    cicdeposit(fold, x + i * incx, ires);
    j++;
  }
  cirenorm(fold, ires);
  float complex res;
  cciconv_sub(fold, ires, &res);
  free(ires);
  return res;
}

void wrap_cicdeposit(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  float complex amax;
  camax_sub(N, x, incx, &amax);
  cicupdate(fold, &amax, c);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= sicapacity()){
      cirenorm(fold, c);
      j = 0;
    }
    cicdeposit(fold, x + i * incx, c);
    j++;
  }
  cirenorm(fold, c);
}

wrap_caugsum wrap_caugsum_func(wrap_caugsum_func_t func) {
  switch(func){
    case wrap_caugsum_RCSUM:
      return wrap_rcsum;
    case wrap_caugsum_RSCASUM:
      return wrap_rdcasum;
    case wrap_caugsum_RSCNRM2:
      return wrap_rdcnrm2;
    case wrap_caugsum_RCDOTU:
      return wrap_rcdotu;
    case wrap_caugsum_RCDOTC:
      return wrap_rcdotc;
    case wrap_caugsum_CICIADD:
      return wrap_rciciadd;
    case wrap_caugsum_CICADD:
      return wrap_rcicadd;
    case wrap_caugsum_CICDEPOSIT:
      return wrap_rcicdeposit;
  }
  return NULL;
}

wrap_ciaugsum wrap_ciaugsum_func(wrap_caugsum_func_t func) {
  switch(func){
    case wrap_caugsum_RCSUM:
      return wrap_cicsum;
    case wrap_caugsum_RSCASUM:
      return wrap_sicasum;
    case wrap_caugsum_RSCNRM2:
      return wrap_sicnrm;
    case wrap_caugsum_RCDOTU:
      return wrap_cicdotu;
    case wrap_caugsum_RCDOTC:
      return wrap_cicdotc;
    case wrap_caugsum_CICIADD:
      return wrap_ciciadd;
    case wrap_caugsum_CICADD:
      return wrap_cicadd;
    case wrap_caugsum_CICDEPOSIT:
      return wrap_cicdeposit;
  }
  return NULL;
}

float wrap_caugsum_bound(int fold, int N, wrap_caugsum_func_t func, float *X, int incX, float *Y, int incY, float res, float ref){
  switch(func){
    case wrap_caugsum_RCSUM:
    case wrap_caugsum_CICIADD:
    case wrap_caugsum_CICADD:
    case wrap_caugsum_CICDEPOSIT:
      {
        float complex amax;
        camax_sub(N, X, incX, &amax);
        return sibound(fold, N, crealf(amax)) + sibound(fold, N, cimagf(amax)) * I;
      }
    case wrap_caugsum_RSCASUM:
      {
        float complex amax;
        camax_sub(N, X, incX, &amax);
        return sibound(fold, N, MAX(crealf(amax), cimagf(amax)));
      }
    case wrap_caugsum_RSCNRM2:
      {
        float amax;
        camaxm_sub(N, X, incX, X, incX, &amax);
        return sibound(fold, N, crealf(amax)) * (crealf(amax) / (res + ref));
      }
    case wrap_caugsum_RCDOTU:
    case wrap_caugsum_RCDOTC:
      {
        float complex amaxm;
        camaxm_sub(N, X, incX, Y, incY, &amaxm);
        return sibound(fold, N, crealf(amaxm)) + sibound(fold, N, cimagf(amaxm)) * I;
      }
  }
}

#endif
