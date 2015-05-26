#ifndef SAUGSUM_WRAPPER_H
#define SAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>

#include "../common/test_util.h"

typedef enum wrap_saugsum_func {
  wrap_saugsum_RSSUM = 0,
  wrap_saugsum_RSASUM,
  wrap_saugsum_RSNRM2,
  wrap_saugsum_RSDOT,
  wrap_saugsum_SISIADD,
  wrap_saugsum_SISADD,
  wrap_saugsum_SISDEPOSIT
} wrap_saugsum_func_t;

typedef float (*wrap_saugsum)(int, int, float*, int, float*, int);
typedef void (*wrap_siaugsum)(int, int, float*, int, float*, int, float_indexed*);
static const int wrap_saugsum_func_n_names = 7;
static const char* wrap_saugsum_func_names[] = {"rssum",
                                                "rsasum",
                                                "rsnrm2",
                                                "rsdot",
                                                "sisiadd",
                                                "sisadd",
                                                "sisdeposit"};
static const char* wrap_saugsum_func_descs[] = {"rssum",
                                                "rsasum",
                                                "rsnrm2",
                                                "rsdot",
                                                "sisiadd",
                                                "sisadd",
                                                "sisdeposit"};

float wrap_rssum(int fold, int N, float *x, int incx, float *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rssum(N, x, incx);
  }else{
    float_indexed *ires = sialloc(fold);
    sisetzero(fold, ires);
    sissum(fold, N, x, incx, ires);
    float res = ssiconv(fold, ires);
    free(ires);
    return res;
  }
}

void wrap_sissum(int fold, int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  (void)y;
  (void)incy;
  sissum(fold, N, x, incx, z);
}

float wrap_rsasum(int fold, int N, float *x, int incx, float *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rsasum(N, x, incx);
  }else{
    float_indexed *ires = sialloc(fold);
    sisetzero(fold, ires);
    sisasum(fold, N, x, incx, ires);
    float res = ssiconv(fold, ires);
    free(ires);
    return res;
  }
}

void wrap_sisasum(int fold, int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  (void)y;
  (void)incy;
  sisasum(fold, N, x, incx, z);
}

float wrap_rsnrm2(int fold, int N, float *x, int incx, float *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DEFAULT_FOLD){
    return rsnrm2(N, x, incx);
  }else{
    float_indexed *ires = sialloc(fold);
    sisetzero(fold, ires);
    float scale = sisnrm(fold, N, x, incx, ires);
    float res = ssiconv(fold, ires);
    free(ires);
    return scale * sqrt(res);
  }
}

void wrap_sisnrm(int fold, int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  (void)y;
  (void)incy;
  sisnrm(fold, N, x, incx, z);
}

float wrap_rsdot(int fold, int N, float *x, int incx, float *y, int incy) {
  if(fold == DEFAULT_FOLD){
    return rsdot(N, x, incx, y, incy);
  }else{
    float_indexed *ires = sialloc(fold);
    sisetzero(fold, ires);
    sisdot(fold, N, x, incx, y, incy, ires);
    float res = ssiconv(fold, ires);
    free(ires);
    return res;
  }
}

void wrap_sisdot(int fold, int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  sisdot(fold, N, x, incx, y, incy, z);
}

float wrap_rsisiadd(int fold, int N, float *x, int incx, float *y, int incy) {
  (void)y;
  (void)incy;
  float_indexed *ires = sialloc(fold);
  float_indexed *itmp = sialloc(fold);
  sisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    sisconv(fold, x[i * incx], itmp);
    sisiadd(fold, itmp, ires);
  }
  float res = ssiconv(fold, ires);
  free(ires);
  free(itmp);
  return res;
}

void wrap_sisiadd(int fold, int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  (void)y;
  (void)incy;
  float_indexed *itmp = sialloc(fold);
  int i;
  for(i = 0; i < N; i++){
    sisconv(fold, x[i * incx], itmp);
    sisiadd(fold, itmp, z);
  }
  free(itmp);
}

float wrap_rsisadd(int fold, int N, float *x, int incx, float *y, int incy) {
  (void)y;
  (void)incy;
  float_indexed *ires = sialloc(fold);
  sisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    sisadd(fold, x[i * incx], ires);
  }
  float res = ssiconv(fold, ires);
  free(ires);
  return res;
}

void wrap_sisadd(int fold, int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  (void)y;
  (void)incy;
  int i;
  for(i = 0; i < N; i++){
    sisadd(fold, x[i * incx], z);
  }
}

float wrap_rsisdeposit(int fold, int N, float *x, int incx, float *y, int incy) {
  (void)y;
  (void)incy;
  float_indexed *ires = sialloc(fold);
  sisetzero(fold, ires);
  float amax = samax(N, x, incx);
  sisupdate(fold, amax, ires);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= sicapacity()){
      sirenorm(fold, ires);
      j = 0;
    }
    sisdeposit(fold, x[i * incx], ires);
    j++;
  }
  sirenorm(fold, ires);
  float res = ssiconv(fold, ires);
  free(ires);
  return res;
}

void wrap_sisdeposit(int fold, int N, float *x, int incx, float *y, int incy, float_indexed *z) {
  (void)y;
  (void)incy;
  float amax = samax(N, x, incx);
  sisupdate(fold, amax, z);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= sicapacity()){
      sirenorm(fold, z);
      j = 0;
    }
    sisdeposit(fold, x[i * incx], z);
    j++;
  }
  sirenorm(fold, z);
}

wrap_saugsum wrap_saugsum_func(wrap_saugsum_func_t func) {
  switch(func){
    case wrap_saugsum_RSSUM:
      return wrap_rssum;
    case wrap_saugsum_RSASUM:
      return wrap_rsasum;
    case wrap_saugsum_RSNRM2:
      return wrap_rsnrm2;
    case wrap_saugsum_RSDOT:
      return wrap_rsdot;
    case wrap_saugsum_SISIADD:
      return wrap_rsisiadd;
    case wrap_saugsum_SISADD:
      return wrap_rsisadd;
    case wrap_saugsum_SISDEPOSIT:
      return wrap_rsisdeposit;
  }
  return NULL;
}

wrap_siaugsum wrap_siaugsum_func(wrap_saugsum_func_t func) {
  switch(func){
    case wrap_saugsum_RSSUM:
      return wrap_sissum;
    case wrap_saugsum_RSASUM:
      return wrap_sisasum;
    case wrap_saugsum_RSNRM2:
      return wrap_sisnrm;
    case wrap_saugsum_RSDOT:
      return wrap_sisdot;
    case wrap_saugsum_SISIADD:
      return wrap_sisiadd;
    case wrap_saugsum_SISADD:
      return wrap_sisadd;
    case wrap_saugsum_SISDEPOSIT:
      return wrap_sisdeposit;
  }
  return NULL;
}

float wrap_saugsum_result(int N, wrap_saugsum_func_t func, util_vec_fill_t FillX, double ScaleX, double CondX, util_vec_fill_t FillY, double ScaleY, double CondY){
  float small = 1.0 / (1024.0 * 4.0); // 2^-12
  float big   = 1024.0 * 8.0;  // 2^13
  switch(func){
    case wrap_saugsum_RSSUM:
    case wrap_saugsum_SISIADD:
    case wrap_saugsum_SISADD:
    case wrap_saugsum_SISDEPOSIT:
      switch(FillX){
        case util_Vec_Constant:
          return N * ScaleX;
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
          return ScaleX/0.0;
        case util_Vec_Pos_Neg_Inf:
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          return 0.0/0.0;
        case util_Vec_Pos_Big:
          return (N - 1) * ScaleX * small + ScaleX * big;
        case util_Vec_Pos_Pos_Big:
          return ((N - 2) * ScaleX * small + ScaleX * big) + ScaleX * big;
        case util_Vec_Pos_Neg_Big:
          return (N - 2) * ScaleX * small;
        case util_Vec_Sine:
          return ScaleX - ScaleX;
        default:
          printf("foobar %d %d\n", FillX, util_Vec_Sine);
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX);
          exit(125);
      }
    case wrap_saugsum_RSASUM:
      switch(FillX){
        case util_Vec_Constant:
          return N * fabs(ScaleX);
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
        case util_Vec_Pos_Neg_Inf:
          return fabs(ScaleX)/0.0;
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          return 0.0/0.0;
        case util_Vec_Pos_Big:
          return (N - 1) * fabs(ScaleX) * small + fabs(ScaleX) * big;
        case util_Vec_Pos_Pos_Big:
        case util_Vec_Pos_Neg_Big:
          return ((N - 2) * fabs(ScaleX) * small + fabs(ScaleX) * big) + fabs(ScaleX) * big;
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX);
          exit(125);
      }
    case wrap_saugsum_RSNRM2:
      switch(FillX){
        case util_Vec_Constant:
          return sqrt(N) * ScaleX;
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
        case util_Vec_Pos_Neg_Inf:
          return fabs(ScaleX)/0.0;
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          return 0.0/0.0;
        case util_Vec_Pos_Big:
          return sqrt((N - 1) * small * small + big * big) * fabs(ScaleX);
        case util_Vec_Pos_Pos_Big:
        case util_Vec_Pos_Neg_Big:
          return sqrt(((N - 2) * small * small + big * big) + big * big) * fabs(ScaleX);
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX);
          exit(125);
      }
    case wrap_saugsum_RSDOT:
      switch(FillX){
        case util_Vec_Constant:
          switch(FillY){
            case util_Vec_Constant:
              return N * ScaleX * ScaleY;
            case util_Vec_Pos_Inf:
            case util_Vec_Pos_Pos_Inf:
              return (ScaleX * ScaleY)/0.0;
            case util_Vec_Pos_Neg_Inf:
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              return 0.0/0.0;
            case util_Vec_Pos_Big:
              return (N - 1) * ScaleX * ScaleY * small + ScaleX * ScaleY * big;
            case util_Vec_Pos_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small + ScaleX * ScaleY * big) + ScaleX * ScaleY * big;
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * ScaleX * ScaleY * small;
            case util_Vec_Sine:
              return ScaleX * ScaleY - ScaleX * ScaleY;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
              exit(125);
          }
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
          switch(FillY){
            case util_Vec_Constant:
            case util_Vec_Pos_Inf:
            case util_Vec_Pos_Pos_Inf:
              return (ScaleX * ScaleY)/0.0;
            case util_Vec_Pos_Neg_Inf:
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              return 0.0/0.0;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
              exit(125);
          }
        case util_Vec_Pos_Neg_Inf:
          switch(FillY){
            case util_Vec_Constant:
            case util_Vec_Pos_Inf:
            case util_Vec_Pos_Pos_Inf:
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              return 0.0/0.0;
            case util_Vec_Pos_Neg_Inf:
              return (ScaleX * ScaleY)/0.0;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
              exit(125);
          }
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          return 0.0/0.0;
        case util_Vec_Pos_Big:
          switch(FillY){
            case util_Vec_Constant:
              return (N - 1) * ScaleX * ScaleY * small + ScaleX * ScaleY * big;
            case util_Vec_Pos_Big:
              return (N - 1) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * small) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Neg_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small - ScaleX * ScaleY * big * small) + ScaleX * ScaleY * big * big;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
              exit(125);
          }
        case util_Vec_Pos_Pos_Big:
          switch(FillY){
            case util_Vec_Constant:
              return ((N - 2) * ScaleX * ScaleY * small + ScaleX * ScaleY * big) + ScaleX * ScaleY * big;
            case util_Vec_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * small) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * big) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * ScaleX * ScaleY * small * small;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
              exit(125);
          }
        case util_Vec_Pos_Neg_Big:
          switch(FillY){
            case util_Vec_Constant:
              return (N - 2) * ScaleX * ScaleY * small;
            case util_Vec_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small - ScaleX * ScaleY * big * small) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * ScaleX * ScaleY * small * small;
            case util_Vec_Pos_Neg_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * big) + ScaleX * ScaleY * big * big;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
              exit(125);
          }
        case util_Vec_Sine:
          switch(FillY){
            case util_Vec_Constant:
              return ScaleX * ScaleY - ScaleX * ScaleY;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
              exit(125);
          }
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_saugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
          exit(125);
      }
  }
}

float wrap_saugsum_bound(int fold, int N, wrap_saugsum_func_t func, float *X, int incX, float *Y, int incY, float res, float ref){
  switch(func){
    case wrap_saugsum_RSSUM:
    case wrap_saugsum_SISIADD:
    case wrap_saugsum_SISADD:
    case wrap_saugsum_SISDEPOSIT:
    case wrap_saugsum_RSASUM:
      return sibound(fold, N, samax(N, X, incX));
    case wrap_saugsum_RSNRM2:
      {
        float amax = samax(N, X, incX);
        return sibound(fold, N, amax) * (amax / (res + ref));
      }
    case wrap_saugsum_RSDOT:
      return sibound(fold, N, samaxm(N, X, incX, Y, incY));
  }
}


#endif
