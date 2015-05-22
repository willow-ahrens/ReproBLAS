#ifndef DAUGSUM_WRAPPER_H
#define DAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>

#include "../common/test_util.h"

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
    double Scale = didnrm(fold, N, x, incx, ires);
    double res = ddiconv(fold, ires);
    free(ires);
    return Scale * sqrt(res);
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

double wrap_daugsum_result(int N, wrap_daugsum_func_t func, util_vec_fill_t fillX, double ScaleX, double CondX, util_vec_fill_t fillY, double ScaleY, double CondY){
  double small = 1.0 / (1024.0 * 1024.0 * 128.0); // 2^-27
  double big   = 1024.0 * 1024.0 * 128.0;         // 2^27
  switch(func){
    case wrap_daugsum_RDSUM:
    case wrap_daugsum_DIDIADD:
    case wrap_daugsum_DIDADD:
    case wrap_daugsum_DIDDEPOSIT:
      switch(fillX){
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
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX);
          exit(125);
      }
    case wrap_daugsum_RDASUM:
      switch(fillX){
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
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX);
          exit(125);
      }
    case wrap_daugsum_RDNRM2:
      switch(fillX){
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
        case util_Vec_Sine:
          //The sum of sin^2 on evenly spaced intervals over the range [0, 2*pi) is N/2
          return sqrt(N / 2.0) * fabs(ScaleX);
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX);
          exit(125);
      }
    case wrap_daugsum_RDDOT:
      switch(fillX){
        case util_Vec_Constant:
          switch(fillY){
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX, util_vec_fill_descs[fillY], ScaleY);
              exit(125);
          }
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
          switch(fillY){
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX, util_vec_fill_descs[fillY], ScaleY);
              exit(125);
          }
        case util_Vec_Pos_Neg_Inf:
          switch(fillY){
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX, util_vec_fill_descs[fillY], ScaleY);
              exit(125);
          }
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          return 0.0/0.0;
        case util_Vec_Pos_Big:
          switch(fillY){
            case util_Vec_Constant:
              return (N - 1) * ScaleX * ScaleY * small + ScaleX * ScaleY * big;
            case util_Vec_Pos_Big:
              return (N - 1) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * small) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Neg_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small - ScaleX * ScaleY * big * small) + ScaleX * ScaleY * big * big;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX, util_vec_fill_descs[fillY], ScaleY);
              exit(125);
          }
        case util_Vec_Pos_Pos_Big:
          switch(fillY){
            case util_Vec_Constant:
              return ((N - 2) * ScaleX * ScaleY * small + ScaleX * ScaleY * big) + ScaleX * ScaleY * big;
            case util_Vec_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * small) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * big) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * ScaleX * ScaleY * small * small;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX, util_vec_fill_descs[fillY], ScaleY);
              exit(125);
          }
        case util_Vec_Pos_Neg_Big:
          switch(fillY){
            case util_Vec_Constant:
              return (N - 2) * ScaleX * ScaleY * small;
            case util_Vec_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small - ScaleX * ScaleY * big * small) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * ScaleX * ScaleY * small * small;
            case util_Vec_Pos_Neg_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * big) + ScaleX * ScaleY * big * big;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX, util_vec_fill_descs[fillY], ScaleY);
              exit(125);
          }
        case util_Vec_Sine:
          switch(fillY){
            case util_Vec_Constant:
              return ScaleX * ScaleY - ScaleX * ScaleY;
            case util_Vec_Sine:
              //The sum of sin^2 on evenly spaced intervals over the range [0, 2*pi) is N/2
              return (N / 2.0) * ScaleX * ScaleY;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX, util_vec_fill_descs[fillY], ScaleY);
              exit(125);
          }
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[fillX], ScaleX, util_vec_fill_descs[fillY], ScaleY);
          exit(125);
      }
  }
}

#endif
