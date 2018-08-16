#ifndef DAUGSUM_WRAPPER_H
#define DAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <binnedBLAS.h>
#include <binned.h>
#include "../../config.h"

#include "../common/test_util.h"

typedef enum wrap_daugsum_func {
  wrap_daugsum_RDSUM = 0,
  wrap_daugsum_RDASUM,
  wrap_daugsum_RDNRM2,
  wrap_daugsum_RDDOT,
  wrap_daugsum_DBDBADD,
  wrap_daugsum_DIDADD,
  wrap_daugsum_DIDDEPOSIT
} wrap_daugsum_func_t;

typedef double (*wrap_daugsum)(int, int, double*, int, double*, int);
typedef void (*wrap_diaugsum)(int, int, double*, int, double*, int, double_binned*);
static const int wrap_daugsum_func_n_names = 7;
static const char* wrap_daugsum_func_names[] = {"rdsum",
                                                "rdasum",
                                                "rdnrm2",
                                                "rddot",
                                                "dbdbadd",
                                                "dbdadd",
                                                "dbddeposit"};
static const char* wrap_daugsum_func_descs[] = {"rdsum",
                                                "rdasum",
                                                "rdnrm2",
                                                "rddot",
                                                "dbdbadd",
                                                "dbdadd",
                                                "dbddeposit"};

double wrap_rdsum(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DIDEFAULTFOLD){
    return reproBLAS_dsum(N, x, incx);
  }else{
    return reproBLAS_rdsum(fold, N, x, incx);
  }
}

void wrap_dbdsum(int fold, int N, double *x, int incx, double *y, int incy, double_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_dbdsum(fold, N, x, incx, z);
}

double wrap_rdasum(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DIDEFAULTFOLD){
    return reproBLAS_dasum(N, x, incx);
  }else{
    return reproBLAS_rdasum(fold, N, x, incx);
  }
}

void wrap_dbdasum(int fold, int N, double *x, int incx, double *y, int incy, double_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_dbdasum(fold, N, x, incx, z);
}

double wrap_rdnrm2(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DIDEFAULTFOLD){
    return reproBLAS_dnrm2(N, x, incx);
  }else{
    return reproBLAS_rdnrm2(fold, N, x, incx);
  }
}

void wrap_dbdssq(int fold, int N, double *x, int incx, double *y, int incy, double_binned *z) {
  (void)y;
  (void)incy;
  binnedBLAS_dbdssq(fold, N, x, incx, 0.0, z);
}

double wrap_rddot(int fold, int N, double *x, int incx, double *y, int incy) {
  if(fold == DIDEFAULTFOLD){
    return reproBLAS_ddot(N, x, incx, y, incy);
  }else{
    return reproBLAS_rddot(fold, N, x, incx, y, incy);
  }
}

void wrap_dbddot(int fold, int N, double *x, int incx, double *y, int incy, double_binned *z) {
  binnedBLAS_dbddot(fold, N, x, incx, y, incy, z);
}

double wrap_rdbdbadd(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  double_binned *ires = binned_dballoc(fold);
  double_binned *itmp = binned_dballoc(fold);
  binned_dbsetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    binned_dbdconv(fold, x[i * incx], itmp);
    binned_dbdbadd(fold, itmp, ires);
  }
  double res = binned_ddbconv(fold, ires);
  free(ires);
  free(itmp);
  return res;
}

void wrap_dbdbadd(int fold, int N, double *x, int incx, double *y, int incy, double_binned *z) {
  (void)y;
  (void)incy;
  double_binned *itmp = binned_dballoc(fold);
  int i;
  for(i = 0; i < N; i++){
    binned_dbdconv(fold, x[i * incx], itmp);
    binned_dbdbadd(fold, itmp, z);
  }
  free(itmp);
}

double wrap_rdbdadd(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  double_binned *ires = binned_dballoc(fold);
  binned_dbsetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    binned_dbdadd(fold, x[i * incx], ires);
  }
  double res = binned_ddbconv(fold, ires);
  free(ires);
  return res;
}

void wrap_dbdadd(int fold, int N, double *x, int incx, double *y, int incy, double_binned *z) {
  (void)y;
  (void)incy;
  int i;
  for(i = 0; i < N; i++){
    binned_dbdadd(fold, x[i * incx], z);
  }
}

double wrap_rdbddeposit(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  double_binned *ires = binned_dballoc(fold);
  binned_dbsetzero(fold, ires);
  double amax = binnedBLAS_damax(N, x, incx);
  binned_dbdupdate(fold, amax, ires);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= binned_DBENDURANCE){
      binned_dbrenorm(fold, ires);
      j = 0;
    }
    binned_dbddeposit(fold, x[i * incx], ires);
    j++;
  }
  binned_dbrenorm(fold, ires);
  double res = binned_ddbconv(fold, ires);
  free(ires);
  return res;
}

void wrap_dbddeposit(int fold, int N, double *x, int incx, double *y, int incy, double_binned *z) {
  (void)y;
  (void)incy;
  double amax = binnedBLAS_damax(N, x, incx);
  binned_dbdupdate(fold, amax, z);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= binned_DBENDURANCE){
      binned_dbrenorm(fold, z);
      j = 0;
    }
    binned_dbddeposit(fold, x[i * incx], z);
    j++;
  }
  binned_dbrenorm(fold, z);
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
    case wrap_daugsum_DBDBADD:
      return wrap_rdbdbadd;
    case wrap_daugsum_DIDADD:
      return wrap_rdbdadd;
    case wrap_daugsum_DIDDEPOSIT:
      return wrap_rdbddeposit;
  }
  return NULL;
}

wrap_diaugsum wrap_diaugsum_func(wrap_daugsum_func_t func) {
  switch(func){
    case wrap_daugsum_RDSUM:
      return wrap_dbdsum;
    case wrap_daugsum_RDASUM:
      return wrap_dbdasum;
    case wrap_daugsum_RDNRM2:
      return wrap_dbdssq;
    case wrap_daugsum_RDDOT:
      return wrap_dbddot;
    case wrap_daugsum_DBDBADD:
      return wrap_dbdbadd;
    case wrap_daugsum_DIDADD:
      return wrap_dbdadd;
    case wrap_daugsum_DIDDEPOSIT:
      return wrap_dbddeposit;
  }
  return NULL;
}

double wrap_daugsum_result(int N, wrap_daugsum_func_t func, util_vec_fill_t FillX, double RealScaleX, double ImagScaleX, util_vec_fill_t FillY, double RealScaleY, double ImagScaleY){
  double small = 1.0 / (1024.0 * 1024.0 * 128.0); // 2^-27
  double big   = 1024.0 * 1024.0 * 128.0;         // 2^27
  switch(func){
    case wrap_daugsum_RDSUM:
    case wrap_daugsum_DBDBADD:
    case wrap_daugsum_DIDADD:
    case wrap_daugsum_DIDDEPOSIT:
      switch(FillX){
        case util_Vec_Constant:
          return N * RealScaleX;
        case util_Vec_Mountain:
          return 0;
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
          return RealScaleX * INFINITY;
        case util_Vec_Pos_Neg_Inf:
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          return NAN;
        case util_Vec_Pos_Big:
          return (N - 1) * (RealScaleX * small) + RealScaleX * big;
        case util_Vec_Pos_Pos_Big:
          return (N - 2) * (RealScaleX * small) + (RealScaleX * big + RealScaleX * big);
        case util_Vec_Pos_Neg_Big:
          return (N - 2) * (RealScaleX * small);
        case util_Vec_Sine:
          return 0.0;
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX);
          exit(125);
      }

    case wrap_daugsum_RDASUM:
      switch(FillX){
        case util_Vec_Constant:
          return N * fabs(RealScaleX);
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
        case util_Vec_Pos_Neg_Inf:
          return fabs(RealScaleX) * INFINITY;
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          return NAN;
        case util_Vec_Pos_Big:
          return (N - 1) * (fabs(RealScaleX) * small) + fabs(RealScaleX) * big;
        case util_Vec_Pos_Pos_Big:
        case util_Vec_Pos_Neg_Big:
          return (N - 2) * (fabs(RealScaleX) * small) + (fabs(RealScaleX) * big + fabs(RealScaleX) * big);
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX);
          exit(125);
      }

    case wrap_daugsum_RDNRM2:
      {
        double new_scale;
        switch(FillX){
          case util_Vec_Constant:
            new_scale = binned_dscale(RealScaleX);
            RealScaleX /= new_scale;
            return sqrt(N * (RealScaleX * RealScaleX)) * new_scale;
          case util_Vec_Pos_Inf:
          case util_Vec_Pos_Pos_Inf:
          case util_Vec_Pos_Neg_Inf:
            return fabs(RealScaleX) * INFINITY;
          case util_Vec_NaN:
          case util_Vec_Pos_Inf_NaN:
          case util_Vec_Pos_Pos_Inf_NaN:
          case util_Vec_Pos_Neg_Inf_NaN:
            return NAN;
          case util_Vec_Pos_Big:
            new_scale = binned_dscale(RealScaleX * big);
            small *= RealScaleX;
            small /= new_scale;
            big *= RealScaleX;
            big /= new_scale;
            return sqrt((N - 1) * (small * small) + big * big) * new_scale;
          case util_Vec_Pos_Pos_Big:
          case util_Vec_Pos_Neg_Big:
            new_scale = binned_dscale(RealScaleX * big);
            small *= RealScaleX;
            small /= new_scale;
            big *= RealScaleX;
            big /= new_scale;
            return sqrt((N - 2) * (small * small) + (big * big + big * big)) * new_scale;
          default:
            fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX);
            exit(125);
        }
      }

    case wrap_daugsum_RDDOT:
      switch(FillX){
        case util_Vec_Mountain:
          switch(FillY){
            case util_Vec_Constant:
              return 0;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
              exit(125);
          }
        case util_Vec_Constant:
          switch(FillY){
            case util_Vec_Constant:
              return N * (RealScaleX * RealScaleY);
            case util_Vec_Mountain:
              return 0;
            case util_Vec_Pos_Inf:
            case util_Vec_Pos_Pos_Inf:
              return (RealScaleX * RealScaleY) * INFINITY;
            case util_Vec_Pos_Neg_Inf:
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              return NAN;
            case util_Vec_Pos_Big:
              return (N - 1) * (RealScaleX * RealScaleY * small) + RealScaleX * RealScaleY * big;
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small) + (RealScaleX * RealScaleY * big + RealScaleX * RealScaleY * big);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small);
            case util_Vec_Sine:
              return 0.0;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
              exit(125);
          }
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
          switch(FillY){
            case util_Vec_Constant:
            case util_Vec_Pos_Inf:
            case util_Vec_Pos_Pos_Inf:
              return (RealScaleX * RealScaleY) * INFINITY;
            case util_Vec_Pos_Neg_Inf:
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              return NAN;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
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
              return NAN;
            case util_Vec_Pos_Neg_Inf:
              return (RealScaleX * RealScaleY) * INFINITY;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
              exit(125);
          }
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          return NAN;
        case util_Vec_Pos_Big:
          switch(FillY){
            case util_Vec_Constant:
              return (N - 1) * (RealScaleX * RealScaleY * small) + RealScaleX * RealScaleY * big;
            case util_Vec_Pos_Big:
              return (N - 1) * (RealScaleX * RealScaleY * small * small) + RealScaleX * RealScaleY * big * big;
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small * small) + (RealScaleX * RealScaleY * big * small + RealScaleX * RealScaleY * big * big);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small * small) - (RealScaleX * RealScaleY * big * small - RealScaleX * RealScaleY * big * big);
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
              exit(125);
          }
        case util_Vec_Pos_Pos_Big:
          switch(FillY){
            case util_Vec_Constant:
              return (N - 2) * (RealScaleX * RealScaleY * small) + (RealScaleX * RealScaleY * big + RealScaleX * RealScaleY * big);
            case util_Vec_Pos_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small * small) + (RealScaleX * RealScaleY * big * small + RealScaleX * RealScaleY * big * big);
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small * small) + (RealScaleX * RealScaleY * big * big + RealScaleX * RealScaleY * big * big);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small * small);
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
              exit(125);
          }
        case util_Vec_Pos_Neg_Big:
          switch(FillY){
            case util_Vec_Constant:
              return (N - 2) * (RealScaleX * RealScaleY * small);
            case util_Vec_Pos_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small * small) - (RealScaleX * RealScaleY * big * small - RealScaleX * RealScaleY * big * big);
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small * small);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * (RealScaleX * RealScaleY * small * small) + (RealScaleX * RealScaleY * big * big + RealScaleX * RealScaleY * big * big);
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
              exit(125);
          }
        case util_Vec_Sine:
          switch(FillY){
            case util_Vec_Constant:
              return 0.0;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
              exit(125);
          }
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
          exit(125);
      }
    default:
      fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_daugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, util_vec_fill_descs[FillY], RealScaleY);
      exit(125);
  }
}

double wrap_daugsum_bound(int fold, int N, wrap_daugsum_func_t func, double *X, int incX, double *Y, int incY, double res, double ref){
  switch(func){
    case wrap_daugsum_RDSUM:
    case wrap_daugsum_DBDBADD:
    case wrap_daugsum_DIDADD:
    case wrap_daugsum_DIDDEPOSIT:
    case wrap_daugsum_RDASUM:
      return binned_dbbound(fold, N, binnedBLAS_damax(N, X, incX), res);
    case wrap_daugsum_RDNRM2:
      {
        double amax = binnedBLAS_damax(N, X, incX);
        double scale = binned_dscale(amax);
        if (amax == 0.0){
          return 0.0;
        }
        return binned_dbbound(fold, N, (amax/scale) * (amax/scale), res) * (scale / (res + ref)) * scale;
      }
    case wrap_daugsum_RDDOT:
      return binned_dbbound(fold, N, binnedBLAS_damaxm(N, X, incX, Y, incY), res);
  }
  fprintf(stderr, "ReproBLAS error: unknown bound for %s\n", wrap_daugsum_func_descs[func]);
  exit(125);
}

#endif
