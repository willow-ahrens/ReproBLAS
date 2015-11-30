#ifndef DAUGSUM_WRAPPER_H
#define DAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <idxdBLAS.h>
#include <idxd.h>
#include "../../config.h"

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
  if(fold == DIDEFAULTFOLD){
    return reproBLAS_dsum(N, x, incx);
  }else{
    return reproBLAS_rdsum(fold, N, x, incx);
  }
}

void wrap_didsum(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  idxdBLAS_didsum(fold, N, x, incx, z);
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

void wrap_didasum(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  idxdBLAS_didasum(fold, N, x, incx, z);
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

void wrap_didssq(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  idxdBLAS_didssq(fold, N, x, incx, 0.0, z);
}

double wrap_rddot(int fold, int N, double *x, int incx, double *y, int incy) {
  if(fold == DIDEFAULTFOLD){
    return reproBLAS_ddot(N, x, incx, y, incy);
  }else{
    return reproBLAS_rddot(fold, N, x, incx, y, incy);
  }
}

void wrap_diddot(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  idxdBLAS_diddot(fold, N, x, incx, y, incy, z);
}

double wrap_rdidiadd(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  double_indexed *ires = idxd_dialloc(fold);
  double_indexed *itmp = idxd_dialloc(fold);
  idxd_disetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    idxd_didconv(fold, x[i * incx], itmp);
    idxd_didiadd(fold, itmp, ires);
  }
  double res = idxd_ddiconv(fold, ires);
  free(ires);
  free(itmp);
  return res;
}

void wrap_didiadd(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  double_indexed *itmp = idxd_dialloc(fold);
  int i;
  for(i = 0; i < N; i++){
    idxd_didconv(fold, x[i * incx], itmp);
    idxd_didiadd(fold, itmp, z);
  }
  free(itmp);
}

double wrap_rdidadd(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  double_indexed *ires = idxd_dialloc(fold);
  idxd_disetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    idxd_didadd(fold, x[i * incx], ires);
  }
  double res = idxd_ddiconv(fold, ires);
  free(ires);
  return res;
}

void wrap_didadd(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  int i;
  for(i = 0; i < N; i++){
    idxd_didadd(fold, x[i * incx], z);
  }
}

double wrap_rdiddeposit(int fold, int N, double *x, int incx, double *y, int incy) {
  (void)y;
  (void)incy;
  double_indexed *ires = idxd_dialloc(fold);
  idxd_disetzero(fold, ires);
  double amax = idxdBLAS_damax(N, x, incx);
  idxd_didupdate(fold, amax, ires);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= idxd_DIENDURANCE){
      idxd_direnorm(fold, ires);
      j = 0;
    }
    idxd_diddeposit(fold, x[i * incx], ires);
    j++;
  }
  idxd_direnorm(fold, ires);
  double res = idxd_ddiconv(fold, ires);
  free(ires);
  return res;
}

void wrap_diddeposit(int fold, int N, double *x, int incx, double *y, int incy, double_indexed *z) {
  (void)y;
  (void)incy;
  double amax = idxdBLAS_damax(N, x, incx);
  idxd_didupdate(fold, amax, z);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= idxd_DIENDURANCE){
      idxd_direnorm(fold, z);
      j = 0;
    }
    idxd_diddeposit(fold, x[i * incx], z);
    j++;
  }
  idxd_direnorm(fold, z);
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
      return wrap_didssq;
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

double wrap_daugsum_result(int N, wrap_daugsum_func_t func, util_vec_fill_t FillX, double RealScaleX, double ImagScaleX, util_vec_fill_t FillY, double RealScaleY, double ImagScaleY){
  double small = 1.0 / (1024.0 * 1024.0 * 128.0); // 2^-27
  double big   = 1024.0 * 1024.0 * 128.0;         // 2^27
  switch(func){
    case wrap_daugsum_RDSUM:
    case wrap_daugsum_DIDIADD:
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
            new_scale = idxd_dscale(RealScaleX);
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
            new_scale = idxd_dscale(RealScaleX * big);
            small *= RealScaleX;
            small /= new_scale;
            big *= RealScaleX;
            big /= new_scale;
            return sqrt((N - 1) * (small * small) + big * big) * new_scale;
          case util_Vec_Pos_Pos_Big:
          case util_Vec_Pos_Neg_Big:
            new_scale = idxd_dscale(RealScaleX * big);
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
    case wrap_daugsum_DIDIADD:
    case wrap_daugsum_DIDADD:
    case wrap_daugsum_DIDDEPOSIT:
    case wrap_daugsum_RDASUM:
      return idxd_dibound(fold, N, idxdBLAS_damax(N, X, incX), res);
    case wrap_daugsum_RDNRM2:
      {
        double amax = idxdBLAS_damax(N, X, incX);
        double scale = idxd_dscale(amax);
        if (amax == 0.0){
          return 0.0;
        }
        return idxd_dibound(fold, N, (amax/scale) * (amax/scale), res) * (scale / (res + ref)) * scale;
      }
    case wrap_daugsum_RDDOT:
      return idxd_dibound(fold, N, idxdBLAS_damaxm(N, X, incX, Y, incY), res);
  }
  fprintf(stderr, "ReproBLAS error: unknown bound for %s\n", wrap_daugsum_func_descs[func]);
  exit(125);
}

#endif
