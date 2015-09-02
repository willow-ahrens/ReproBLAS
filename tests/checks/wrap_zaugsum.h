#ifndef ZAUGSUM_WRAPPER_H
#define ZAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <idxdBLAS.h>
#include <idxd.h>
#include "../../config.h"

#include "../common/test_util.h"

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
  if(fold == DIDEFAULTFOLD){
    double complex res;
    reproBLAS_zsum_sub(N, x, incx, &res);
    return res;
  }else{
    double complex res;
    reproBLAS_rzsum_sub(fold, N, x, incx, &res);
    return res;
  }
}

void wrap_zizsum(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  idxdBLAS_zizsum(fold, N, x, incx, z);
}

double complex wrap_rdzasum(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DIDEFAULTFOLD){
    return reproBLAS_dzasum(N, x, incx);
  }else{
    return reproBLAS_rdzasum(fold, N, x, incx);
  }
}

void wrap_dizasum(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  idxdBLAS_dmzasum(fold, N, x, incx, z, 2, z + 2 * fold, 2);
}

double complex wrap_rdznrm2(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == DIDEFAULTFOLD){
    return reproBLAS_dznrm2(N, x, incx);
  }else{
    return reproBLAS_rdznrm2(fold, N, x, incx);
  }
}

void wrap_dizssq(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  idxdBLAS_dmzssq(fold, N, x, incx, 0.0, z, 2, z + 2 * fold, 2);
}

double complex wrap_rzdotu(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  if(fold == DIDEFAULTFOLD){
    double complex res;
    reproBLAS_zdotu_sub(N, x, incx, y, incy, &res);
    return res;
  }else{
    double complex res;
    reproBLAS_rzdotu_sub(fold, N, x, incx, y, incy, &res);
    return res;
  }
}

void wrap_zizdotu(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  idxdBLAS_zizdotu(fold, N, x, incx, y, incy, z);
}

double complex wrap_rzdotc(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  if(fold == DIDEFAULTFOLD){
    double complex res;
    reproBLAS_rzdotc_sub(fold, N, x, incx, y, incy, &res);
    return res;
  }else{
    double complex res;
    reproBLAS_rzdotc_sub(fold, N, x, incx, y, incy, &res);
    return res;
  }
}

void wrap_zizdotc(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  idxdBLAS_zizdotc(fold, N, x, incx, y, incy, z);
}

double complex wrap_rziziadd(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  double_complex_indexed *ires = idxd_zialloc(fold);
  double_complex_indexed *itmp = idxd_zialloc(fold);
  idxd_zisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    idxd_zizconv(fold, x + i * incx, itmp);
    idxd_ziziadd(fold, itmp, ires);
  }
  double complex res;
  idxd_zziconv_sub(fold, ires, &res);
  free(ires);
  free(itmp);
  return res;
}

void wrap_ziziadd(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  double_complex_indexed *itmp = idxd_zialloc(fold);
  int i;
  for(i = 0; i < N; i++){
    idxd_zizconv(fold, x + i * incx, itmp);
    idxd_ziziadd(fold, itmp, z);
  }
  free(itmp);
}

double complex wrap_rzizadd(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  double_complex_indexed *ires = idxd_zialloc(fold);
  idxd_zisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    idxd_zizadd(fold, x + i * incx, ires);
  }
  double complex res;
  idxd_zziconv_sub(fold, ires, &res);
  free(ires);
  return res;
}

void wrap_zizadd(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  int i;
  for(i = 0; i < N; i++){
    idxd_zizadd(fold, x + i * incx, z);
  }
}

double complex wrap_rzizdeposit(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  (void)y;
  (void)incy;
  double_complex_indexed *ires = idxd_zialloc(fold);
  idxd_zisetzero(fold, ires);
  double complex amax;
  idxdBLAS_zamax_sub(N, x, incx, &amax);
  idxd_zizupdate(fold, &amax, ires);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= idxd_DIENDURANCE){
      idxd_zirenorm(fold, ires);
      j = 0;
    }
    idxd_zizdeposit(fold, x + i * incx, ires);
    j++;
  }
  idxd_zirenorm(fold, ires);
  double complex res;
  idxd_zziconv_sub(fold, ires, &res);
  free(ires);
  return res;
}

void wrap_zizdeposit(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  double complex amax;
  idxdBLAS_zamax_sub(N, x, incx, &amax);
  idxd_zizupdate(fold, &amax, z);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= idxd_DIENDURANCE){
      idxd_zirenorm(fold, z);
      j = 0;
    }
    idxd_zizdeposit(fold, x + i * incx, z);
    j++;
  }
  idxd_zirenorm(fold, z);
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
      return wrap_dizssq;
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

double complex wrap_zaugsum_result(int N, wrap_zaugsum_func_t func, util_vec_fill_t FillX, double RealScaleX, double ImagScaleX, util_vec_fill_t FillY, double RealScaleY, double ImagScaleY){
  double small = 1.0 / (1024.0 * 1024.0 * 128.0); // 2^-27
  double big   = 1024.0 * 1024.0 * 128.0;         // 2^27
  double complex ScaleX     = RealScaleX + ImagScaleX * I;
  double complex ScaleY     = RealScaleY + ImagScaleY * I;
  double complex tmpX0;
  double *tmpX0_base = (double*)&tmpX0;
  double complex tmpX1;
  double *tmpX1_base = (double*)&tmpX1;
  double complex tmpY0;
  double *tmpY0_base = (double*)&tmpY0;
  double complex tmpY1;
  double *tmpY1_base = (double*)&tmpY1;
  switch(func){
    case wrap_zaugsum_RZSUM:
    case wrap_zaugsum_ZIZIADD:
    case wrap_zaugsum_ZIZADD:
    case wrap_zaugsum_ZIZDEPOSIT:
      switch(FillX){
        case util_Vec_Constant:
          return N * ScaleX;
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
        case util_Vec_Pos_Neg_Inf:
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          switch(FillX){
            case util_Vec_Pos_Inf:
            case util_Vec_Pos_Pos_Inf:
              tmpX0_base[0] = INFINITY * RealScaleX;
              tmpX0_base[1] = INFINITY * ImagScaleX;
              break;
            case util_Vec_Pos_Neg_Inf:
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              tmpX0_base[0] = NAN;
              tmpX0_base[1] = NAN;
              break;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
          if(RealScaleX == 0.0){
            tmpX0_base[0] = 0.0;
          }
          if(ImagScaleX == 0.0){
            tmpX0_base[1] = 0.0;
          }
          return tmpX0;
        case util_Vec_Pos_Big:
          return (N - 1) * ScaleX * small + ScaleX * big;
        case util_Vec_Pos_Pos_Big:
          return (N - 2) * ScaleX * small + (ScaleX * big + ScaleX * big);
        case util_Vec_Pos_Neg_Big:
          return (N - 2) * ScaleX * small;
        case util_Vec_Sine:
          return 0.0;
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX);
          exit(125);
      }

    case wrap_zaugsum_RDZASUM:
      switch(FillX){
        case util_Vec_Constant:
          return N * (fabs(RealScaleX) + fabs(ImagScaleX));
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
        case util_Vec_Pos_Neg_Inf:
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          switch(FillX){
            case util_Vec_Pos_Inf:
            case util_Vec_Pos_Pos_Inf:
            case util_Vec_Pos_Neg_Inf:
              tmpX0_base[0] = INFINITY;
              tmpX0_base[1] = INFINITY;
              break;
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              tmpX0_base[0] = NAN;
              tmpX0_base[1] = NAN;
              break;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
          if(RealScaleX == 0.0){
            tmpX0_base[0] = 0.0;
          }
          if(ImagScaleX == 0.0){
            tmpX0_base[1] = 0.0;
          }
          return tmpX0_base[0] + tmpX0_base[1];
        case util_Vec_Pos_Big:
          return (N - 1) * (fabs(RealScaleX) + fabs(ImagScaleX)) * small + (fabs(RealScaleX) + fabs(ImagScaleX)) * big;
        case util_Vec_Pos_Pos_Big:
        case util_Vec_Pos_Neg_Big:
          return (N - 2) * (fabs(RealScaleX) + fabs(ImagScaleX)) * small + ((fabs(RealScaleX) + fabs(ImagScaleX)) * big + (fabs(RealScaleX) + fabs(ImagScaleX)) * big);
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX);
          exit(125);
      }

    case wrap_zaugsum_RDZNRM2:
      {
        double new_scale = MAX(fabs(RealScaleX), fabs(ImagScaleX));
        double big_real;
        double big_imag;
        double small_real;
        double small_imag;
        switch(FillX){
          case util_Vec_Constant:
            new_scale = idxd_dscale(new_scale);
            RealScaleX /= new_scale;
            ImagScaleX /= new_scale;
            return sqrt(N * (RealScaleX * RealScaleX + ImagScaleX * ImagScaleX)) * new_scale;
          case util_Vec_Pos_Inf:
          case util_Vec_Pos_Pos_Inf:
          case util_Vec_Pos_Neg_Inf:
          case util_Vec_NaN:
          case util_Vec_Pos_Inf_NaN:
          case util_Vec_Pos_Pos_Inf_NaN:
          case util_Vec_Pos_Neg_Inf_NaN:
            switch(FillX){
              case util_Vec_Pos_Inf:
              case util_Vec_Pos_Pos_Inf:
              case util_Vec_Pos_Neg_Inf:
                tmpX0_base[0] = INFINITY;
                tmpX0_base[1] = INFINITY;
                break;
              case util_Vec_NaN:
              case util_Vec_Pos_Inf_NaN:
              case util_Vec_Pos_Pos_Inf_NaN:
              case util_Vec_Pos_Neg_Inf_NaN:
                tmpX0_base[0] = NAN;
                tmpX0_base[1] = NAN;
                break;
              default:
                fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
                exit(125);
            }
            if(RealScaleX == 0.0){
              tmpX0_base[0] = 0.0;
            }
            if(ImagScaleX == 0.0){
              tmpX0_base[1] = 0.0;
            }
            return tmpX0_base[0] + tmpX0_base[1];
          case util_Vec_Pos_Big:
            new_scale = idxd_dscale(new_scale * big);
            big_real = (big * RealScaleX)/new_scale;
            big_imag = (big * ImagScaleX)/new_scale;
            small_real = (small * RealScaleX)/new_scale;
            small_imag = (small * ImagScaleX)/new_scale;
            return sqrt((N - 1) * (small_real * small_real + small_imag * small_imag) + big_real * big_real + big_imag * big_imag) * new_scale;
          case util_Vec_Pos_Pos_Big:
          case util_Vec_Pos_Neg_Big:
            new_scale = idxd_dscale(new_scale * big);
            big_real = (big * RealScaleX)/new_scale;
            big_imag = (big * ImagScaleX)/new_scale;
            small_real = (small * RealScaleX)/new_scale;
            small_imag = (small * ImagScaleX)/new_scale;
            return sqrt((N - 2) * (small_real * small_real + small_imag * small_imag) + (big_real * big_real + big_imag * big_imag + big_real * big_real + big_imag * big_imag)) * new_scale;
          default:
            fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX);
            exit(125);
        }
      }

    case wrap_zaugsum_RZDOTC:
      ScaleX = RealScaleX - ImagScaleX * I;
      ImagScaleX = -1 * ImagScaleX;

    case wrap_zaugsum_RZDOTU:
      switch(FillX){
        case util_Vec_Constant:
          switch(FillY){
            case util_Vec_Constant:
              return N * (ScaleX * ScaleY);
            case util_Vec_Pos_Inf:
            case util_Vec_Pos_Pos_Inf:
            case util_Vec_Pos_Neg_Inf:
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              switch(FillY){
                case util_Vec_Pos_Inf:
                case util_Vec_Pos_Pos_Inf:
                  tmpY0_base[0] = INFINITY * RealScaleY;
                  tmpY0_base[1] = INFINITY * ImagScaleY;
                  break;
                case util_Vec_Pos_Neg_Inf:
                case util_Vec_NaN:
                case util_Vec_Pos_Inf_NaN:
                case util_Vec_Pos_Pos_Inf_NaN:
                case util_Vec_Pos_Neg_Inf_NaN:
                  tmpY0_base[0] = NAN;
                  tmpY0_base[1] = NAN;
                  break;
                default:
                  fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
                  exit(125);
              }
              if(RealScaleY == 0.0){
                tmpY0_base[0] = 0.0;
              }
              if(ImagScaleY == 0.0){
                tmpY0_base[1] = 0.0;
              }
              return zmul(tmpY0, ScaleX);
            case util_Vec_Pos_Big:
              return (N - 1) * (ScaleX * ScaleY * small) + ScaleX * ScaleY * big;
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * (ScaleX * ScaleY * small) + (ScaleX * ScaleY * big + ScaleX * ScaleY * big);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * (ScaleX * ScaleY * small);
            case util_Vec_Sine:
              return 0.0;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
        case util_Vec_Pos_Neg_Inf:
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
          switch(FillX){
            case util_Vec_Pos_Inf:
              tmpX0_base[0] = INFINITY * RealScaleX;
              tmpX0_base[1] = INFINITY * ImagScaleX;
              tmpX1_base[0] = 0.0;
              tmpX1_base[1] = 0.0;
            case util_Vec_Pos_Pos_Inf:
              tmpX0_base[0] = INFINITY * RealScaleX;
              tmpX0_base[1] = INFINITY * ImagScaleX;
              tmpX1_base[0] = INFINITY * RealScaleX;
              tmpX1_base[1] = INFINITY * ImagScaleX;
              break;
            case util_Vec_Pos_Neg_Inf:
              tmpX0_base[0] = INFINITY * RealScaleX;
              tmpX0_base[1] = INFINITY * ImagScaleX;
              tmpX1_base[0] = -1.0 * INFINITY * RealScaleX;
              tmpX1_base[1] = -1.0 * INFINITY * ImagScaleX;
              break;
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              tmpX0_base[0] = NAN;
              tmpX0_base[1] = NAN;
              tmpX1_base[0] = 0.0;
              tmpX1_base[1] = 0.0;
              break;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
          if(RealScaleX == 0.0){
            tmpX0_base[0] = 0.0;
            tmpX1_base[0] = 0.0;
          }
          if(ImagScaleX == 0.0){
            tmpX0_base[1] = 0.0;
            tmpX1_base[1] = 0.0;
          }
          switch(FillY){
            case util_Vec_Constant:
              return zmul(tmpX0, ScaleY) + zmul(tmpX1, ScaleY);
            case util_Vec_Pos_Inf:
            case util_Vec_Pos_Pos_Inf:
            case util_Vec_Pos_Neg_Inf:
            case util_Vec_NaN:
            case util_Vec_Pos_Inf_NaN:
            case util_Vec_Pos_Pos_Inf_NaN:
            case util_Vec_Pos_Neg_Inf_NaN:
              switch(FillY){
                case util_Vec_Pos_Inf:
                  tmpY0_base[0] = INFINITY * RealScaleY;
                  tmpY0_base[1] = INFINITY * ImagScaleY;
                  tmpY1_base[0] = 0.0;
                  tmpY1_base[1] = 0.0;
                case util_Vec_Pos_Pos_Inf:
                  tmpY0_base[0] = INFINITY * RealScaleY;
                  tmpY0_base[1] = INFINITY * ImagScaleY;
                  tmpY1_base[0] = INFINITY * RealScaleY;
                  tmpY1_base[1] = INFINITY * ImagScaleY;
                  break;
                case util_Vec_Pos_Neg_Inf:
                  tmpY0_base[0] = INFINITY * RealScaleY;
                  tmpY0_base[1] = INFINITY * ImagScaleY;
                  tmpY1_base[0] = -1.0 * INFINITY * RealScaleY;
                  tmpY1_base[1] = -1.0 * INFINITY * ImagScaleY;
                  break;
                case util_Vec_NaN:
                case util_Vec_Pos_Inf_NaN:
                case util_Vec_Pos_Pos_Inf_NaN:
                case util_Vec_Pos_Neg_Inf_NaN:
                  tmpY0_base[0] = NAN;
                  tmpY0_base[1] = NAN;
                  tmpY1_base[0] = 0.0;
                  tmpY1_base[1] = 0.0;
                  break;
                default:
                  fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
                  exit(125);
              }
              if(RealScaleY == 0.0){
                tmpY0_base[0] = 0.0;
                tmpY1_base[0] = 0.0;
              }
              if(ImagScaleY == 0.0){
                tmpY0_base[1] = 0.0;
                tmpY1_base[1] = 0.0;
              }
              return zmul(tmpX0, tmpY0) + zmul(tmpX1, tmpY1);
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
        case util_Vec_Pos_Big:
          switch(FillY){
            case util_Vec_Constant:
              return (N - 1) * (ScaleX * ScaleY * small) + ScaleX * ScaleY * big;
            case util_Vec_Pos_Big:
              return (N - 1) * (ScaleX * ScaleY * small * small) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * (ScaleX * ScaleY * small * small) + (ScaleX * ScaleY * big * small + ScaleX * ScaleY * big * big);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * (ScaleX * ScaleY * small * small) - (ScaleX * ScaleY * big * small - ScaleX * ScaleY * big * big);
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
        case util_Vec_Pos_Pos_Big:
          switch(FillY){
            case util_Vec_Constant:
              return (N - 2) * (ScaleX * ScaleY * small) + (ScaleX * ScaleY * big + ScaleX * ScaleY * big);
            case util_Vec_Pos_Big:
              return (N - 2) * (ScaleX * ScaleY * small * small) + (ScaleX * ScaleY * big * small + ScaleX * ScaleY * big * big);
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * (ScaleX * ScaleY * small * small) + (ScaleX * ScaleY * big * big + ScaleX * ScaleY * big * big);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * (ScaleX * ScaleY * small * small);
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
        case util_Vec_Pos_Neg_Big:
          switch(FillY){
            case util_Vec_Constant:
              return (N - 2) * (ScaleX * ScaleY * small);
            case util_Vec_Pos_Big:
              return (N - 2) * (ScaleX * ScaleY * small * small) - (ScaleX * ScaleY * big * small - ScaleX * ScaleY * big * big);
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * (ScaleX * ScaleY * small * small);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * (ScaleX * ScaleY * small * small) + (ScaleX * ScaleY * big * big + ScaleX * ScaleY * big * big);
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
        case util_Vec_Sine:
          switch(FillY){
            case util_Vec_Constant:
              return 0.0;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
          exit(125);
      }
    default:
      fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
      exit(125);
  }
}

double complex wrap_zaugsum_bound(int fold, int N, wrap_zaugsum_func_t func, double complex *X, int incX, double complex *Y, int incY, double complex res, double complex ref){
  switch(func){
    case wrap_zaugsum_RZSUM:
    case wrap_zaugsum_ZIZIADD:
    case wrap_zaugsum_ZIZADD:
    case wrap_zaugsum_ZIZDEPOSIT:
      {
        double complex amax;
        double complex bound;
        double *bound_base = (double*)&bound;
        idxdBLAS_zamax_sub(N, X, incX, &amax);
        bound_base[0] = idxd_dibound(fold, N, creal(amax), creal(res));
        bound_base[1] = idxd_dibound(fold, N, cimag(amax), cimag(res));
        return bound;
      }
    case wrap_zaugsum_RDZASUM:
      {
        double complex amax2;
        double amax;
        idxdBLAS_zamax_sub(N, X, incX, &amax2);
        amax = MAX(creal(amax2), cimag(amax2));
        return idxd_dibound(fold, N, amax, creal(res));
      }
    case wrap_zaugsum_RDZNRM2:
      {
        double complex amax2;
        double amax;
        double scale;
        idxdBLAS_zamax_sub(N, X, incX, &amax2);
        amax = MAX(creal(amax2), cimag(amax2));
        scale = idxd_dscale(amax);
        if (amax == 0.0){
          return 0.0;
        }
        return idxd_dibound(fold, N, (amax/scale) * (amax/scale), (creal(res)/scale) * (creal(res)/scale)) * (scale / (creal(res) + creal(ref))) * scale;
      }
    case wrap_zaugsum_RZDOTU:
    case wrap_zaugsum_RZDOTC:
      {
        double complex amaxm;
        double complex bound;
        double *bound_base = (double*)&bound;
        idxdBLAS_zamaxm_sub(N, X, incX, Y, incY, &amaxm);
        bound_base[0] = idxd_dibound(fold, N, creal(amaxm), creal(res));
        bound_base[1] = idxd_dibound(fold, N, cimag(amaxm), cimag(res));
        return bound;
      }
  }
  fprintf(stderr, "ReproBLAS error: unknown bound for %s\n", wrap_zaugsum_func_descs[func]);
  exit(125);
}

#endif
