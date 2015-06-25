#ifndef ZAUGSUM_WRAPPER_H
#define ZAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>
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
  if(fold == DIDEFAULTFOLD){
    return rdzasum(N, x, incx);
  }else{
    double_indexed *ires = dialloc(fold);
    disetzero(fold, ires);
    dizasum(fold, N, x, incx, ires);
    double res = ddiconv(fold, ires);
    free(ires);
    return (double complex)res;
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
  if(fold == DIDEFAULTFOLD){
    return rdznrm2(N, x, incx);
  }else{
    double_indexed *ires = dialloc(fold);
    disetzero(fold, ires);
    double scale = dizssq(fold, N, x, incx, 0.0, ires);
    double res = ddiconv(fold, ires);
    free(ires);
    return (double complex)(scale * sqrt(res));
  }
}

void wrap_dizssq(int fold, int N, double complex *x, int incx, double complex *y, int incy, double_complex_indexed *z) {
  (void)y;
  (void)incy;
  dmzssq(fold, N, x, incx, 0.0, z, 2, z + 2 * fold, 2);
}

double complex wrap_rzdotu(int fold, int N, double complex *x, int incx, double complex *y, int incy) {
  if(fold == DIDEFAULTFOLD){
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
  if(fold == DIDEFAULTFOLD){
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
    if(j >= DIENDURANCE){
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
    if(j >= DIENDURANCE){
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
            new_scale = dscale(new_scale);
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
            new_scale = dscale(new_scale * big);
            big_real = (big * RealScaleX)/new_scale;
            big_imag = (big * ImagScaleX)/new_scale;
            small_real = (small * RealScaleX)/new_scale;
            small_imag = (small * ImagScaleX)/new_scale;
            return sqrt((N - 1) * (small_real * small_real + small_imag * small_imag) + big_real * big_real + big_imag * big_imag) * new_scale;
          case util_Vec_Pos_Pos_Big:
          case util_Vec_Pos_Neg_Big:
            new_scale = dscale(new_scale * big);
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
              return tmpY0 * ScaleX;
            case util_Vec_Pos_Big:
              return (N - 1) * ScaleX * ScaleY * small + ScaleX * ScaleY * big;
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * ScaleX * ScaleY * small + (ScaleX * ScaleY * big + ScaleX * ScaleY * big);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * ScaleX * ScaleY * small;
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
              return tmpX0 * ScaleY + tmpX1 * ScaleY;
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
              return tmpX0 * tmpY0 + tmpX1 * tmpY1;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
        case util_Vec_Pos_Pos_Big:
          switch(FillY){
            case util_Vec_Constant:
              return (N - 2) * ScaleX * ScaleY * small + (ScaleX * ScaleY * big + ScaleX * ScaleY * big);
            case util_Vec_Pos_Big:
              return ((N - 2) * ScaleX * ScaleY * small * small + ScaleX * ScaleY * big * small) + ScaleX * ScaleY * big * big;
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * ScaleX * ScaleY * small * small + (ScaleX * ScaleY * big * big + ScaleX * ScaleY * big * big);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * ScaleX * ScaleY * small * small;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_zaugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
              return (N - 2) * ScaleX * ScaleY * small * small + (ScaleX * ScaleY * big * big + ScaleX * ScaleY * big * big);
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

double complex wrap_zaugsum_bound(int fold, int N, wrap_zaugsum_func_t func, double complex *X, int incX, double complex *Y, int incY, double res, double ref){
  switch(func){
    case wrap_zaugsum_RZSUM:
    case wrap_zaugsum_ZIZIADD:
    case wrap_zaugsum_ZIZADD:
    case wrap_zaugsum_ZIZDEPOSIT:
      {
        double complex amax;
        double complex bound;
        double *bound_base = (double*)&bound;
        zamax_sub(N, X, incX, &amax);
        bound_base[0] = dibound(fold, N, creal(amax));
        bound_base[1] = dibound(fold, N, cimag(amax));
        return bound;
      }
    case wrap_zaugsum_RDZASUM:
      {
        double complex amax2;
        double amax;
        zamax_sub(N, X, incX, &amax2);
        amax = MAX(creal(amax2), cimag(amax2));
        return dibound(fold, N, amax);
      }
    case wrap_zaugsum_RDZNRM2:
      {
        double complex amax2;
        double amax;
        zamax_sub(N, X, incX, &amax2);
        amax = MAX(creal(amax2), cimag(amax2));
        if (amax == 0.0){
          return 0.0;
        }
        return dibound(fold, N, amax) * (amax / (res + ref));
      }
    case wrap_zaugsum_RZDOTU:
    case wrap_zaugsum_RZDOTC:
      {
        double complex amaxm;
        double complex bound;
        double *bound_base = (double*)&bound;
        zamaxm_sub(N, X, incX, Y, incY, &amaxm);
        bound_base[0] = dibound(fold, N, creal(amaxm));
        bound_base[1] = dibound(fold, N, cimag(amaxm));
        return bound;
      }
  }
  fprintf(stderr, "ReproBLAS error: unknown bound for %s\n", wrap_zaugsum_func_descs[func]);
  exit(125);
}

#endif
