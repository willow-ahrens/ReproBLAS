#ifndef CAUGSUM_WRAPPER_H
#define CAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <idxdBLAS.h>
#include <idxd.h>
#include "../../config.h"

#include "../common/test_util.h"

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
  if(fold == SIDEFAULTFOLD){
    float complex res;
    reproBLAS_csum_sub(N, x, incx, &res);
    return res;
  }else{
    float complex res;
    reproBLAS_rcsum_sub(fold, N, x, incx, &res);
    return res;
  }
}

void wrap_cicsum(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  idxdBLAS_cicsum(fold, N, x, incx, c);
}

float complex wrap_rscasum(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == SIDEFAULTFOLD){
    return reproBLAS_scasum(N, x, incx);
  }else{
    return reproBLAS_rscasum(fold, N, x, incx);
  }
}

void wrap_sicasum(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  idxdBLAS_smcasum(fold, N, x, incx, c, 2, c + 2 * fold, 2);
}

float complex wrap_rscnrm2(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  if(fold == SIDEFAULTFOLD){
    return reproBLAS_scnrm2(N, x, incx);
  }else{
    return reproBLAS_rscnrm2(fold, N, x, incx);
  }
}

void wrap_sicssq(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  idxdBLAS_smcssq(fold, N, x, incx, 0.0, c, 2, c + 2 * fold, 2);
}

float complex wrap_rcdotu(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  if(fold == SIDEFAULTFOLD){
    float complex res;
    reproBLAS_cdotu_sub(N, x, incx, y, incy, &res);
    return res;
  }else{
    float complex res;
    reproBLAS_rcdotu_sub(fold, N, x, incx, y, incy, &res);
    return res;
  }
}

void wrap_cicdotu(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  idxdBLAS_cicdotu(fold, N, x, incx, y, incy, c);
}

float complex wrap_rcdotc(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  if(fold == SIDEFAULTFOLD){
    float complex res;
    reproBLAS_cdotc_sub(N, x, incx, y, incy, &res);
    return res;
  }else{
    float complex res;
    reproBLAS_rcdotc_sub(fold, N, x, incx, y, incy, &res);
    return res;
  }
}

void wrap_cicdotc(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  idxdBLAS_cicdotc(fold, N, x, incx, y, incy, c);
}

float complex wrap_rciciadd(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  float_complex_indexed *ires = idxd_cialloc(fold);
  float_complex_indexed *itmp = idxd_cialloc(fold);
  idxd_cisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    idxd_cicconv(fold, x + i * incx, itmp);
    idxd_ciciadd(fold, itmp, ires);
  }
  float complex res;
  idxd_cciconv_sub(fold, ires, &res);
  free(ires);
  free(itmp);
  return res;
}

void wrap_ciciadd(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  float_complex_indexed *itmp = idxd_cialloc(fold);
  int i;
  for(i = 0; i < N; i++){
    idxd_cicconv(fold, x + i * incx, itmp);
    idxd_ciciadd(fold, itmp, c);
  }
  free(itmp);
}

float complex wrap_rcicadd(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  float_complex_indexed *ires = idxd_cialloc(fold);
  idxd_cisetzero(fold, ires);
  int i;
  for(i = 0; i < N; i++){
    idxd_cicadd(fold, x + i * incx, ires);
  }
  float complex res;
  idxd_cciconv_sub(fold, ires, &res);
  free(ires);
  return res;
}

void wrap_cicadd(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  int i;
  for(i = 0; i < N; i++){
    idxd_cicadd(fold, x + i * incx, c);
  }
}

float complex wrap_rcicdeposit(int fold, int N, float complex *x, int incx, float complex *y, int incy) {
  (void)y;
  (void)incy;
  float_complex_indexed *ires = idxd_cialloc(fold);
  idxd_cisetzero(fold, ires);
  float complex amax;
  idxdBLAS_camax_sub(N, x, incx, &amax);
  idxd_cicupdate(fold, &amax, ires);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= idxd_SIENDURANCE){
      idxd_cirenorm(fold, ires);
      j = 0;
    }
    idxd_cicdeposit(fold, x + i * incx, ires);
    j++;
  }
  idxd_cirenorm(fold, ires);
  float complex res;
  idxd_cciconv_sub(fold, ires, &res);
  free(ires);
  return res;
}

void wrap_cicdeposit(int fold, int N, float complex *x, int incx, float complex *y, int incy, float_complex_indexed *c) {
  (void)y;
  (void)incy;
  float complex amax;
  idxdBLAS_camax_sub(N, x, incx, &amax);
  idxd_cicupdate(fold, &amax, c);
  int i;
  int j = 0;
  for(i = 0; i < N; i++){
    if(j >= idxd_SIENDURANCE){
      idxd_cirenorm(fold, c);
      j = 0;
    }
    idxd_cicdeposit(fold, x + i * incx, c);
    j++;
  }
  idxd_cirenorm(fold, c);
}

wrap_caugsum wrap_caugsum_func(wrap_caugsum_func_t func) {
  switch(func){
    case wrap_caugsum_RCSUM:
      return wrap_rcsum;
    case wrap_caugsum_RSCASUM:
      return wrap_rscasum;
    case wrap_caugsum_RSCNRM2:
      return wrap_rscnrm2;
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
      return wrap_sicssq;
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

float complex wrap_caugsum_result(int N, wrap_caugsum_func_t func, util_vec_fill_t FillX, double RealScaleX, double ImagScaleX, util_vec_fill_t FillY, double RealScaleY, double ImagScaleY){
  float small              = 1.0 / (1024.0 * 4.0); // 2^-12
  float big                = 1024.0 * 8.0;         // 2^13
  float complex ScaleX     = RealScaleX + ImagScaleX * I;
  float complex ScaleY     = RealScaleY + ImagScaleY * I;
  float complex tmpX0;
  float *tmpX0_base = (float*)&tmpX0;
  float complex tmpX1;
  float *tmpX1_base = (float*)&tmpX1;
  float complex tmpY0;
  float *tmpY0_base = (float*)&tmpY0;
  float complex tmpY1;
  float *tmpY1_base = (float*)&tmpY1;
  switch(func){
    case wrap_caugsum_RCSUM:
    case wrap_caugsum_CICIADD:
    case wrap_caugsum_CICADD:
    case wrap_caugsum_CICDEPOSIT:
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
          return (N - 1) * (ScaleX * small) + ScaleX * big;
        case util_Vec_Pos_Pos_Big:
          return (N - 2) * (ScaleX * small) + (ScaleX * big + ScaleX * big);
        case util_Vec_Pos_Neg_Big:
          return (N - 2) * (ScaleX * small);
        case util_Vec_Sine:
          return 0.0;
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX);
          exit(125);
      }

    case wrap_caugsum_RSCASUM:
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
          return (N - 1) * ((fabs(RealScaleX) + fabs(ImagScaleX)) * small) + (fabs(RealScaleX) + fabs(ImagScaleX)) * big;
        case util_Vec_Pos_Pos_Big:
        case util_Vec_Pos_Neg_Big:
          return (N - 2) * ((fabs(RealScaleX) + fabs(ImagScaleX)) * small) + ((fabs(RealScaleX) + fabs(ImagScaleX)) * big + (fabs(RealScaleX) + fabs(ImagScaleX)) * big);
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX);
          exit(125);
      }

    case wrap_caugsum_RSCNRM2:
      {
        double new_scale = MAX(fabs(RealScaleX), fabs(ImagScaleX));
        double big_real;
        double big_imag;
        double small_real;
        double small_imag;
        switch(FillX){
          case util_Vec_Constant:
            new_scale = idxd_sscale(new_scale);
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
                fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
            new_scale = idxd_sscale(new_scale * big);
            big_real = (big * RealScaleX)/new_scale;
            big_imag = (big * ImagScaleX)/new_scale;
            small_real = (small * RealScaleX)/new_scale;
            small_imag = (small * ImagScaleX)/new_scale;
            return sqrt((N - 1) * (small_real * small_real + small_imag * small_imag) + big_real * big_real + big_imag * big_imag) * new_scale;
          case util_Vec_Pos_Pos_Big:
          case util_Vec_Pos_Neg_Big:
            new_scale = idxd_sscale(new_scale * big);
            big_real = (big * RealScaleX)/new_scale;
            big_imag = (big * ImagScaleX)/new_scale;
            small_real = (small * RealScaleX)/new_scale;
            small_imag = (small * ImagScaleX)/new_scale;
            return sqrt((N - 2) * (small_real * small_real + small_imag * small_imag) + (big_real * big_real + big_imag * big_imag + big_real * big_real + big_imag * big_imag)) * new_scale;
          default:
            fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX);
            exit(125);
        }
      }

    case wrap_caugsum_RCDOTC:
      ScaleX = RealScaleX - ImagScaleX * I;
      ImagScaleX = -1 * ImagScaleX;

    case wrap_caugsum_RCDOTU:
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
                  fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
                  exit(125);
              }
              if(RealScaleY == 0.0){
                tmpY0_base[0] = 0.0;
              }
              if(ImagScaleY == 0.0){
                tmpY0_base[1] = 0.0;
              }
              return cmul(tmpY0, ScaleX);
            case util_Vec_Pos_Big:
              return (N - 1) * (ScaleX * ScaleY * small) + ScaleX * ScaleY * big;
            case util_Vec_Pos_Pos_Big:
              return (N - 2) * (ScaleX * ScaleY * small) + (ScaleX * ScaleY * big + ScaleX * ScaleY * big);
            case util_Vec_Pos_Neg_Big:
              return (N - 2) * (ScaleX * ScaleY * small);
            case util_Vec_Sine:
              return 0.0;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
              return cmul(tmpX0, ScaleY) + cmul(tmpX1, ScaleY);
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
                  fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
              return cmul(tmpX0, tmpY0) + cmul(tmpX1, tmpY1);
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
        case util_Vec_Sine:
          switch(FillY){
            case util_Vec_Constant:
              return 0.0;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
              exit(125);
          }
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
          exit(125);
      }
    default:
      fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * (%g + %gi), %s * (%g + %gi))\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], RealScaleX, ImagScaleX, util_vec_fill_descs[FillY], RealScaleY, ImagScaleY);
      exit(125);
  }
}

float complex wrap_caugsum_bound(int fold, int N, wrap_caugsum_func_t func, float complex *X, int incX, float complex *Y, int incY, float complex res, float complex ref){
  switch(func){
    case wrap_caugsum_RCSUM:
    case wrap_caugsum_CICIADD:
    case wrap_caugsum_CICADD:
    case wrap_caugsum_CICDEPOSIT:
      {
        float complex amax;
        float complex bound;
        float *bound_base = (float*)&bound;
        idxdBLAS_camax_sub(N, X, incX, &amax);
        bound_base[0] = idxd_sibound(fold, N, crealf(amax), crealf(res));
        bound_base[1] = idxd_sibound(fold, N, cimagf(amax), cimagf(res));
        return bound;
      }
    case wrap_caugsum_RSCASUM:
      {
        float complex amax2;
        float amax;
        idxdBLAS_camax_sub(N, X, incX, &amax2);
        amax = MAX(crealf(amax2), cimagf(amax2));
        return idxd_sibound(fold, N, amax, crealf(res));
      }
    case wrap_caugsum_RSCNRM2:
      {
        float complex amax2;
        float amax;
        float scale;
        idxdBLAS_camax_sub(N, X, incX, &amax2);
        amax = MAX(crealf(amax2), cimagf(amax2));
        scale = idxd_sscale(amax);
        if (amax == 0.0){
          return 0.0;
        }
        return idxd_sibound(fold, N, (amax/scale) * (amax/scale), (crealf(res)/scale) * (crealf(res)/scale)) * (scale / (crealf(res) + crealf(ref))) * scale;
      }
    case wrap_caugsum_RCDOTU:
    case wrap_caugsum_RCDOTC:
      {
        float complex amaxm;
        float complex bound;
        float *bound_base = (float*)&bound;
        idxdBLAS_camaxm_sub(N, X, incX, Y, incY, &amaxm);
        bound_base[0] = idxd_sibound(fold, 2*N, crealf(amaxm), crealf(res));
        bound_base[1] = idxd_sibound(fold, 2*N, cimagf(amaxm), cimagf(res));
        return bound;
      }
  }
  fprintf(stderr, "ReproBLAS error: unknown bound for %s\n", wrap_caugsum_func_descs[func]);
  exit(125);
}

#endif
