#ifndef CAUGSUM_WRAPPER_H
#define CAUGSUM_WRAPPER_H

#include <reproBLAS.h>
#include <indexedBLAS.h>
#include <indexed.h>

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

float wrap_caugsum_result(int N, wrap_caugsum_func_t func, util_vec_fill_t FillX, double ScaleX, double CondX, util_vec_fill_t FillY, double ScaleY, double CondY){
  float small = 1.0 / (1024.0 * 4.0); // 2^-12
  float big   = 1024.0 * 8.0;  // 2^13
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
        case util_Vec_Constant_Real_Imag:
          return N * ScaleX * (1.0 + I);
        case util_Vec_Pos_Inf_Real_Imag:
        case util_Vec_Pos_Pos_Inf_Real_Imag:
          return ScaleX * (1.0 + I)/0.0;
        case util_Vec_Pos_Neg_Inf_Real_Imag:
        case util_Vec_NaN_Real_Imag:
        case util_Vec_Pos_Inf_NaN_Real_Imag:
        case util_Vec_Pos_Pos_Inf_NaN_Real_Imag:
        case util_Vec_Pos_Neg_Inf_NaN_Real_Imag:
          return 0.0/0.0 * (1.0 + I);
        case util_Vec_Pos_Big_Real_Imag:
          return ((N - 1) * ScaleX * small + ScaleX * big) * (1.0 + I);
        case util_Vec_Pos_Pos_Big_Real_Imag:
          return (((N - 2) * ScaleX * small + ScaleX * big) + ScaleX * big) * (1.0 + I);
        case util_Vec_Pos_Neg_Big_Real_Imag:
          return ((N - 2) * ScaleX * small) * (1.0 + I);
        case util_Vec_Sine_Real_Imag:
          return (ScaleX - ScaleX) * (1.0 + I);
        case util_Vec_Constant_Imag:
          return N * ScaleX * I;
        case util_Vec_Pos_Inf_Imag:
        case util_Vec_Pos_Pos_Inf_Imag:
          return ScaleX * I/0.0;
        case util_Vec_Pos_Neg_Inf_Imag:
        case util_Vec_NaN_Imag:
        case util_Vec_Pos_Inf_NaN_Imag:
        case util_Vec_Pos_Pos_Inf_NaN_Imag:
        case util_Vec_Pos_Neg_Inf_NaN_Imag:
          return 0.0/0.0 * I;
        case util_Vec_Pos_Big_Imag:
          return ((N - 1) * ScaleX * small + ScaleX * big) * I;
        case util_Vec_Pos_Pos_Big_Imag:
          return (((N - 2) * ScaleX * small + ScaleX * big) + ScaleX * big) * I;
        case util_Vec_Pos_Neg_Big_Imag:
          return ((N - 2) * ScaleX * small) * I;
        case util_Vec_Sine_Imag:
          return (ScaleX - ScaleX) * I;
        default:
          printf("foobar %d %d\n", FillX, util_Vec_Sine);
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX);
          exit(125);
      }
    case wrap_caugsum_RSCASUM:
      switch(FillX){
        case util_Vec_Constant:
        case util_Vec_Constant_Imag:
          return N * fabs(ScaleX);
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
        case util_Vec_Pos_Neg_Inf:
        case util_Vec_Pos_Inf_Real_Imag:
        case util_Vec_Pos_Pos_Inf_Real_Imag:
        case util_Vec_Pos_Neg_Inf_Real_Imag:
        case util_Vec_Pos_Inf_Imag:
        case util_Vec_Pos_Pos_Inf_Imag:
        case util_Vec_Pos_Neg_Inf_Imag:
          return fabs(ScaleX)/0.0;
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
        case util_Vec_NaN_Real_Imag:
        case util_Vec_Pos_Inf_NaN_Real_Imag:
        case util_Vec_Pos_Pos_Inf_NaN_Real_Imag:
        case util_Vec_Pos_Neg_Inf_NaN_Real_Imag:
        case util_Vec_NaN_Imag:
        case util_Vec_Pos_Inf_NaN_Imag:
        case util_Vec_Pos_Pos_Inf_NaN_Imag:
        case util_Vec_Pos_Neg_Inf_NaN_Imag:
          return 0.0/0.0;
        case util_Vec_Pos_Big:
        case util_Vec_Pos_Big_Imag:
          return (N - 1) * fabs(ScaleX) * small + fabs(ScaleX) * big;
        case util_Vec_Pos_Pos_Big:
        case util_Vec_Pos_Neg_Big:
        case util_Vec_Pos_Pos_Big_Imag:
        case util_Vec_Pos_Neg_Big_Imag:
          return ((N - 2) * fabs(ScaleX) * small + fabs(ScaleX) * big) + fabs(ScaleX) * big;
        case util_Vec_Constant_Real_Imag:
          return 2.0 * N * fabs(ScaleX);
        case util_Vec_Pos_Big_Real_Imag:
          return 2.0 * ((N - 1) * fabs(ScaleX) * small + fabs(ScaleX) * big);
        case util_Vec_Pos_Pos_Big_Real_Imag:
        case util_Vec_Pos_Neg_Big_Real_Imag:
          return 2.0 * (((N - 2) * fabs(ScaleX) * small + fabs(ScaleX) * big) + fabs(ScaleX) * big);
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX);
          exit(125);
      }
    case wrap_caugsum_RSCNRM2:
      switch(FillX){
        case util_Vec_Constant:
        case util_Vec_Constant_Imag:
          return sqrt(N) * ScaleX;
        case util_Vec_Pos_Inf:
        case util_Vec_Pos_Pos_Inf:
        case util_Vec_Pos_Neg_Inf:
        case util_Vec_Pos_Inf_Real_Imag:
        case util_Vec_Pos_Pos_Inf_Real_Imag:
        case util_Vec_Pos_Neg_Inf_Real_Imag:
        case util_Vec_Pos_Inf_Imag:
        case util_Vec_Pos_Pos_Inf_Imag:
        case util_Vec_Pos_Neg_Inf_Imag:
          return fabs(ScaleX)/0.0;
        case util_Vec_NaN:
        case util_Vec_Pos_Inf_NaN:
        case util_Vec_Pos_Pos_Inf_NaN:
        case util_Vec_Pos_Neg_Inf_NaN:
        case util_Vec_NaN_Real_Imag:
        case util_Vec_Pos_Inf_NaN_Real_Imag:
        case util_Vec_Pos_Pos_Inf_NaN_Real_Imag:
        case util_Vec_Pos_Neg_Inf_NaN_Real_Imag:
        case util_Vec_NaN_Imag:
        case util_Vec_Pos_Inf_NaN_Imag:
        case util_Vec_Pos_Pos_Inf_NaN_Imag:
        case util_Vec_Pos_Neg_Inf_NaN_Imag:
          return 0.0/0.0;
        case util_Vec_Pos_Big:
        case util_Vec_Pos_Big_Imag:
          return sqrt((N - 1) * small * small + big * big) * fabs(ScaleX);
        case util_Vec_Pos_Pos_Big:
        case util_Vec_Pos_Neg_Big:
        case util_Vec_Pos_Pos_Big_Imag:
        case util_Vec_Pos_Neg_Big_Imag:
          return sqrt(((N - 2) * small * small + big * big) + big * big) * fabs(ScaleX);
        case util_Vec_Constant_Real_Imag:
          return sqrt(2 * N) * ScaleX;
        case util_Vec_Pos_Big_Real_Imag:
          return sqrt(2.0 * ((N - 1) * small * small + big * big)) * fabs(ScaleX);
        case util_Vec_Pos_Pos_Big_Real_Imag:
        case util_Vec_Pos_Neg_Big_Real_Imag:
          return sqrt(2.0 * (((N - 2) * small * small + big * big) + big * big)) * fabs(ScaleX);
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX);
          exit(125);
      }
    case wrap_caugsum_RCDOTU:
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
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
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
              exit(125);
          }
        case util_Vec_Sine:
          switch(FillY){
            case util_Vec_Constant:
              return ScaleX * ScaleY - ScaleX * ScaleY;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
              exit(125);
          }
        default:
          fprintf(stderr, "ReproBLAS error: unknown result for %s(%s * %g, %s * %g)\n", wrap_caugsum_func_descs[func], util_vec_fill_descs[FillX], ScaleX, util_vec_fill_descs[FillY], ScaleY);
          exit(125);
      }
  }
}

float wrap_caugsum_bound(int fold, int N, wrap_caugsum_func_t func, float *X, int incX, float *Y, int incY, float res, float ref){
  switch(func){
    case wrap_caugsum_RSSUM:
    case wrap_caugsum_SISIADD:
    case wrap_caugsum_SISADD:
    case wrap_caugsum_SISDEPOSIT:
      {
        float complex amax = camax(N, X, incX);
        return sibound(fold, N, crealf(amax)) + sibound(fold, N, cimagf(amax)) * I;
      }
    case wrap_caugsum_RSCASUM:
      {
        float complex amax = camax(N, X, incX);
        return sibound(fold, N, MAX(crealf(amax), cimagf(amax)));
      }
    case wrap_caugsum_RSCNRM2:
      {
        float amax = crealf(camaxm(N, X, incX, X, incX));
        return sibound(fold, N, amax) * (amax / (res + ref));
      }
    case wrap_caugsum_RCDOTU:
    case wrap_caugsum_RCDOTC:
      {
        float complex amaxm = crealf(camaxm(N, X, incX, Y, incY));
        return sibound(fold, N, crealf(amaxm)) + sibound(fold, N, cimagf(amaxm)) * I;
      }
  }
}


#endif
