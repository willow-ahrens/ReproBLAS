#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_perf.h"

#include "bench_matvec_fill_header.h"

static opt_option incY   = {._int.header.type       = opt_int,
                            ._int.header.short_name = 'y',
                            ._int.header.long_name  = "incY",
                            ._int.header.help       = "Y vector increment",
                            ._int.required          = 0,
                            ._int.min               = 1,
                            ._int.max               = INT_MAX,
                            ._int.value             = 1};

static opt_option FillY  = {._named.header.type       = opt_named,
                            ._named.header.short_name = 'j',
                            ._named.header.long_name  = "FillY",
                            ._named.header.help       = "Y fill type",
                            ._named.required          = 0,
                            ._named.n_names           = (int)util_vec_fill_n_names,
                            ._named.names             = (char**)util_vec_fill_names,
                            ._named.descs             = (char**)util_vec_fill_descs,
                            ._named.value             = 0};

static opt_option ScaleY = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'v',
                            ._double.header.long_name  = "ScaleY",
                            ._double.header.help       = "Y scale",
                            ._double.required          = 0,
                            ._double.min               = 0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1.0};

static opt_option CondY  = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'e',
                            ._double.header.long_name  = "CondY",
                            ._double.header.help       = "Y condition number",
                            ._double.required          = 0,
                            ._double.min               = 1.0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1e3};

static opt_option alpha  = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'l',
                            ._double.header.long_name  = "alpha",
                            ._double.header.help       = "alpha",
                            ._double.required          = 0,
                            ._double.min               = 1.0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1e3};

static opt_option beta   = {._double.header.type       = opt_double,
                            ._double.header.short_name = 'm',
                            ._double.header.long_name  = "beta",
                            ._double.header.help       = "beta",
                            ._double.required          = 0,
                            ._double.min               = 1.0,
                            ._double.max               = DBL_MAX,
                            ._double.value             = 1e3};

int bench_matvec_fill_desc(void){
  /*
  char *op_names[] = {"d_add", "d_orb"};
  int op_counts[] = {7, 3};
  perf_output_desc(2, op_names, op_counts);
  */
  return 0;
}

int bench_matvec_fill_show_help(void){
  opt_show_option(incY);
  opt_show_option(FillY);
  opt_show_option(ScaleY);
  opt_show_option(CondY);
  opt_show_option(alpha);
  opt_show_option(beta);
  return 0;
}

const char* bench_matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rdgemv]");
  return name_buffer;
}

int bench_matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, int FillA, double ScaleA, double CondA, int lda, int FillX, double ScaleX, double CondX, int incX, int trials){
  int rc = 0;

  opt_eval_option(argc, argv, &incY);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &ScaleY);
  opt_eval_option(argc, argv, &CondY);
  opt_eval_option(argc, argv, &alpha);
  opt_eval_option(argc, argv, &beta);

  util_random_seed();

  rblas_order_t o;
  rblas_transpose_t t;
  int NX;
  int NY;
  char NTransA;
  switch(Order){
    case 'r':
    case 'R':
      o = rblas_Row_Major;
      break;
    default:
      o = rblas_Col_Major;
      break;
  }
  switch(TransA){
    case 'n':
    case 'N':
      t = rblas_No_Trans;
      NX = N;
      NY = M;
      NTransA = 'T';
      break;
    default:
      NX = M;
      NY = N;
      NTransA = 'N';
      t = rblas_Trans;
      break;
  }

  double *A  = util_dmat_alloc(Order, M, N, lda);
  double *X  = util_dvec_alloc(NX, incX);
  double *Y  = util_dvec_alloc(NY, incY._int.value);

  util_dmat_fill(Order, 'n', M, N, A, lda, FillA, ScaleA, CondA);
  util_dvec_fill(NX, X, incX, FillX, ScaleX, CondX);
  util_dvec_fill(NY, Y, incY._int.value, FillY._named.value, ScaleY._double.value, CondY._double.value);
  double *res  = (double*)malloc(NY * incY._int.value * sizeof(double));
  double time = 0.0;

  for(int i = 0; i < trials; i++){
    memcpy(res, Y, NY * incY._int.value * sizeof(double));
    time_tic();
    rdgemv(o, t, M, N, A, lda, X, incX, res, incY._int.value);
    time_toc();
    time += time_read();
  }

  perf_output_perf(time, N * M, trials);

  free(X);
  return rc;
}
