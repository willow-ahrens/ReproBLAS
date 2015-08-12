#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "../../config.h"

#include "acc_vecvec_fill_header.h"

static opt_option fold;

static void acc_rcsum_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = SIMAXFOLD;
  fold._int.value             = SIDEFAULTFOLD;
}

int acc_vecvec_fill_show_help(void){
  acc_rcsum_options_initialize();

  opt_show_option(fold);

  return 0;
}

const char* acc_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];

  acc_rcsum_options_initialize();
  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Accuracy [rcsum] (fold = %d)", fold._int.value);
  return name_buffer;
}

int acc_vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  (void)argc;
  (void)argv;
  (void)FillY;
  (void)RealScaleY;
  (void)ImagScaleY;
  (void)incY;
  int rc = 0;
  int i;
  int j;
  float complex res = 0.0;
  float_complex_indexed *ires;
  double s;
  float ref;

  acc_rcsum_options_initialize();
  opt_eval_option(argc, argv, &fold);

  util_random_seed();

  float complex *X = util_cvec_alloc(N, incX);
  float *ratios = util_svec_alloc(2 * N, 1);

  for(i = 0; i < trials; i++){
    util_cvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
    ires = idxd_cialloc(fold._int.value);
    idxd_cisetzero(fold._int.value, ires);
    idxdBLAS_cicsum(fold._int.value, N, X, incX, ires);
    idxd_cciconv_sub(fold._int.value, ires, &res);
    free(ires);

    util_svec_sort(N, (float*)X, incX * 2, NULL, 0, util_Decreasing_Magnitude);
    s = 0.0;
    for(j = 0; j < N; j++){
      s += crealf(X[j * incX]);
    }
    ref = s;
    ratios[2 * i] = fabsf(crealf(res) - ref)/MAX(fabsf(ref), FLT_MIN);

    util_svec_sort(N, ((float*)X) + 1, incX * 2, NULL, 0, util_Decreasing_Magnitude);
    s = 0.0;
    for(j = 0; j < N; j++){
      s += cimagf(X[j * incX]);
    }
    ref = s;
    ratios[2 * i + 1] = fabsf(cimagf(res) - ref)/MAX(fabsf(ref), FLT_MIN);
  }

  util_svec_sort(2 * N, ratios, 1, NULL, 0, util_Increasing);
  metric_load_float("min_ratio", ratios[0]);
  metric_load_float("med_ratio", ratios[N]);
  metric_load_float("max_ratio", ratios[2 * N - 1]);
  metric_load_float("e", FLT_EPSILON);
  metric_load_double("trials", (double)trials);
  metric_dump();

  free(X);
  free(ratios);
  return rc;
}
