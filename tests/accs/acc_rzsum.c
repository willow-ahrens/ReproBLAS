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

static void acc_rzsum_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = MAX_FOLD;
  fold._int.value             = DEFAULT_FOLD;
}

int acc_vecvec_fill_show_help(void){
  acc_rzsum_options_initialize();

  opt_show_option(fold);

  return 0;
}

const char* acc_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];

  acc_rzsum_options_initialize();
  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Accuracy [rzsum] (fold = %d)", fold._int.value);
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
  double complex res = 0.0;
  double_complex_indexed *ires;
  double ref;
  double s[2];

  acc_rzsum_options_initialize();
  opt_eval_option(argc, argv, &fold);

  util_random_seed();

  double complex *X = util_zvec_alloc(N, incX);
  double *ratios = util_dvec_alloc(2 * N, 1);

  for(i = 0; i < trials; i++){
    util_zvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
    ires = zialloc(fold._int.value);
    zisetzero(fold._int.value, ires);
    zizsum(fold._int.value, N, X, incX, ires);
    zziconv_sub(fold._int.value, ires, &res);
    free(ires);

    util_dvec_sort(N, (double*)X, incX * 2, NULL, 0, util_Decreasing_Magnitude);
    s[0] = 0.0;
    s[1] = 0.0;
    for(j = 0; j < N; j++){
      util_ddpd(s, creal(X[j * incX]));
    }
    ref = s[0] + s[1];
    ratios[2 * i] = fabs(creal(res) - ref)/MAX(fabs(ref), DBL_MIN);

    util_dvec_sort(N, ((double*)X) + 1, incX * 2, NULL, 0, util_Decreasing_Magnitude);
    s[0] = 0.0;
    s[1] = 0.0;
    for(j = 0; j < N; j++){
      util_ddpd(s, cimag(X[j * incX]));
    }
    ref = s[0] + s[1];
    ratios[2 * i + 1] = fabs(cimag(res) - ref)/MAX(fabs(ref), DBL_MIN);
  }

  util_dvec_sort(2 * N, ratios, 1, NULL, 0, util_Increasing);
  metric_load_double("min_ratio", ratios[0]);
  metric_load_double("med_ratio", ratios[N]);
  metric_load_double("max_ratio", ratios[2 * N - 1]);
  metric_load_double("e", DBL_EPSILON);
  metric_load_double("trials", (double)trials);
  metric_dump();

  free(X);
  free(ratios);
  return rc;
}
