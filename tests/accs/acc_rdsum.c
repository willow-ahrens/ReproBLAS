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

static void acc_rdsum_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 1;
  fold._int.max               = MAX_FOLD;
  fold._int.value             = DEFAULT_FOLD;
}

int acc_vecvec_fill_show_help(void){
  acc_rdsum_options_initialize();

  opt_show_option(fold);

  return 0;
}

const char* acc_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];

  acc_rdsum_options_initialize();
  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Accuracy [rdsum] (fold = %d)", fold._int.value);
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
  double res = 0.0;
  double_indexed *ires;
  double ref;
  double ratio = 0.0;

  acc_rdsum_options_initialize();
  opt_eval_option(argc, argv, &fold);

  util_random_seed();

  double *X = util_dvec_alloc(N, incX);

  for(i = 0; i < trials; i++){
    util_dvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
    ires = dialloc(fold._int.value);
    disetzero(fold._int.value, ires);
    didsum(fold._int.value, N, X, incX, ires);
    res = ddiconv(fold._int.value, ires);
    free(ires);

    util_dvec_sort(N, (double*)X, incX, NULL, 1, util_Decreasing_Magnitude);
    double s[2] = {0.0, 0.0};
    for(j = 0; j < N; j++){
      util_ddpd(s, X[j * incX]);
    }
    ref = s[0] + s[1];
    ratio += fabs(res - ref)/MAX(fabs(ref), DBL_MIN);
  }

  metric_load_double("ratio", ratio);
  metric_load_double("e", DBL_EPSILON);
  metric_load_double("trials", (double)trials);
  metric_dump();

  free(X);
  return rc;
}
