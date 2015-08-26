#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "../../config.h"

#include "acc_vecvec_fill_header.h"

int acc_vecvec_fill_show_help(void){
  return 0;
}

const char* acc_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Accuracy [csum]");
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
  double s;
  float ref;

  util_random_seed();

  float complex *X = util_cvec_alloc(N, incX);
  float *ratios = util_svec_alloc(2 * N, 1);

  for(i = 0; i < trials; i++){
    util_cvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
    for(j = 0; j < N; j++){
      res += X[j * incX];
    }

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
