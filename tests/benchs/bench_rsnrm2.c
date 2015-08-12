#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <idxdBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "../../config.h"

#include "bench_vecvec_fill_header.h"

static opt_option fold;

static void bench_rsnrm2_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 0;
  fold._int.max               = SIMAXFOLD;
  fold._int.value             = SIDEFAULTFOLD;
}

int bench_vecvec_fill_show_help(void){
  bench_rsnrm2_options_initialize();

  opt_show_option(fold);

  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  bench_rsnrm2_options_initialize();
  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rsnrm2] (fold = %d)", fold._int.value);
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  (void)FillY;
  (void)RealScaleY;
  (void)ImagScaleY;
  (void)incY;
  int rc = 0;
  int i;
  int j;
  float res = 0.0;
  float_indexed *ires;
  float scale = 0.0;

  bench_rsnrm2_options_initialize();
  opt_eval_option(argc, argv, &fold);

  util_random_seed();

  float *X = util_svec_alloc(N, incX);

  //fill x
  util_svec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);

  if(fold._int.value == SIDEFAULTFOLD){
    time_tic();
    for(i = 0; i < trials; i++){
      res = rsnrm2(N, X, incX);
    }
    time_toc();
  }else if(fold._int.value == 0){
    time_tic();
    for(j = 2; j <= SIMAXFOLD; j++){
      ires = idxd_sialloc(j);
      idxd_sisetzero(j, ires);
      for(i = 0; i < trials; i++){
        scale = idxdBLAS_sisssq(j, N, X, incX, 0.0, ires);
      }
      res = scale * sqrt(idxd_ssiconv(j, ires));
      free(ires);
    }
    time_toc();
  }else{
    time_tic();
    ires = idxd_sialloc(fold._int.value);
    idxd_sisetzero(fold._int.value, ires);
    for(i = 0; i < trials; i++){
      scale = idxdBLAS_sisssq(fold._int.value, N, X, incX, 0.0, ires);
    }
    res = scale * sqrt(idxd_ssiconv(fold._int.value, ires));
    free(ires);
    time_toc();
  }

  metric_load_double("time", time_read());
  metric_load_float("res", res);
  metric_load_double("trials", (double)trials);
  metric_load_double("input", (double)N);
  metric_load_double("output", (double)1);
  if(fold._int.value != 0){
    metric_load_double("s_mul", (double)2 * N);
    metric_load_double("s_add", (double)(3 * fold._int.value - 2) * N);
    metric_load_double("s_orb", (double)fold._int.value * N);
  }
  metric_dump();

  free(X);
  return rc;
}
