#include <stdio.h>
#include "../common/test_opt.h"

#include "../common/test_matvec_fill_header.h"

int bench_matvec_fill_show_help(void);
const char* bench_matvec_fill_name(int argc, char** argv);
int bench_matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials);

static opt_option trials;

static void bench_matvec_fill_options_initialize(void){
  trials._int.header.type       = opt_int;
  trials._int.header.short_name = 'a';
  trials._int.header.long_name  = "trials";
  trials._int.header.help       = "number of trials";
  trials._int.required          = 0;
  trials._int.min               = 1;
  trials._int.max               = INT_MAX;
  trials._int.value             = 10;
}

int matvec_fill_show_help(void){
  bench_matvec_fill_options_initialize();

  opt_show_option(trials);
  return bench_matvec_fill_show_help();
}

const char* matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  bench_matvec_fill_options_initialize();

  opt_eval_option(argc, argv, &trials);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%d trials)", bench_matvec_fill_name(argc, argv), trials._int.value);
  return name_buffer;
}

int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY){
  bench_matvec_fill_options_initialize();

  opt_eval_option(argc, argv, &trials);
  int rc = bench_matvec_fill_test(argc, argv, Order, TransA, M, N, RealAlpha, ImagAlpha, FillA, RealScaleA, ImagScaleA, lda, FillX, RealScaleX, ImagScaleX, incX, RealBeta, ImagBeta, FillY, RealScaleY, ImagScaleY, incY, trials._int.value);
  return rc;
}
