#include <stdio.h>
#include "../common/test_opt.h"

#include "../common/test_matmat_fill_header.h"

int bench_matmat_fill_show_help(void);
const char* bench_matmat_fill_name(int argc, char** argv);
int bench_matmat_fill_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillB, double RealScaleB, double ImagScaleB, int ldb, double RealBeta, double ImagBeta, int FillC, double RealScaleC, double ImagScaleC, int ldc, int trials);

static opt_option trials;

static void bench_matmat_fill_options_initialize(void){
  trials._int.header.type       = opt_int;
  trials._int.header.short_name = 'a';
  trials._int.header.long_name  = "trials";
  trials._int.header.help       = "number of trials";
  trials._int.required          = 0;
  trials._int.min               = 1;
  trials._int.max               = INT_MAX;
  trials._int.value             = 100;
}

int matmat_fill_show_help(void){
  bench_matmat_fill_options_initialize();

  opt_show_option(trials);
  return bench_matmat_fill_show_help();
}

const char* matmat_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  bench_matmat_fill_options_initialize();

  opt_eval_option(argc, argv, &trials);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%d trials)", bench_matmat_fill_name(argc, argv), trials._int.value);
  return name_buffer;
}

int matmat_fill_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillB, double RealScaleB, double ImagScaleB, int ldb, double RealBeta, double ImagBeta, int FillC, double RealScaleC, double ImagScaleC, int ldc){
  bench_matmat_fill_options_initialize();

  opt_eval_option(argc, argv, &trials);
  int rc = bench_matmat_fill_test(argc, argv, Order, TransA, TransB, M, N, K, RealAlpha, ImagAlpha, FillA, RealScaleA, ImagScaleA, lda, FillB, RealScaleB, ImagScaleB, ldb, RealBeta, ImagBeta, FillC, RealScaleC, ImagScaleC, ldc, trials._int.value);
  return rc;
}
