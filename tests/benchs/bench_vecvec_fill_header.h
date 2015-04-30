#include <stdio.h>
#include "../common/test_opt.h"
#include "../common/test_perf.h"

#include "../common/test_vecvec_fill_header.h"

int bench_vecvec_fill_show_help(void);
const char* bench_vecvec_fill_name(int argc, char** argv);
int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY, int trials);

static opt_option trials = {._int.header.type       = opt_int,
                            ._int.header.short_name = 'a',
                            ._int.header.long_name  = "trials",
                            ._int.header.help       = "number of trials",
                            ._int.required          = 0,
                            ._int.min               = 1,
                            ._int.max               = INT_MAX,
                            ._int.value             = 10000};

int vecvec_fill_show_help(void){
  opt_show_option(trials);
  return bench_vecvec_fill_show_help();
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &trials);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%d trials)", bench_vecvec_fill_name(argc, argv), trials._int.value);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY){
  opt_eval_option(argc, argv, &trials);
  int rc = bench_vecvec_fill_test(argc, argv, N, FillX, ScaleX, CondX, incX, FillY, ScaleY, CondY, incY, trials._int.value);
  return rc;
}
