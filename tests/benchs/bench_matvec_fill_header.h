#include <stdio.h>
#include "../common/test_opt.h"
#include "../common/test_perf.h"

#include "../common/test_matvec_fill_header.h"

int bench_matvec_fill_desc(void);
int bench_matvec_fill_show_help(void);
const char* bench_matvec_fill_name(int argc, char** argv);
int bench_matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, int FillA, double ScaleA, double CondA, int lda, int FillX, double ScaleX, double CondX, int incX, int trials);

static opt_option trials = {._int.header.type       = opt_int,
                            ._int.header.short_name = 'a',
                            ._int.header.long_name  = "trials",
                            ._int.header.help       = "number of trials",
                            ._int.required          = 0,
                            ._int.min               = 1,
                            ._int.max               = INT_MAX,
                            ._int.value             = 10000};

static opt_option desc   = {._flag.header.type       = opt_flag,
                            ._flag.header.short_name = 'd',
                            ._flag.header.long_name  = "desc",
                            ._flag.header.help       = "show benchmark description"};

int matvec_fill_show_help(void){
  opt_show_option(trials);
  return bench_matvec_fill_show_help();
}

const char* matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &trials);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%d trials)", bench_matvec_fill_name(argc, argv), trials._int.value);
  return name_buffer;
}

int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, int FillA, double ScaleA, double CondA, int lda, int FillX, double ScaleX, double CondX, int incX){
  opt_eval_option(argc, argv, &desc);
  if(desc._flag.exists){
    return bench_matvec_fill_desc();
  }

  opt_eval_option(argc, argv, &trials);
  int rc = bench_matvec_fill_test(argc, argv, Order, TransA, M, N, FillA, ScaleA, CondA, lda, FillX, ScaleX, CondX, incX, trials._int.value);
  return rc;
}
