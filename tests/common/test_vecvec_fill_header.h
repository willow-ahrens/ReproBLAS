#include <stdio.h>
#include "test_opt.h"

#include "test_vecvec_header.h"

int vecvec_fill_show_help(void);
const char* vecvec_fill_name(int argc, char** argv);
int vecvec_fill_test(int argc, char** argv, int N, int incX, int incY, int vec_fill_type, double scale, double cond);

static opt_option vec_fill_type = {._named.header.type       = opt_named,
                                   ._named.header.short_name = 'f',
                                   ._named.header.long_name  = "f_type",
                                   ._named.header.help       = "vector fill type",
                                   ._named.required          = 0,
                                   ._named.n_names           = (int)vec_fill_n_names,
                                   ._named.names             = (char**)vec_fill_names,
                                   ._named.descs             = (char**)vec_fill_descs,
                                   ._named.value             = 0};

static opt_option scale         = {._double.header.type       = opt_double,
                                   ._double.header.short_name = 's',
                                   ._double.header.long_name  = "scale",
                                   ._double.header.help       = "vector scale",
                                   ._double.required          = 0,
                                   ._double.min               = 0,
                                   ._double.max               = DBL_MAX,
                                   ._double.value             = 1.0};

static opt_option cond          = {._double.header.type       = opt_double,
                                   ._double.header.short_name = 'c',
                                   ._double.header.long_name  = "cond",
                                   ._double.header.help       = "condition number",
                                   ._double.required          = 0,
                                   ._double.min               = 1.0,
                                   ._double.max               = DBL_MAX,
                                   ._double.value             = 1e3};

int vecvec_show_help(void){
  opt_show_option(vec_fill_type);
  opt_show_option(scale);
  opt_show_option(cond);
  return vecvec_fill_show_help();
}

const char* vecvec_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &vec_fill_type);
  opt_eval_option(argc, argv, &scale);
  opt_eval_option(argc, argv, &cond);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%s)", vecvec_fill_name(argc, argv), vec_fill_names[vec_fill_type._named.value]);
  return name_buffer;
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY){
  int rc;

  opt_eval_option(argc, argv, &vec_fill_type);
  opt_eval_option(argc, argv, &scale);
  opt_eval_option(argc, argv, &cond);
  rc = vecvec_fill_test(argc, argv, N, incX, incY, vec_fill_type._named.value, scale._double.value, cond._double.value);
  return rc;
}
