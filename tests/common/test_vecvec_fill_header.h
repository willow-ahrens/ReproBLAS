#include <stdio.h>
#include "test_opt.h"
#include "test_vec.h"

#include "test_vecvec_header.h"

const char* vecvec_fill_name(int argc, char** argv);
int vecvec_fill_test(int argc, char** argv, int N, int incX, int incY, int vec_fill_type, double scale, double cond);

const char* vecvec_name(int argc, char** argv){
  opt_option vec_fill_type;
  opt_option scale;
  opt_option cond;
  static char name_buffer[MAX_LINE];

  vec_fill_type.header.type       = opt_named;
  vec_fill_type.header.short_name = 'f';
  vec_fill_type.header.long_name  = "f_type";
  vec_fill_type.header.help       = "vector fill type";
  vec_fill_type._named.required   = 0;
  vec_fill_type._named.n_names    = (int)vec_fill_n_names;
  vec_fill_type._named.names      = (char**)vec_fill_names;
  vec_fill_type._named.descs      = (char**)vec_fill_descs;
  vec_fill_type._named.value      = 0;

  scale.header.type       = opt_double;
  scale.header.short_name = 's';
  scale.header.long_name  = "scale";
  scale.header.help       = "vector scale";
  scale._double.required  = 0;
  scale._double.min       = 0;
  scale._double.max       = DBL_MAX;
  scale._double.value     = 1.0;

  cond.header.type       = opt_double;
  cond.header.short_name = 'c';
  cond.header.long_name  = "cond";
  cond.header.help       = "condition number";
  cond._double.required  = 0;
  cond._double.min       = 1.0;
  cond._double.max       = DBL_MAX;
  cond._double.value     = 1e3;

  if(help._flag.exists){
    opt_show_option(vec_fill_type);
    opt_show_option(scale);
    opt_show_option(cond);
    vecvec_fill_name(argc, argv);
  }

  opt_eval_option(argc, argv, &vec_fill_type);
  opt_eval_option(argc, argv, &scale);
  opt_eval_option(argc, argv, &cond);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%s)", vecvec_fill_name(argc, argv), vec_fill_names[vec_fill_type._named.value]);
  return name_buffer;
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY){
  opt_option vec_fill_type;
  opt_option scale;
  opt_option cond;
  int rc;

  vec_fill_type.header.type       = opt_named;
  vec_fill_type.header.short_name = 'f';
  vec_fill_type.header.long_name  = "f_type";
  vec_fill_type.header.help       = "vector fill type";
  vec_fill_type._named.required   = 0;
  vec_fill_type._named.n_names    = (int)vec_fill_n_names;
  vec_fill_type._named.names      = (char**)vec_fill_names;
  vec_fill_type._named.descs      = (char**)vec_fill_descs;
  vec_fill_type._named.value      = 0;

  scale.header.type       = opt_double;
  scale.header.short_name = 's';
  scale.header.long_name  = "scale";
  scale.header.help       = "vector scale";
  scale._double.required  = 0;
  scale._double.min       = 0;
  scale._double.max       = DBL_MAX;
  scale._double.value     = 1.0;

  cond.header.type       = opt_double;
  cond.header.short_name = 'c';
  cond.header.long_name  = "cond";
  cond.header.help       = "condition number";
  cond._double.required  = 0;
  cond._double.min       = 1.0;
  cond._double.max       = DBL_MAX;
  cond._double.value     = 1e3;

  opt_eval_option(argc, argv, &vec_fill_type);
  rc = vecvec_fill_test(argc, argv, N, incX, incY, vec_fill_type._named.value, scale._double.value, cond._double.value);
  return rc;
}
