#include <stdio.h>
#include "test_opt.h"
#include "test_vec.h"

#include "test_vecvec_header.h"

const char* vecvec_fill_name(int argc, char** argv);
int vecvec_fill_test(int argc, char** argv, int N, int incX, int incY, int vec_fill_type);

const char* vecvec_name(int argc, char** argv){
  opt_option vec_fill_type;
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

  if(help._flag.exists){
    opt_show_option(vec_fill_type);
    vecvec_fill_name(argc, argv);
  }

  opt_eval_option(argc, argv, &vec_fill_type);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%s)", vecvec_fill_name(argc, argv), vec_fill_names[vec_fill_type._named.value]);
  return name_buffer;
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY){
  opt_option vec_fill_type;
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

  opt_eval_option(argc, argv, &vec_fill_type);
  rc = vecvec_fill_test(argc, argv, N, incX, incY, vec_fill_type._named.value);
  return rc;
}
