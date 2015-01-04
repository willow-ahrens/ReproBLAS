#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "test_vec.h"
#include "test_opt.h"
#include "test_file.h"
#include "test_limits.h"

int main (int argc, char** argv) {
  opt_option help;
  opt_option N;
  opt_option output_file;
  opt_option data_type;
  opt_option vec_fill_type;
  //opt_option mat_fill_type;
  opt_option cond;
  opt_option scale;

  help.header.type       = opt_flag;
  help.header.short_name = 'h';
  help.header.long_name  = "help";
  help.header.help       = "show help";

  N.header.type       = opt_int;
  N.header.short_name = 'N';
  N.header.long_name  = "N_dim";
  N.header.help       = "N dimension size";
  N._int.required     = 0;
  N._int.min          = 0;
  N._int.max          = INT_MAX;
  N._int.value        = 2048;

  output_file.header.type       = opt_string;
  output_file.header.short_name = 'o';
  output_file.header.long_name  = "output";
  output_file.header.help       = "output file name";
  output_file._string.required  = 0;
  output_file._string.value     = "";

  data_type.header.type       = opt_named;
  data_type.header.short_name = 'd';
  data_type.header.long_name  = "d_type";
  data_type.header.help       = "data type";
  data_type._named.required   = 0;
  data_type._named.n_names    = (int)opt_data_type_n_names;
  data_type._named.names      = (char**)opt_data_type_names;
  data_type._named.descs      = (char**)opt_data_type_descs;
  data_type._named.value      = 0;

  vec_fill_type.header.type       = opt_named;
  vec_fill_type.header.short_name = 'f';
  vec_fill_type.header.long_name  = "f_type";
  vec_fill_type.header.help       = "vector fill type";
  vec_fill_type._named.required   = 0;
  vec_fill_type._named.n_names    = (int)vec_fill_n_names;
  vec_fill_type._named.names      = (char**)vec_fill_names;
  vec_fill_type._named.descs      = (char**)vec_fill_descs;
  vec_fill_type._named.value      = 0;
  /*
  mat_fill_type.header.type       = opt_named;
  mat_fill_type.header.short_name = 'g';
  mat_fill_type.header.long_name  = "g_type";
  mat_fill_type.header.help       = "matrix fill type";
  mat_fill_type._named.required   = 0;
  mat_fill_type._named.n_names    = vec_fill_type_n_names;
  mat_fill_type._named.names      = vec_fill_type_names;
  mat_fill_type._named.descs      = vec_fill_type_descs;
  */
  cond.header.type       = opt_double;
  cond.header.short_name = 'c';
  cond.header.long_name  = "cond";
  cond.header.help       = "condition number";
  cond._double.required  = 0;
  cond._double.min       = 0;
  cond._double.max       = DBL_MAX;
  cond._double.value     = 1e3;

  scale.header.type       = opt_double;
  scale.header.short_name = 's';
  scale.header.long_name  = "scale";
  scale.header.help       = "scaling factor";
  scale._double.required  = 0;
  scale._double.min       = DBL_MIN;
  scale._double.max       = DBL_MIN;
  scale._double.value     = 1;

  opt_eval_option(argc, argv, &help);
  if (help._flag.exists) {
    opt_show_option(help);
    opt_show_option(N);
    opt_show_option(output_file);
    opt_show_option(data_type);
    opt_show_option(vec_fill_type);
    //opt_show_option(mat_fill_type);
    opt_show_option(cond);
    opt_show_option(scale);
    return 0;
  }
  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &output_file);
  opt_eval_option(argc, argv, &data_type);
  opt_eval_option(argc, argv, &vec_fill_type);
  //opt_eval_option(argc, argv, &mat_fill_type);
  opt_eval_option(argc, argv, &cond);
  opt_eval_option(argc, argv, &scale);

  char output_file_name[MAX_NAME];
  char data_type_char;
  int i;
  switch(data_type._named.value){
    case 0: data_type_char = 'd'; break;
    case 1: data_type_char = 's'; break;
    case 2: data_type_char = 'z'; break;
    case 3: data_type_char = 'c'; break;
  }
  snprintf(output_file_name, MAX_NAME, "%c_%s_N%d", data_type_char, vec_fill_names[vec_fill_type._named.value], N._int.value);
  for(i = 0; output_file_name[i] != '\0'; i++){
    switch(output_file_name[i]){
      case '/':
      case '\\':
      case '?':
      case '%':
      case '*':
      case ':':
      case '|':
      case '"':
      case '<':
      case '>':
      case '.':
      case ' ':
      output_file_name[i] = '_';
    }
  }
  strncat(output_file_name, ".dat", MAX_NAME);
  if(strcmp(output_file._string.value, "") == 0){
    output_file._string.value = output_file_name;
  }

  if (file_exists(output_file._string.value)) {
    fprintf(stderr, "error: file %s already existed, please choose a different name.\n", output_file._string.value);
    exit(1);
  }

  vec_random_seed();

  switch(data_type._named.value){
    case 0: {
      double* data = dvec_alloc(N._int.value, 1);
      dvec_fill(N._int.value, data, 1, vec_fill_type._named.value, scale._double.value, cond._double.value);
      file_write_vector(output_file._string.value, N._int.value, data, sizeof(double));
      free(data);
      break;
    }
    case 1: {
      float* data = svec_alloc(N._int.value, 1);
      svec_fill(N._int.value, data, 1, vec_fill_type._named.value, scale._double.value, cond._double.value);
      file_write_vector(output_file._string.value, N._int.value, data, sizeof(float));
      free(data);
      break;
    }
    case 2: {
      double complex* data = zvec_alloc(N._int.value, 1);
      zvec_fill(N._int.value, data, 1, vec_fill_type._named.value, scale._double.value, cond._double.value);
      file_write_vector(output_file._string.value, N._int.value, data, sizeof(double complex));
      free(data);
      break;
    }
    case 3: {
      float complex* data = cvec_alloc(N._int.value, 1);
      cvec_fill(N._int.value, data, 1, vec_fill_type._named.value, scale._double.value, cond._double.value);
      file_write_vector(output_file._string.value, N._int.value, data, sizeof(float complex));
      free(data);
      break;
    }
  }
  return 0;
}
