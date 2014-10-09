#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "test_vec.h"
#include "test_opt.h"
#include "test_file.h"

#define MAX_LINE 128
#define MAX_PATH 1024
#define MAX_FILE 1024
int main (int argc, char** args) {
  if (argc < 2 || opt_find_option(argc, args, "--help") >= 0) {
    int i;
    char line[MAX_LINE];
    fprintf(stdout, "Usage: %s [options] \n", args[0]);
    opt_show_option("-N", "vector size[default 1000]");
    opt_show_option("-o", "output file name");
    opt_show_option("-d", "data type  [default 0] (0: double,");
    opt_show_option_extended("                        1: float,");
    opt_show_option_extended("                        2: complex double,");
    opt_show_option_extended("                        3: complex float)");
    snprintf(line, MAX_LINE, "fill [default 1] (0: %s,", vec_fill_name(0));
    opt_show_option("-t", line);
    for(i = 1; i < vec_fill_MAX - 1; i++){
      snprintf(line, MAX_LINE, "                %3d: %s,", i, vec_fill_name(i));
      opt_show_option_extended(line);
    }
    snprintf(line, MAX_LINE, "                %3d: %s)", i, vec_fill_name(i));
    opt_show_option_extended(line);
    opt_show_option("-K", "Condition number [default: 1e3]");
    opt_show_option("-s", "scale [default: 1]");
    return 0;
  }

  int fill     = opt_read_int(argc, args, "-t", 1);
  int N        = opt_read_int(argc, args, "-N", 1000);
  int type     = opt_read_int(argc, args, "-d", 0);
  double K     = opt_read_double(argc, args, "-K", 1e3);
  double scale = opt_read_double(argc, args, "-s", 1.0);
  char fname[MAX_FILE];
  char type_name_char;
  int i;
  switch(type){
    case 0: type_name_char = 'd'; break;
    case 1: type_name_char = 's'; break;
    case 2: type_name_char = 'z'; break;
    case 3: type_name_char = 'c'; break;
    default:
      fprintf(stderr, "error: invalid datatype\n");
      exit(125);
  }
  snprintf(fname, MAX_FILE, "%c_%s_N%d", type_name_char, vec_fill_name(fill), N);
  for(i = 0; fname[i] != '\0'; i++){
    switch(fname[i]){
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
      fname[i] = '_';
    }
  }
  strncat(fname, ".dat", MAX_FILE);
  strncpy(fname, opt_read_string(argc, args, "-o", fname), MAX_FILE);

  if (file_exists(fname)) {
    fprintf(stderr, "error: file %s already existed, please choose a different name.\n", fname);
    exit(1);
  }

  vec_random_seed();

  switch(type){
    case 0: {
      double* data = dvec_alloc(N, 1);
      dvec_fill(N, data, 1, fill, scale, K);
      int i;
      file_write_vector(fname, N, data, sizeof(double));
      free(data);
      break;
    }
    case 1: {
      float* data = svec_alloc(N, 1);
      svec_fill(N, data, 1, fill, scale, K);
      file_write_vector(fname, N, data, sizeof(float));
      free(data);
      break;
    }
    case 2: {
      double complex* data = zvec_alloc(N, 1);
      zvec_fill(N, data, 1, fill, scale, K);
      file_write_vector(fname, N, data, sizeof(double complex));
      free(data);
      break;
    }
    case 3: {
      float complex* data = cvec_alloc(N, 1);
      cvec_fill(N, data, 1, fill, scale, K);
      file_write_vector(fname, N, data, sizeof(float complex));
      free(data);
      break;
    }
    default:
      fprintf(stderr, "error: invalid datatype\n");
      exit(125);
  }

  return 0;
}
