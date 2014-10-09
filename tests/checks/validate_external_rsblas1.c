#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/test_vec.h"
#include "../common/test_opt.h"
#include "../common/test_file.h"
#include "rsblas1_wrapper.h"

#define MAX_FILE 1024
#define MAX_NAME 100

char namebuf[MAX_NAME];

const char* file_name(int argc, char** argv) {
  int func = opt_read_int(argc, argv, "-f", 0);
  snprintf(namebuf, MAX_NAME, "Validate %s external", wrap_rsblas1_name(func));
  return namebuf;
}

int file_check(int argc, char** argv, char *fname) {
  int i;
  int N;
  char ref_fname[MAX_FILE];
  char Iref_fname[MAX_FILE];
  int func = opt_read_int(argc, argv, "-f", 0);

  float *x;
  float *y;

  float ref;
  I_float Iref;

  float res;
  I_float Ires;

  file_read_vector(fname, &N, (void**)&x, sizeof(float));
  y = svec_alloc(N, 1);
  //fill y with 1 where necessary
  svec_fill(N, y, 1, vec_fill_CONSTANT, 1.0, 1.0);

  ((char*)file_ext(fname))[0] = '\0';
  snprintf(ref_fname, MAX_FILE, "%s__%s.dat", fname, wrap_rsblas1_name(func));
  snprintf(Iref_fname, MAX_FILE, "%s__I%s.dat", fname, wrap_rsblas1_name(func));

  res = (wrap_rsblas1_func(func))(N, x, 1, y, 1);
  Ires = (wrap_Isblas1_func(func))(N, x, 1, y, 1);

  if(opt_find_option(argc, argv, "-r") >= 0){
    ref = res;
    Iref = Ires;

    file_write_vector(ref_fname, 1, &ref, sizeof(ref));
    file_write_vector(Iref_fname, 1, &Iref, sizeof(Iref));
  } else {
    void *data;
    int unused0;
    file_read_vector(ref_fname, &unused0, &data, sizeof(ref));
    ref = *(float*)data;
    free(data);
    file_read_vector(Iref_fname, &unused0, &data, sizeof(Iref));
    Iref = *(I_float*)data;
    free(data);
    if(ref != res){
      printf("%s(%s) = %g != %g\n", wrap_rsblas1_name(func), fname, res, ref);
      return 1;
    }
    if(memcmp(&Iref, &Ires, sizeof(Iref)) != 0){
      printf("I%s(%s) = %g != %g\n", wrap_rsblas1_name(func), fname, Iconv2f(Ires), Iconv2f(Iref));
      printf("Ref I_float:\n");
      sIprint(Iref);
      printf("\nRes I_float:\n");
      sIprint(Ires);
      printf("\n");
      return 1;
    }
  }

  free(x);
  free(y);
  return 0;
}
