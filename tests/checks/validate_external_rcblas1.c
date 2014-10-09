#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/test_vec.h"
#include "../common/test_opt.h"
#include "../common/test_file.h"
#include "rcblas1_wrapper.h"

#define MAX_FILE 1024
#define MAX_NAME 100

char namebuf[MAX_NAME];

const char* file_name(int argc, char** argv) {
  int func = opt_read_int(argc, argv, "-f", 0);
  snprintf(namebuf, MAX_NAME, "Validate %s external", wrap_rcblas1_name(func));
  return namebuf;
}

int file_check(int argc, char** argv, char *fname) {
  int i;
  int N;
  char ref_fname[MAX_FILE];
  char Iref_fname[MAX_FILE];
  int func = opt_read_int(argc, argv, "-f", 0);

  float complex *x;
  float complex *y;

  float complex ref;
  I_float_Complex Iref;

  float complex res;
  I_float_Complex Ires;

  file_read_vector(fname, &N, (void**)&x, sizeof(float complex));
  y = cvec_alloc(N, 1);
  //fill y with -i where necessary
  cvec_fill(N, y, 1, vec_fill_CONSTANT, -_Complex_I, 1.0);

  ((char*)file_ext(fname))[0] = '\0';
  snprintf(ref_fname, MAX_FILE, "%s__%s.dat", fname, wrap_rcblas1_name(func));
  snprintf(Iref_fname, MAX_FILE, "%s__I%s.dat", fname, wrap_rcblas1_name(func));

  res = (wrap_rcblas1_func(func))(N, x, 1, y, 1);
  Ires = (wrap_Icblas1_func(func))(N, x, 1, y, 1);

  if(opt_find_option(argc, argv, "-r") >= 0){
    ref = res;
    Iref = Ires;

    file_write_vector(ref_fname, 1, &ref, sizeof(ref));
    file_write_vector(Iref_fname, 1, &Iref, sizeof(Iref));
  } else {
    void *data;
    int unused0;
    file_read_vector(ref_fname, &unused0, &data, sizeof(ref));
    ref = *(float complex*)data;
    free(data);
    file_read_vector(Iref_fname, &unused0, &data, sizeof(Iref));
    Iref = *(I_float_Complex*)data;
    free(data);
    if(ref != res){
      printf("%s(%s) = %g + %gi != %g + %gi\n", wrap_rcblas1_name(func), fname, CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
      return 1;
    }
    if(memcmp(&Iref, &Ires, sizeof(Iref)) != 0){
      printf("I%s(%s) = %g + %gi != %g + %gi\n", wrap_rcblas1_name(func), fname, CREAL_(Iconv2c(Ires)), CIMAG_(Iconv2c(Ires)), CREAL_(Iconv2c(Iref)), CIMAG_(Iconv2c(Iref)));
      printf("Ref I_float_Complex:\n");
      cIprint(Iref);
      printf("\nRes I_float_Complex:\n");
      cIprint(Ires);
      printf("\n");
      return 1;
    }
  }

  free(x);
  free(y);
  return 0;
}
