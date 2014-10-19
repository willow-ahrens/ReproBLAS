#include <stdio.h>
#include "../common/test_opt.h"
#include "test.h"
#define MAX_NAME 100

const char *file_name(int argc, char** argv);
int file_test(int argc, char** argv, char *fname);

const char* name(int argc, char** argv){
  static char namebuf[MAX_NAME];
  char *fname = opt_read_string(argc, argv, "-i", "");
  snprintf(namebuf, MAX_NAME, "%s (%s)", file_name(argc, argv), fname);
  return namebuf;
}

int test(int argc, char** argv){
  char *fname = opt_read_string(argc, argv, "-i", "");
  int rc = file_test(argc, argv, fname);
  return rc;
}
