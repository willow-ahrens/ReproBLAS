#include <stdio.h>
#include "../common/test_opt.h"
#define MAX_NAME 100

extern const char *file_name(int argc, char** argv);
extern int file_check(int argc, char** argv, char *fname);

static char namebuf[MAX_NAME];

const char* name(int argc, char** argv){
  char *fname = opt_read_string(argc, argv, "-i", "");
  snprintf(namebuf, MAX_NAME, "%s (%s)", file_name(argc, argv), fname);
  return namebuf;
}

int check(int argc, char** argv){
  char *fname = opt_read_string(argc, argv, "-i", "");
  int rc = file_check(argc, argv, fname);
  return rc;
}
