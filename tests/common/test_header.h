#include <stdio.h>
#include "test_opt.h"

const char *name(int argc, char** argv);
int test(int argc, char** argv);

int main(int argc, char** argv){
  if(opt_find_option(argc, argv, "-p") >= 0){
    printf("%s\n", name(argc, argv));
  } else {
    return test(argc, argv);
  }
}
