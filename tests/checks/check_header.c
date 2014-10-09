#include <stdio.h>
#include "../common/test_opt.h"

extern const char *name(int argc, char** argv);
extern const char check(int argc, char** argv);

int main(int argc, char** argv){
  if(opt_find_option(argc, argv, "-p") >= 0){
    printf("%s\n", name(argc, argv));
  } else {
    return check(argc, argv);
  }
}
