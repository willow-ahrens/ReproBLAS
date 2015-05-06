#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"

int main(int argc, char **argv){
  (void)argc;
  (void)argv;
  float_indexed *foo = sialloc(3);
  sisetzero(3, foo);
  float bar = ldexp(0.5, FLT_MAX_EXP);
  printf("bar %g\n", bar);
  sisadd(3, bar, foo);
  printf("foo %g\n", ssiconv(3, foo));
  return 0;
}
