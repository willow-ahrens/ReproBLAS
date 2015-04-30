#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"

int main(int argc, char **argv){
  double_indexed *foo = dialloc(3);
  disetzero(3, foo);
  double bar = ldexp(-0.5, DBL_MAX_EXP);
  printf("bar %g\n", bar);
  didadd(3, bar, foo);
  //didconv(3, bar, foo);
  printf("foo %g\n", ddiconv(3, foo));
  printf("max %d %d\n", dindex(ldexp(0.5, 987)), (int)log2(dbound(dindex(1.0/0.0) + 1)/1.5));
  return 0;
}
