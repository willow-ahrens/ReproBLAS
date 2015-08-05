#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "../common/test_header.h"
#include "../common/test_util.h"

int show_help(void){
  return 0;
}

const char* name(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return "Validate DIGRANULARITY internally";
}

int test(int argc, char** argv) {
  (void)argc;
  (void)argv;

  util_random_seed();

  double_indexed *above = dialloc(2);
  double_indexed *at = dialloc(2);
  double_indexed *below = dialloc(2);

  didconv(2, DIGRANULARITY * 2, above);
  didconv(2, DIGRANULARITY, at);
  didconv(2, DIGRANULARITY * 0.5, below);

  if(ddiconv(2, above) == 0.0){
    printf("ddiconv(didconv(DIGRANULARITY * 2)) == 0\n");
    printf("ddiconv(didconv(%g)) == 0\n", DIGRANULARITY * 2);
    printf("ddiconv(didconv(2^%g)) == 0\n", log2(DIGRANULARITY * 2));
    diprint(2, above);
    return 1;
  }

  if(ddiconv(2, at) == 0.0){
    printf("ddiconv(didconv(DIGRANULARITY)) == 0\n");
    printf("ddiconv(didconv(%g)) == 0\n", DIGRANULARITY);
    printf("ddiconv(didconv(2^%g)) == 0\n", log2(DIGRANULARITY));
    diprint(2, at);
    return 1;
  }

  if(ddiconv(2, below) != 0.0){
    printf("ddiconv(didconv(DIGRANULARITY * 0.5)) != 0\n");
    printf("ddiconv(didconv(%g)) != 0\n", DIGRANULARITY * 0.5);
    printf("ddiconv(didconv(2^%g)) != 0\n", log2(DIGRANULARITY * 0.5));
    diprint(2, below);
    return 1;
  }
  return 0;
}
