#include <indexedBLAS.h>
#include <idxd.h>
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

  double_indexed *above = idxd_dialloc(2);
  double_indexed *at = idxd_dialloc(2);
  double_indexed *below = idxd_dialloc(2);

  idxd_didconv(2, DIGRANULARITY * 2, above);
  idxd_didconv(2, DIGRANULARITY, at);
  idxd_didconv(2, DIGRANULARITY * 0.5, below);

  if(idxd_ddiconv(2, above) == 0.0){
    printf("idxd_ddiconv(idxd_didconv(DIGRANULARITY * 2)) == 0\n");
    printf("idxd_ddiconv(idxd_didconv(%g)) == 0\n", DIGRANULARITY * 2);
    printf("idxd_ddiconv(idxd_didconv(2^%g)) == 0\n", log2(DIGRANULARITY * 2));
    idxd_diprint(2, above);
    return 1;
  }

  if(idxd_ddiconv(2, at) == 0.0){
    printf("idxd_ddiconv(idxd_didconv(DIGRANULARITY)) == 0\n");
    printf("idxd_ddiconv(idxd_didconv(%g)) == 0\n", DIGRANULARITY);
    printf("idxd_ddiconv(idxd_didconv(2^%g)) == 0\n", log2(DIGRANULARITY));
    idxd_diprint(2, at);
    return 1;
  }

  if(idxd_ddiconv(2, below) != 0.0){
    printf("idxd_ddiconv(idxd_didconv(DIGRANULARITY * 0.5)) != 0\n");
    printf("idxd_ddiconv(idxd_didconv(%g)) != 0\n", DIGRANULARITY * 0.5);
    printf("idxd_ddiconv(idxd_didconv(2^%g)) != 0\n", log2(DIGRANULARITY * 0.5));
    idxd_diprint(2, below);
    return 1;
  }
  return 0;
}
