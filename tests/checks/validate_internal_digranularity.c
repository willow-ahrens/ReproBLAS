#include <binnedBLAS.h>
#include <binned.h>
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

  double_binned *above = binned_dballoc(2);
  double_binned *at = binned_dballoc(2);
  double_binned *below = binned_dballoc(2);

  binned_dbdconv(2, DIGRANULARITY * 2, above);
  binned_dbdconv(2, DIGRANULARITY, at);
  binned_dbdconv(2, DIGRANULARITY * 0.5, below);

  if(binned_ddbconv(2, above) == 0.0){
    printf("binned_ddbconv(binned_dbdconv(DIGRANULARITY * 2)) == 0\n");
    printf("binned_ddbconv(binned_dbdconv(%g)) == 0\n", DIGRANULARITY * 2);
    printf("binned_ddbconv(binned_dbdconv(2^%g)) == 0\n", log2(DIGRANULARITY * 2));
    binned_dbprint(2, above);
    return 1;
  }

  if(binned_ddbconv(2, at) == 0.0){
    printf("binned_ddbconv(binned_dbdconv(DIGRANULARITY)) == 0\n");
    printf("binned_ddbconv(binned_dbdconv(%g)) == 0\n", DIGRANULARITY);
    printf("binned_ddbconv(binned_dbdconv(2^%g)) == 0\n", log2(DIGRANULARITY));
    binned_dbprint(2, at);
    return 1;
  }

  if(binned_ddbconv(2, below) != 0.0){
    printf("binned_ddbconv(binned_dbdconv(DIGRANULARITY * 0.5)) != 0\n");
    printf("binned_ddbconv(binned_dbdconv(%g)) != 0\n", DIGRANULARITY * 0.5);
    printf("binned_ddbconv(binned_dbdconv(2^%g)) != 0\n", log2(DIGRANULARITY * 0.5));
    binned_dbprint(2, below);
    return 1;
  }
  return 0;
}
