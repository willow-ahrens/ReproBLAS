#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_csum_sub(const int N, const void* X, const int incX, void *sum) {
  reproBLAS_rcsum_sub(SIDEFAULTFOLD, N, X, incX, sum);
}
