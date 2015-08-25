#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_zsum_sub(const int N, const void* X, const int incX, void *sum) {
  reproBLAS_rzsum_sub(DIDEFAULTFOLD, N, X, incX, sum);
}
