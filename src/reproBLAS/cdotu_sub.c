#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_cdotu_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  reproBLAS_rcdotu_sub(SIDEFAULTFOLD, N, X, incX, Y, incY, dotu);
}
