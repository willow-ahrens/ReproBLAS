#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_zdotu_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  reproBLAS_rzdotu_sub(DIDEFAULTFOLD, N, X, incX, Y, incY, dotu);
}
