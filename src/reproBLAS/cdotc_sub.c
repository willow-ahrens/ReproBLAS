#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_cdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  reproBLAS_rcdotc_sub(SIDEFAULTFOLD, N, X, incX, Y, incY, dotc);
}
