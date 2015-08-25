#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_zdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  reproBLAS_rzdotc_sub(DIDEFAULTFOLD, N, X, incX, Y, incY, dotc);
}
