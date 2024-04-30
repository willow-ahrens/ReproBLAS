#include <string.h>

#include <binned.h>

/**
 * @brief Set binned double precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_dbsetzero(const int fold, double_binned *X){
  memset(X, 0, binned_dbsize(fold));
}
