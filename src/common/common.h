/*
 *  Created   13/10/25   H.D. Nguyen
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <stdint.h>
#include <float.h>

#ifndef MAX
#define MAX(A,B) ((A)>(B)?(A):(B))
#endif

#ifndef MIN
#define MIN(A,B) ((A)<(B)?(A):(B))
#endif

typedef union long_double_ {
  double   d;
  uint64_t l;
} long_double;

typedef union int_float_ {
  float    f;
  uint32_t i;
} int_float;

inline float UFP(double X) {
  long_double tmp_UFP;
  tmp_UFP.d = X;
  tmp_UFP.l &= (2ull * DBL_MAX_EXP - 1) << (DBL_MANT_DIG - 1);
  return tmp_UFP.d;
}

inline float UFPF(float X) {
  int_float tmp_UFPF;
  tmp_UFPF.f = X;
  tmp_UFPF.i &= (2ul * FLT_MAX_EXP - 1) << (FLT_MANT_DIG - 1);
  return tmp_UFPF.f;
}

inline int EXP(double X) {
  long_double tmp_EXP;
  tmp_EXP.d = X;
  return (tmp_EXP.l >> (DBL_MANT_DIG - 1)) & (2 * DBL_MAX_EXP - 1);
}

#define EXP_BIAS (DBL_MAX_EXP - 2)

inline int EXPF(float X) {
  int_float tmp_EXPF;
  tmp_EXPF.f = X;
  return (tmp_EXPF.i >> (FLT_MANT_DIG - 1)) & (2 * FLT_MAX_EXP - 1);
}

#define EXPF_BIAS (FLT_MAX_EXP - 2)

inline int ISNANINF(double X) {
  long_double tmp_ISNANINF;
  tmp_ISNANINF.d = X;
  return (tmp_ISNANINF.l & ((2ull * DBL_MAX_EXP - 1) << (DBL_MANT_DIG - 1))) == ((2ull * DBL_MAX_EXP - 1) << (DBL_MANT_DIG - 1));
}

inline int ISNANINFF(float X) {
  int_float tmp_ISNANINFF;
  tmp_ISNANINFF.f = X;
  return (tmp_ISNANINFF.i & ((2ul * FLT_MAX_EXP - 1) << (FLT_MANT_DIG - 1))) == ((2ul * FLT_MAX_EXP - 1) << (FLT_MANT_DIG - 1));
}

#endif
