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

#define UFP(X) ({ \
  long_double tmp_UFP; \
  tmp_UFP.d = X; \
  tmp_UFP.l &= (2ull * DBL_MAX_EXP - 1) << (DBL_MANT_DIG - 1); \
  tmp_UFP.d; \
  })

#define UFPF(X) ({ \
  int_float tmp_UFPF; \
  tmp_UFPF.f = X; \
  tmp_UFPF.i &= (2ul * FLT_MAX_EXP - 1) << (FLT_MANT_DIG - 1); \
  tmp_UFPF.f; \
  })

#define EXP(X) ({ \
  long_double tmp_EXP; \
  tmp_EXP.d = X; \
  (int)((tmp_EXP.l >> (DBL_MANT_DIG - 1)) & (2 * DBL_MAX_EXP - 1));\
  })

#define EXP_BIAS (DBL_MAX_EXP - 2)

#define EXPF(X) ({ \
  int_float tmp_EXPF; \
  tmp_EXPF.f = X; \
  (int)((tmp_EXPF.i >> (FLT_MANT_DIG - 1)) & (2 * FLT_MAX_EXP - 1));\
  })

#define EXPF_BIAS (FLT_MAX_EXP - 2)

#define ISNANINF(X) ({ \
  long_double tmp_ISNANINF; \
  tmp_ISNANINF.d = X; \
  (tmp_ISNANINF.l & ((2ull * DBL_MAX_EXP - 1) << (DBL_MANT_DIG - 1))) == ((2ull * DBL_MAX_EXP - 1) << (DBL_MANT_DIG - 1));\
  })

#define ISNANINFF(X) ({ \
  int_float tmp_ISNANINFF; \
  tmp_ISNANINFF.f = X; \
  (tmp_ISNANINFF.i & ((2ul * FLT_MAX_EXP - 1) << (FLT_MANT_DIG - 1))) == ((2ul * FLT_MAX_EXP - 1) << (FLT_MANT_DIG - 1)); \
  })

#endif
