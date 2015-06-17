/*
 *  Created   13/10/25   H.D. Nguyen
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <stdint.h>

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
  tmp_UFP.l &= 0x7FF0000000000000ull; \
  tmp_UFP.d; \
  })

#define UFPF(X) ({ \
  int_float tmp_UFPF; \
  tmp_UFPF.f = X; \
  tmp_UFPF.i &= 0x7F800000ul; \
  tmp_UFPF.f; \
  })

#define EXP(X) ({ \
  long_double tmp_EXP; \
  tmp_EXP.d = X; \
  (int)((tmp_EXP.l >> 52) & 0x7FF);\
  })

#define EXP_BIAS 0x3FE

#define EXPF(X) ({ \
  int_float tmp_EXPF; \
  tmp_EXPF.f = X; \
  (int)((tmp_EXPF.i >> 23) & 0xFF);\
  })

#define EXPF_BIAS 0x7E

#define ISNANINF(X) ({ \
  long_double tmp_ISNANINF; \
  tmp_ISNANINF.d = X; \
  tmp_ISNANINF.l & 0x7FF0000000000000ull == 0x7FF0000000000000ull;\
  })

#define ISNANINFF(X) ({ \
  int_float tmp_ISNANINFF; \
  tmp_ISNANINFF.f = X; \
  tmp_ISNANINFF.i & 0x7F800000ul == 0x7F800000ul; \
  })

#endif
