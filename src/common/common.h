/*
 *  Created   13/10/25   H.D. Nguyen
 */

#ifndef COMMON_H_
#define COMMON_H_

#ifndef MAX
#define MAX(A,B) ((A)>(B)?(A):(B))
#endif

#ifndef MIN
#define MIN(A,B) ((A)<(B)?(A):(B))
#endif

typedef union long_double_ {
  double d;
  long   l;
} long_double;

typedef union int_float_ {
  float f;
  int   i;
} int_float;

#endif
