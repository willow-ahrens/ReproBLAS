#ifndef _UNIT_WRAPPER_H
#define _UNIT_WRAPPER_H

#include <stdio.h>

#define UNIT_FLOPS         0
#define UNIT_HERTZ         1
#define UNIT_PERCENT_PEAK  2
#define PREC_SINGLE        0
#define PREC_DOUBLE        1

const char* unit_name(int unit){
  switch(unit){
    case UNIT_FLOPS:
      return "FLOPS";
    case UNIT_HERTZ:
      return "Hertz";
  }
  return "";
}

void print_performance(double time, int N, int trials, int flop_per_N, int unit, int prec){
  switch(unit){
    case UNIT_HERTZ:
      printf("%e\n", (double)((long double)N * trials)/time);
      break;
    case UNIT_FLOPS:
      printf("%e\n", (double)((long double)N * flop_per_N * trials)/time);
      break;
    case UNIT_PERCENT_PEAK:
      #if defined( __AVX__ )
        int SIMD_FACTOR = 8 / prec;
      #elif defined( __SSE__ )
        int SIMD_FACTOR = 4 / prec;
      #else
        int SIMD_FACTOR = 1;
      #endif
      printf("%e\n", (double)((long double)N * flop_per_N * trials)/time);
      break;
    default:
      printf("Error: unknown unit.\n");
      exit(1);
  }
}

#endif
