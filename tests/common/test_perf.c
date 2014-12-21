#include "test_perf.h"
#include <stdio.h>
#include <stdlib.h>
#define BUFFER_SIZE 124

const char* perf_unit_name(int perf_unit){
  switch(perf_unit){
    case perf_unit_HERTZ:
      return "Hertz";
    case perf_unit_FLOPS:
      return "FLOPS";
    case perf_unit_PEAK:
      return "% Peak";
  }
  return "";
}

double perf_cpu_freq(){
  char buffer[BUFFER_SIZE];
  buffer[BUFFER_SIZE - 1] = '\0';
  int i = 0;
  #if defined( __MACH__ )
    FILE *f = popen("sysctl hw.cpufrequency", "r");
  #elif defined( __LINUX__ )
    FILE *f = popen("grep -m 1 ^'cpu MHz' /proc/cpuinfo", "r");
  #else
  #endif
  fgets(buffer, BUFFER_SIZE - 1, f);
  while(buffer[i] != ':' && i < BUFFER_SIZE - 1){
    i++;
  }
  pclose(f);
  return (double)atol(buffer + i + 1);
}

/*
ifeq ($(shell uname), Darwin)
system_info:
	@echo '  CPU FREQUENCY: ' ${patsubst hw.cpufrequency:,,$(shell sysctl hw.cpufrequency)}
	@echo '  L1 CACHE     : ' ${patsubst hw.l1dcachesize:,,$(shell sysctl hw.l1dcachesize)}
	@echo '  L2 CACHE     : ' ${patsubst hw.l2cachesize:,,$(shell sysctl hw.l2cachesize)}
	@echo '  L3 CACHE     : ' ${patsubst hw.l3cachesize:,,$(shell sysctl hw.l3cachesize)}
else
ifeq ($(shell uname), Linux)
system_info:
	@echo ' ' $(shell grep -m 1 ^'cpu MHz' /proc/cpuinfo)
	@echo ' L1 CACHE: ' $(shell cat /sys/devices/system/cpu/cpu0/cache/index1/size)
	@echo ' L2 CACHE: ' $(shell cat /sys/devices/system/cpu/cpu0/cache/index2/size)
	@echo ' L3 CACHE: ' $(shell cat /sys/devices/system/cpu/cpu0/cache/index3/size)
else
system_info:
    
  #
}
*/

double perf_peak_flops(int perf_prec){
  int simd;
  //TODO: use cpuid
  #if defined( __AVX__ )
    simd = 8 / perf_prec;
  #elif defined( __SSE__ )
    simd = 4 / perf_prec;
  #else
    simd = 1;
  #endif
  return simd * 2 * perf_cpu_freq();
}

double perf_output(double time, int N, int trials, int flop_per_N, int unit, int perf_prec){
  switch(unit){
    case perf_unit_HERTZ:
      return ((long double)N * trials) / time;
    case perf_unit_FLOPS:
      return ((long double)N * flop_per_N * trials) / time;
    case perf_unit_PEAK:
      return (((long double)N * flop_per_N * trials * 100) / time) / perf_peak_flops(perf_prec);
  }
  return 0;
}
