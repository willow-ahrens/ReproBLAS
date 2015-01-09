#include "test_perf.h"
#include <stdio.h>
#include <stdlib.h>
#include "test_limits.h"

double perf_cpu_freq(){
  char buffer[MAX_LINE];
  buffer[MAX_LINE - 1] = '\0';
  int i = 0;
  #if defined( __MACH__ )
    FILE *f = popen("sysctl hw.cpufrequency", "r");
  #elif defined( __LINUX__ )
    FILE *f = popen("grep -m 1 ^'cpu MHz' /proc/cpuinfo", "r");
  #else
  #endif
  fgets(buffer, MAX_LINE - 1, f);
  while(buffer[i] != ':' && i < MAX_LINE - 1){
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

const char* perf_vec(){
  //TODO: use cpuid
  #if defined( __AVX__ )
    return "AVX";
  #elif defined( __SSE__ )
    return "SSE";
  #else
    return "SIMD";
  #endif
}

void perf_output_perf(double time, int n_elements, int trials){
  printf("%e\n", (double)(((long double)n_elements * trials) / time));
}

void perf_output_desc(int n_ops, char** op_names, int* op_counts){
  int i;
  printf("{");
  for(i = 0; i < n_ops; i++){
    printf("\"%s\":%d, ", op_names[i], op_counts[i]);
  }
  printf("\"cpu_freq\":%e, \"vec\":\"%s\"}\n", perf_cpu_freq(), perf_vec());
}
