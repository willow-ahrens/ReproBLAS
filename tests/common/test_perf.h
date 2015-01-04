#ifndef __TEST_PERF_H
#define __TEST_PERF_H

#define perf_unit_HERTZ  0
#define perf_unit_FLOPS  1
#define perf_unit_PEAK   2
#define perf_prec_SINGLE 1
#define perf_prec_DOUBLE 2

const char* perf_unit_name(int perf_unit);

double perf_cpu_freq();

double perf_peak_flops(int perf_prec);

double perf_output(double time, int N, int trials, int flop_per_N, int unit, int perf_prec);

#endif
