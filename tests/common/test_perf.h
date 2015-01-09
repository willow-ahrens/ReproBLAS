#ifndef __TEST_PERF_H
#define __TEST_PERF_H

#define perf_prec_SINGLE 1
#define perf_prec_DOUBLE 2

double perf_cpu_freq();

double perf_peak_flops(int perf_prec);

const char* perf_vec();

void perf_output_perf(double time, int n_elements, int trials);

void perf_output_desc(int n_adds, int n_muls, int n_or_bits, int perf_prec);

#endif
