import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--bench_data", default="bench_data.csv", help="input data file path")
args = parser.parse_args()

f = open(args.bench_data)
r = csv.reader(f)
tables = [[]]
for row in r:
  row = [entry.strip(" ()") for entry in row]
  newrow = []
  for entry in row:
    try:
      newrow.append(float(entry.strip()))
    except ValueError:
      newrow.append(entry.strip())
  row = newrow
  if row:
    tables[-1].append(row)
  else:
    tables.append([])
tables = tables[:-1]

plt.close("all")

#Easy Vs Hard Sums
n_groups = 2

index = np.arange(2)
bar_width = 0.35

opacity = 0.8

rects1 = plt.bar(index, np.array(tables[1][1][-2:]) * 1000000.0, bar_width,
                 alpha=opacity,
                 color='red',
                 label='Uniform Random [0, 1) Distribution')

rects2 = plt.bar(index + bar_width, np.array(tables[1][2][-2:]) * 1000000.0, bar_width,
                 alpha=opacity,
                 hatch='//',
                 color='blue',
                 label='Half Exponent Range Distribution')

plt.xlabel('Reproducible Summation Function')
plt.ylabel('Time (microseconds)')
plt.title('Reproducible Summation Time Vs. Distribution ({} Values)'.format(int(tables[1][1][0])))
plt.xticks(index + bar_width, ('rdsum (double)', 'rssum (float)'))
plt.legend(loc="upper right")

plt.tight_layout()
plt.savefig("easy_vs_hard.eps", format='eps', dpi=1200)
plt.close("all")

t = 2

N = np.array([row[0] for row in tables[t][1:]])
blas_sum_peak = np.array([row[-2] for row in tables[t][1:]])
rblas_sum_peak = np.array([row[-1] for row in tables[t][1:]])/blas_sum_peak
blas_sum = np.array([row[-2] for row in tables[t + 1][1:]])/blas_sum_peak
rblas_sum = np.array([row[-1] for row in tables[t + 1][1:]])/blas_sum_peak
blas_sum_peak /= blas_sum_peak
plt.plot(N, blas_sum_peak, ':', label = "for loop theoretical peak")
plt.plot(N, rblas_sum_peak, '-.', label = "rdsum theoretical peak")
plt.plot(N, blas_sum, '--', label = "for loop")
plt.plot(N, rblas_sum, '-', label = "rdsum")
plt.title("Summation Time (Normalized To Peak) Vs. N")
plt.xlabel("N (Vector Size N)")
plt.ylabel("Time (Normalized to Theoretical Recursive Summation Peak)")
plt.legend(loc="best")
plt.savefig("sum_comparison.eps", format='eps', dpi=1200)
plt.close("all")

t = 4

N = np.array([row[0] for row in tables[t][1:]])
blas_sum_peak = np.array([row[-2] for row in tables[t][1:]])
rblas_sum_peak = np.array([row[-1] for row in tables[t][1:]])/blas_sum_peak
blas_sum = np.array([row[-2] for row in tables[t + 1][1:]])/blas_sum_peak
rblas_sum = np.array([row[-1] for row in tables[t + 1][1:]])/blas_sum_peak
blas_sum_peak /= blas_sum_peak
plt.plot(N, blas_sum_peak, ':', label = "ddot theoretical peak")
plt.plot(N, rblas_sum_peak, '-.', label = "rddot theoretical peak")
plt.plot(N, blas_sum, '--', label = "ddot")
plt.plot(N, rblas_sum, '-', label = "rddot")
plt.title("Dot Product Time (Normalized To Peak) Vs. N")
plt.xlabel("N (Vector Size N)")
plt.ylabel("Time (Normalized to Theoretical Dot Product Peak)")
plt.legend(loc="best")
plt.savefig("dot_comparison.eps", format='eps', dpi=1200)
plt.close("all")

t = 6

N = np.array([row[0] for row in tables[t][1:]])
blas_sum_peak = np.array([row[-2] for row in tables[t][1:]])
rblas_sum_peak = np.array([row[-1] for row in tables[t][1:]])/blas_sum_peak
blas_sum = np.array([row[-2] for row in tables[t + 1][1:]])/blas_sum_peak
rblas_sum = np.array([row[-1] for row in tables[t + 1][1:]])/blas_sum_peak
blas_sum_peak /= blas_sum_peak
plt.plot(N, blas_sum_peak, ':', label = "dgemv theoretical peak")
plt.plot(N, rblas_sum_peak, '-.', label = "rdgemv theoretical peak")
plt.plot(N, blas_sum, '--', label = "dgemv")
plt.plot(N, rblas_sum, '-', label = "rdgemv")
plt.title("Matrix-Vector Product Time (Normalized To Peak) Vs. N")
plt.xlabel("N (Matrix Size NxN)")
plt.ylabel("Time (Normalized to Theoretical Matrix-Vector Product Peak)")
plt.legend(loc="best")
plt.savefig("gemv_comparison.eps", format='eps', dpi=1200)
plt.close("all")

t = 8

N = np.array([row[0] for row in tables[t][1:]])
blas_sum_peak = np.array([row[-2] for row in tables[t][1:]])
rblas_sum_peak = np.array([row[-1] for row in tables[t][1:]])/blas_sum_peak
blas_sum = np.array([row[-2] for row in tables[t + 1][1:]])/blas_sum_peak
rblas_sum = np.array([row[-1] for row in tables[t + 1][1:]])/blas_sum_peak
blas_sum_peak /= blas_sum_peak
plt.plot(N, blas_sum_peak, ':', label = "dgemv theoretical peak")
plt.plot(N, rblas_sum_peak, '-.', label = "rdgemv theoretical peak")
plt.plot(N, blas_sum, '--', label = "dgemv")
plt.plot(N, rblas_sum, '-', label = "rdgemv")
plt.title("Transposed Matrix-Vector Product Time (Normalized To Peak) Vs. N")
plt.xlabel("N (Matrix Size NxN)")
plt.ylabel("Time (Normalized to Theoretical Matrix-Vector Product Peak)")
plt.legend(loc="best")
plt.savefig("gemv_trans_comparison.eps", format='eps', dpi=1200)
plt.close("all")

t = 10

N = np.array([row[0] for row in tables[t][1:]])
blas_sum_peak = np.array([row[-2] for row in tables[t][1:]])
rblas_sum_peak = np.array([row[-1] for row in tables[t][1:]])/blas_sum_peak
blas_sum = np.array([row[-2] for row in tables[t + 1][1:]])/blas_sum_peak
rblas_sum = np.array([row[-1] for row in tables[t + 1][1:]])/blas_sum_peak
blas_sum_peak /= blas_sum_peak
plt.plot(N, blas_sum_peak, ':', label = "dgemm theoretical peak")
plt.plot(N, rblas_sum_peak, '-.', label = "rdgemm theoretical peak")
plt.plot(N, blas_sum, '--', label = "dgemm")
plt.plot(N, rblas_sum, '-', label = "rdgemm")
plt.title("Matrix-Matrix Product Time (Normalized To Peak) Vs. N")
plt.xlabel("N (Matrix Size NxN)")
plt.ylabel("Time (Normalized to Theoretical Matrix-Matrix Product Peak)")
plt.legend(loc="best")
plt.savefig("gemm_comparison.eps", format='eps', dpi=1200)
plt.close("all")

f.close()
