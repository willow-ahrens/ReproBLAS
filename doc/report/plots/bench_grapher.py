import matplotlib
matplotlib.use("PDF")
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

#Easy Vs Hard Sums
n_groups = 2

#fig, ax = plt.subplots()

index = np.arange(2)
bar_width = 0.35

opacity = 0.4

rects1 = plt.bar(index, np.array(tables[0][1][-2:]), bar_width,
                 alpha=opacity,
                 hatch='//',
                 color='b',
                 label='Uniform Random [0, 1) Distribution')

rects2 = plt.bar(index + bar_width, np.array(tables[0][2][-2:]), bar_width,
                 alpha=opacity,
                 hatch='',
                 color='r',
                 label='Full Exponent Range Distribution')

plt.xlabel('Reproducible Summation Function')
plt.ylabel('Time (seconds)')
plt.title('Reproducible Summation Time Vs. Distribution ({} Values)'.format(int(tables[0][1][0])))
plt.xticks(index + bar_width, ('rdsum (double)', 'rssum (float)'))
plt.legend(loc="lower right")

plt.tight_layout()
plt.savefig("easy_vs_hard")
plt.close("all")

t = 1

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
plt.title("Summation Time (Normalized To Theoretical Peak) Vs. N")
plt.xlabel("N (Vector Size N)")
plt.ylabel("Time (Normalized to Theoretical Recursive Summation Peak")
plt.legend(loc="best")
plt.savefig("sum_comparison")
plt.close("all")

t = 3

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
plt.title("Dot Product Time (Normalized To Theoretical Peak) Vs. N")
plt.xlabel("N (Vector Size N)")
plt.ylabel("Time (Normalized to Theoretical Dot Product Peak")
plt.legend(loc="best")
plt.savefig("dot_comparison")
plt.close("all")

t = 5

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
plt.title("Square Matrix-Vector Product Time (Normalized To Theoretical Peak) Vs. N")
plt.xlabel("N (Matrix Size NxN)")
plt.ylabel("Time (Normalized to Theoretical Matrix-Vector Product Peak")
plt.legend(loc="best")
plt.savefig("gemv_comparison")
plt.close("all")

t = 7

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
plt.title("Square Transposed Matrix-Vector Product Time (Normalized To Theoretical Peak) Vs. N")
plt.xlabel("N (Matrix Size NxN)")
plt.ylabel("Time (Normalized to Theoretical Matrix-Vector Product Peak")
plt.legend(loc="best")
plt.savefig("gemv_trans_comparison")
plt.close("all")

t = 9

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
plt.title("Square Matrix-Matrix Product Time (Normalized To Theoretical Peak) Vs. N")
plt.xlabel("N (Matrix Size NxN)")
plt.ylabel("Time (Normalized to Theoretical Matrix-Matrix Product Peak")
plt.legend(loc="best")
plt.savefig("gemm_comparison")
plt.close("all")

f.close()
