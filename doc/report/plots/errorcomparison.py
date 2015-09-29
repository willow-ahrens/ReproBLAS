import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plot

eps = 2.0**(-53.0)
K = 3.0
W = 40.0

def samples(low, high, n):
  for i in range(n):
    yield low + i * ((low + high)/float(n))

def indexed_careful_relerr(cond):
  return 2.0**(W * (1.0 - K) - 1.0) * cond + 7 * eps

def indexed_naive_relerr(cond):
  return (2.0**(W * (1.0 - K) - 1.0) + ((2.0 * K - 1.0) * eps)) * cond

def standard_relerr(n):
  def foo(cond):
    return n * eps * cond
  return foo

relerrs = [standard_relerr(1000), standard_relerr(100), indexed_naive_relerr, indexed_careful_relerr]
names = ["Standard Summation n=1000", "Standard Summation n=100", "Indexed Summation (Naive)", "Indexed Summation (Careful)"]
conds = [10.0**e for e in samples(0, 20, 10000)]

plot.plot(conds, [1.0 for cond in conds], label="", linestyle="dotted", color="black")
for (relerr, name) in zip(relerrs, names):
  plot.plot(conds, [relerr(cond) for cond in conds], label="{}".format(name))


plot.xscale("log")
plot.yscale("log")
plot.xlabel("Condition Number")
plot.ylabel("Relative Error Bound")
plot.legend(loc='best')
plot.savefig("errorcomparison")
