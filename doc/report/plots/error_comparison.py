import matplotlib
matplotlib.use("svg")
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

def one(n):
  return 1.0

relerrs = [standard_relerr(1000), standard_relerr(100), indexed_naive_relerr, indexed_careful_relerr, one]
linestyles = ["-", "-", "-", "-", "--"]
names = ["Standard Sum n=1000", "Standard Sum n=100", "Indexed Sum (Naive) n=100,1000", "Indexed Sum (Careful) n=100,1000", "Relative Error 1"]
conds = [10.0**e for e in samples(0, 22, 10000)]
for (relerr, name, style) in zip(relerrs, names, linestyles):
  plot.plot(conds, [relerr(cond) for cond in conds], label="{}".format(name), linestyle=style)

plot.xscale("log")
plot.yscale("log")
plot.xlabel("Condition Number")
plot.ylabel("Relative Error Bound")
plot.legend(loc='upper left')
plot.savefig("error_comparison.svg", format='svg', dpi=1200)
