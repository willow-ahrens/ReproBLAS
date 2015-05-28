import multiprocessing
import subprocess
import sys
import time

def status(i, n):
  i += 1
  width = 80
  done = (i * (width - 2))//n
  remaining = width - 2 - done
  sys.stdout.write("\r[{}{}] {:7.3f}%".format(done * "#", remaining * " ", (100.0*i)/n))
  sys.stdout.flush()

def execute(command, verbose):
  if verbose == "true":
    print(command)
  rc = 0
  try:
    out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True).decode(sys.stdout.encoding)
  except subprocess.CalledProcessError as e:
    rc = e.returncode
    out = e.output.decode(sys.stdout.encoding)
  if verbose == "true":
    print(out)
  return (rc, out)

def run(command_list, verbose="false"):
  """
  A function that runs the command list (a list of string commands) on the
  target and returns a list of their results as a list of tuples of (return
  code, output)
  """
  result = []
  for (i, command) in enumerate(command_list):
    result.append(execute(command, verbose))
    if verbose != "true":
      status(i, len(command_list))
  if verbose != "true":
    print()
  return result

def execute_parallel(command_completed_verbose):
  command, completed, verbose = command_completed_verbose
  (rc, out) = execute(command, verbose)
  completed.put(command)
  return (rc, out)

def run_parallel(command_list, verbose="false"):
  """
  A function that runs the command list (a list of string commands) on the
  target (possibly in parallel) and returns a list of their results as a list
  of tuples of (return code, output)
  """
  p = multiprocessing.Pool(multiprocessing.cpu_count())
  m = multiprocessing.Manager()
  completed = m.Queue()
  result = p.map_async(execute_parallel, [(command, completed, verbose) for command in command_list])
  while not result.ready():
    if verbose != "true":
      status(completed.qsize(), len(command_list))
    time.sleep(1);
  if verbose != "true":
    print()
  return list(result.get())

def peak_time(data):
  """
    data is a dictionary containing the following keys:
      d_add - number of double precision additions
      d_mul - number of double precision multiplications
      d_fma - number of single precision fused multiplication and additions
      d_cmp - number of double precision comparisons
      d_orb  - number of double precision bitwise or
      s_add - number of single precision additions
      s_mul - number of single precision multiplications
      s_fma - number of single precision fused multiplication and additions
      s_cmp - number of single precision comparisons
      s_orb  - number of single precision bitwise or
      freq  - frequency of cpu
      vec   - best vectorization available ("AVX", "SSE", "SISD")
    peak_time(data) returns the theoretical best time in which the cpu could
    complete the given instructions (in any order)
  """
  if data["vec"] == "SISD":
    vec_d_ops = 1
    vec_s_ops = 1
  elif data["vec"] == "SSE":
    vec_d_ops = 2
    vec_s_ops = 4
  elif data["vec"] == "AVX":
    vec_d_ops = 4
    vec_s_ops = 8
  d_ops = data["d_add"] + data["d_mul"] + data["d_fma"] + data["d_orb"]
  s_ops = data["s_add"] + data["s_mul"] + data["s_fma"] + data["s_orb"]
#  d_ops = max(data["d_add"] + data["d_mul"] + data["d_fma"], data["d_orb"])
#  s_ops = max(data["s_add"] + data["s_mul"] + data["s_fma"], data["s_orb"])
  print(d_ops)
  return float(d_ops/vec_d_ops + s_ops/vec_s_ops)/data["freq"];

version = "0.0.0"
