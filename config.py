import multiprocessing
import subprocess
import sys

def run(command_list):
  """
  A function that runs the command list (a list of string commands) on the
  target and returns a list of their results as a list of tuples of (return
  code, output)
  """
  def callsafe(command):
    print command
    rc = 0
    try:
      out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True).decode(sys.stdout.encoding)
    except subprocess.CalledProcessError as e:
      rc = e.returncode
      out = e.output.decode(sys.stdout.encoding)
    print out
    return (rc, out)
  return list(map(callsafe, command_list))

def run_parallel(command_list):
  """
  A function that runs the command list (a list of string commands) on the
  target (possibly in parallel) and returns a list of their results as a list 
  of tuples of (return code, output)
  """
  def callsafe(command):
    print command
    rc = 0
    try:
      out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True).decode(sys.stdout.encoding)
    except subprocess.CalledProcessError as e:
      rc = e.returncode
      out = e.output.decode(sys.stdout.encoding)
    print out
    return (rc, out)
  p = multiprocessing.Pool(multiprocessing.cpu_count())
  return p.map(callsafe, command_list)

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
  d_ops = max(data["d_add"] + data["d_mul"] + data["d_fma"], data["d_orb"])
  s_ops = max(data["s_add"] + data["s_mul"] + data["s_fma"], data["s_orb"])
  print(d_ops)
  return float(d_ops/vec_d_ops + s_ops/vec_s_ops)/data["freq"];

version = "0.0.0"
