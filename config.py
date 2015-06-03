##
#  @file  config.py
#  @brief config.py is the configuration file for python code in ReproBLAS.
#
#  Feel free to modify anything here, but be sure to read the comments describing necessary functionality
#

import itertools
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

def execute(command_verbose):
  (command, verbose) = command_verbose
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

##
#  @brief Run commands on the target machine
#
#  A function that runs the command list (a list of string commands) on the target and returns a list of their results
#
#  @param command_list list of string commands to be run on the target machine. It may be of use to know that these string commands are always in the form @c"./<executable> <args>", where @c<executable> is the absolute path to the executable and @c<args> are the arguments to be passed to the executable.
#  @param verbose be verbose if this string is equal to @c"true"
#  @return a list of tuples such that the ith tuple is @c(<rc>, <output>) where @c<rc> is the return code and @c<output> is the output of running the ith command on the host machine, respectively
#
#  @author Peter Ahrens
#  @date   28 May 2015
#
def run(command_list, verbose="false"):
  result_list = []
  if verbose != "true":
    status(0, len(command_list));
  for (i, result) in enumerate(itertools.imap(execute, itertools.izip(command_list, itertools.repeat(verbose)))):
    if verbose != "true":
      status(i, len(command_list));
    result_list.append(result)
  if verbose != "true":
    print("")
  return result_list

##
#  @brief Run commands on the target machine
#
#  A function that runs the command list (a list of string commands) on the target (optionally in parallel) and returns a list of their results
#
#  @param command_list list of string commands to be run on the target machine. It may be of use to know that these string commands are always in the form @c"./<executable> <args>", where @c<executable> is the absolute path to the executable and @c<args> are the arguments to be passed to the executable.
#  @param verbose be verbose if this string is equal to @c"true"
#  @return a list of tuples such that the ith tuple is @c(<rc>, <output>) where @c<rc> is the return code and @c<output> is the output of running the ith command on the host machine, respectively
#
#  @author Peter Ahrens
#  @date   28 May 2015
#
def run_parallel(command_list, verbose="false"):
  p = multiprocessing.Pool(multiprocessing.cpu_count())
  result_list = []
  if verbose != "true":
    status(0, len(command_list));
  for (i, result) in enumerate(p.imap(execute, [(command, verbose) for command in command_list], chunksize=multiprocessing.cpu_count() * 16)):
    if verbose != "true":
      status(i, len(command_list));
    result_list.append(result)
  if verbose != "true":
    print()
  return result_list

##
#  @brief theoretical time to execute a set of instructions sequentially on the host machine
#
#  Note that a few implementations are provided, depending on whether or not the cpu can perform fused multiply additions or process integer operations in parallel with floating point operations, etc.
#
#  @param data a dictionary containing the following keys:
#                d_add - number of double precision additions
#                d_mul - number of double precision multiplications
#                d_fma - number of single precision fused multiply additions
#                d_cmp - number of double precision comparisons
#                d_orb  - number of double precision bitwise or
#                s_add - number of single precision additions
#                s_mul - number of single precision multiplications
#                s_fma - number of single precision fused multiply additions
#                s_cmp - number of single precision comparisons
#                s_orb  - number of single precision bitwise or
#                freq  - frequency of cpu
#                vec   - best vectorization available ("AVX", "SSE", "SISD")
#  @return idealized theoretical time in which the cpu could complete the given instructions (in any order)
#
#  @author Peter Ahrens
#  @date   28 May 2015
#
def peak_time(data):
  if data["vec"] == "SISD":
    vec_d_ops = 1
    vec_s_ops = 1
  elif data["vec"] == "SSE":
    vec_d_ops = 2
    vec_s_ops = 4
  elif data["vec"] == "AVX":
    vec_d_ops = 4
    vec_s_ops = 8
#  data["d_add"] += data["d_fma"]
#  data["d_mul"] += data["d_fma"]
#  data["d_fma"] = 0
#  data["s_add"] += data["s_fma"]
#  data["s_mul"] += data["s_fma"]
#  data["s_fma"] = 0
  d_ops = data["d_add"] + data["d_mul"] + data["d_fma"] + data["d_orb"]
  s_ops = data["s_add"] + data["s_mul"] + data["s_fma"] + data["s_orb"]
#  d_ops = max(data["d_add"] + data["d_mul"] + data["d_fma"], data["d_orb"])
#  s_ops = max(data["s_add"] + data["s_mul"] + data["s_fma"], data["s_orb"])
  return float(d_ops/vec_d_ops + s_ops/vec_s_ops)/data["freq"];

##
#  @brief version number
#
#  @author Peter Ahrens
#  @date   28 May 2015
#
version = "2.0.0"
