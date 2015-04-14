import argparse
import multiprocessing
import os
import subprocess
import sys

def callsafe(command):
#  print command
  rc = 0
  try:
    out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True).decode(sys.stdout.encoding)
  except subprocess.CalledProcessError as e:
    rc = e.returncode
    out = e.output.decode(sys.stdout.encoding)
#  print out
  return (rc, out)

def call(command):
#  print command
  return subprocess.check_output(command, shell=True).decode(sys.stdout.encoding)

top = call("make top").split()[-1]

def make_clean(location):
  call("cd {0}; make clean".format(os.path.join(top, location)))

def make(executable, args = None, id = None, remake = False):
  executable_dir = os.path.join(top, os.path.split(executable)[0])
  executable_name = os.path.split(executable)[1]
  if executable_dir not in make.build_dir:
    make.build_dir[executable_dir] = call("cd {0}; make pbd".format(executable_dir)).split()[-1]
  build_dir = make.build_dir[executable_dir]
  build_name = executable_name
  if id:
    build_name = "{}__{}{}".format(os.path.splitext(executable_name)[1], id, os.path.splitext(executable_name[0]))
  build = os.path.join(build_dir, build_name)
  if not os.path.isfile(build) or remake:
    result = os.path.join(build_dir, executable_name)
    callsafe("rm {}".format(result))
    env = ""
    if args:
      env = "ARGS={}".format(args)
    callsafe("make -j {} {} {}".format(multiprocessing.cpu_count(), result, env))
    assert os.path.isfile(result), "Error: make unsuccessful."
    if id:
      call("cp {} {}".format(result, build))
    assert os.path.isfile(build), "Error: make unsuccessful."
  return build
make.build_dir = {}

def flags(params, args):
  return " ".join(['-{0} "{1}"'.format(param, arg) if len(param) == 1 else '--{0} "{1}"'.format(param, arg) for (param, arg) in zip(params, args)])

def get_vectorization():
  make("scripts/get_vectorization.c", remake = True)
  return call(make("scripts/get_vectorization", remake = True)).split()[0]

