#Harness.py
#the intended test harness for all tests here.
#By peter ahrens
#comments/code niceness to come


import argparse
import os
import subprocess

OUT_LEN = 80

divider = "+" + (OUT_LEN - 2) * "-" + "+"
prev    = "feed"
built   = set()
passed  = 0
failed  = 0
not_run = 0

def callsafe(command):
  rc = 0
  try:
    out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
  except subprocess.CalledProcessError as e:
    rc = e.returncode
    out = e.output
  return (rc, out)

def call(command):
  return subprocess.check_output(command, shell=True)

def clean(test):
  call("cd {0}; make clean".format(os.path.split(test)[0]))

def run(test, args):
  global built
  test = os.path.realpath(test)
  if test not in built:
    #assert not os.path.isfile(test), "Error: make clean unsuccessful."
    call("cd {0}; make {1}".format(os.path.split(test)[0], os.path.split(test)[1]))
    assert os.path.isfile(test), "Error: make unsuccessful."
    built.add(test)
  return callsafe("{0} {1}".format(test, args))

def settings(params):
  if params:
    for value in params[0][1]:
      for setting in settings(params[1:]):
        yield [(params[0][0], value)] + setting 
  else:
    yield []

def flags(setting):
  return " ".join(['-{0} "{1}"'.format(*s) for s in setting])

def benchmark(tests, params):
  global divider
  global prev
  global names
  BENCH_LEN = 8

  if prev != "feed":
    print(divider)

  param_result  = ""
  param_divider = ""
  test_result   = ""
  test_divider  = ""
  tests_per_row = 0
  for param in params:
    param_result  += "|{{: <{0}}}".format(len(str(max(param[1]))))
    param_divider += "+" + ("-" * len(str(max(param[1]))))
  while tests_per_row < len(tests) and len(param_divider + "+" + test_divider) + BENCH_LEN + 1 < OUT_LEN: 
    test_result   += "{{: <{0}}}|".format(BENCH_LEN)
    test_divider  += ("-" * BENCH_LEN) + "+"
    tests_per_row += 1
  middle_result = ((OUT_LEN - len(param_divider + "|" + test_divider)) * " ") + "|"
  middle_divider = ((OUT_LEN - len(param_divider + "+" + test_divider)) * "-") + "+"
  result = param_result + middle_result + test_result 
  divider = param_divider + middle_divider + test_divider

  assert tests_per_row >= 1, "Error: params too large."
  for test in tests:
    assert len(names[test]) < BENCH_LEN, "Error: test name too long."

def script(tests, params):
  global divider
  global prev

  result = "|{{0: <{0}}}|".format(OUT_LEN - 2)
  if prev != "script":
    print(divider)
  for test in tests:
    for setting in settings(params):
      name = "{0} {1}".format(test, flags(setting))
      print(result.format(name))
      (rc, out) = run(test, flags(setting))
      if rc != 0:
        for line in out.split("\n"):
          while line:
            print(result.format(line[:OUT_LEN - 2]))
            line = line[OUT_LEN - 2:]
  prev = "script"

def check(tests, params):
  global divider
  global prev
  global names
  global passed
  global failed
  global not_run

  assert tests, "Error: no tests to check."

  if prev not in {"feed", "check"}:
    print(divider)
  divider = "{0:-<{1}}+----+".format("+", OUT_LEN - 6)
  result = "|{{0: <{0}}}|{{1: <4}}|".format(OUT_LEN - 7)
  error = "|{{0: <{0}}}|".format(OUT_LEN - 2)
  if prev != "check":
    print(divider)
    print(result.format("NAME", "RES"))
    print(divider)

  for test in tests:
    for setting in settings(params):
      name = run(test, "{0} -p".format(flags(setting)))[1].replace("\n", "")
      assert len(name) < OUT_LEN - 7, "Error: check name too long."
      (rc, out) = run(test, flags(setting))
      if rc == 0:
        print(result.format(name, "PASS"))
        passed += 1
      elif rc == 125:
        print(result.format(name, "N/A"))
        not_run += 1
      else:
        print(result.format(name, "FAIL"))
        failed += 1
      if rc != 0:
        for line in out.split("\n"):
          while line:
            print(error.format(line[:OUT_LEN - 2]))
            line = line[OUT_LEN - 2:]
  prev = "check"

def feed():
  global divider
  global prev
  if prev != "feed":
    print(divider)
  print("")
  divider = "+" + (OUT_LEN - 2) * "-" + "+"
  prev = "feed"

def title(desc):
  feed()
  section(desc)

def section(desc):
  global divider
  global prev
  assert len(desc) < OUT_LEN-2, "Error: desc too long."
  if prev != "section":
    print(divider)
  print("|{0: <{1}}|".format(desc, OUT_LEN - 2))
  prev = "section"

def suite(name):
  assert os.path.isfile(name), "Error: suite does not exist."
  assert os.path.splitext(name)[1] == ".suite", "Error: incorrect suite file type."
  f = open(name, 'r')
  parentcwd = os.getcwd()
  os.chdir(os.path.dirname(os.path.realpath(name)))
  for line in f.readlines():
    line = line.strip()
    exec(line)
  os.chdir(parentcwd)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Test Harness")
  parser.add_argument('-s', '--suites', default="check.suite", nargs='+', type=str, help='test suites to run')

  args = parser.parse_args()
  for name in args.suites:
    suite(name)
  for test in built:
    clean(test)
  feed()
  if (passed + failed + not_run) > 0:
    print(divider)
    emoticon = ":D"
    if (failed + not_run) > 0:
      emoticon = ":("
    print("|{0: <{1}} {2}  |".format("TESTS PASSED: {0}   TESTS FAILED: {1}   TESTS NOT RUN: {2}".format(passed, failed, not_run), OUT_LEN - 7, emoticon))
    print(divider)
  if not_run > 0:
    exit(125)
  elif failed > 0:
    exit(1)
  else:
    exit(0)
