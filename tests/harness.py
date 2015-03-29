#Harness.py
#the intended test harness for all tests here.
#By peter ahrens
#comments/code niceness to come


import argparse
import os, subprocess, sys
import re
import decimal

OUT_LEN = 80

divider = "+" + (OUT_LEN - 2) * "-" + "+"
prev    = "feed"
built   = {}
passed  = 0
failed  = 0
not_run = 0
style   = None

class Style:
  @classmethod
  def create_template(cls, left_align, right_align):
    template = ""
    left_expanded_columns = ["{{:{0}<{1}}}".format(fill, size) for (fill, size) in left_align]
    right_expanded_columns = ["{{:{0}<{1}}}".format(fill, size) for (fill, size) in right_align]
    template = "{0}{1}{0}".format(cls.edge, cls.inner.join(left_expanded_columns + right_expanded_columns))
    example_output = template.format(*["" for _ in left_align + right_align])
    if cls.min_len != 0 and len(example_output) < cls.min_len:
      if(left_expanded_columns):
        left_expanded_columns = left_expanded_columns + [""]
      if(right_expanded_columns):
        right_expanded_columns = [""] + right_expanded_columns
      left_template = "{0}{1}".format(cls.edge, cls.inner.join(left_expanded_columns))
      right_template = "{0}{1}".format(cls.inner.join(right_expanded_columns), cls.edge)
      incomplete_example_output = (left_template + right_template).format(*["" for _ in left_align + right_align])
      template = "{0}{1}{2}".format(left_template, cls.fill * (cls.min_len - len(incomplete_example_output)), right_template)
    return template

  @classmethod
  def create_divider(cls, left_align, right_align):
    if cls.divider.inner == "" and cls.divider.edge == "" and cls.divider.fill == "":
      return ""
    return cls.divider.create_template(left_align, right_align).format(*[cls.divider.fill * size for (fill, size) in left_align + right_align])

class CSVDivider(Style):
  inner   = ""
  edge    = ""
  fill    = ""
  min_len = 0

class CSV(CSVDivider):
  inner = ", "
  edge  = ""
  fill  = ""
  divider = CSVDivider

class TermDivider(Style):
  inner   = "+"
  edge    = "+"
  fill    = "-"
  min_len = OUT_LEN

class Term(TermDivider):
  inner = "|"
  edge  = "|"
  fill  = " "
  divider = TermDivider

styles = {"term":Term, "csv":CSV}

def callsafe(command):
  rc = 0
  try:
    out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True).decode(sys.stdout.encoding)
  except subprocess.CalledProcessError as e:
    rc = e.returncode
    out = e.output.decode(sys.stdout.encoding)
  return (rc, out)

def call(command):
  return subprocess.check_output(command, shell=True).decode(sys.stdout.encoding)

def run(test, args):
  global built
  test = os.path.realpath(test)
  if test not in built:
    built[test] = os.path.join(call("cd {0}; make pbd".format(os.path.split(test)[0])).split()[-1], os.path.split(test)[1])
    call("make {0}".format(built[test]))
    assert os.path.isfile(built[test]), "Error: make unsuccessful."
  return callsafe("{0} {1}".format(built[test], args))

def settings(params):
  if params:
    for value in params[0][1]:
      for setting in settings(params[1:]):
        yield [(params[0][0], value)] + setting
  else:
    yield []

def flags(setting):
  return " ".join(['-{0} "{1}"'.format(*s) if len(s[0]) == 1 else '--{0} "{1}"'.format(*s) for s in setting])

def engineer(f, d):
  f = decimal.Decimal(f)
  f = f.normalize().to_eng_string()
  if (f.find("E+") == -1):
    m = f
    e = ""
    #e = "e0"
  else:
    (m, e) = f.split("E+")
    e = "e" + e
  m = "{{0:0<{0}.{0}}}".format(d - len(e)).format(m)
  return m + e

def benchmark(tests, params):
  global divider
  global prev
  BENCH_LEN = 9

  if prev != "feed":
    print(divider)

  result_template = style.create_template([(" ", len(str(max(param[1])))) for param in params], [(" ", BENCH_LEN) for test in tests])
  divider = style.create_divider([(" ", len(str(max(param[1])))) for param in params], [(" ", BENCH_LEN) for test in tests])

  labels = list(list(zip(*params))[0])
  for test in tests:
    if test != "":
      try:
        name = re.findall("\[(.*)\]", run(test, "{0} -p")[1])[0]
      except IndexError:
        assert False, "Error: benchmark name does not match format."
      assert len(name) < BENCH_LEN, "Error: test name too long."
      labels.append(name)
    else:
      labels.append("")
  print(divider)
  print(result_template.format(*labels))
  print(divider)
  for setting in settings(params):
    results = list(list(zip(*setting))[1])
    for test in tests:
      if test != "":
        try:
          results.append(engineer(re.findall("([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)", run(test, flags(setting))[1])[0][0], BENCH_LEN))
        except (ValueError, IndexError):
          assert False, "Error: benchmark output does not match format"
      else:
        results.append("")
    print(result_template.format(*results))
  prev = "benchmark"

def script(tests, params):
  global divider
  global prev

  result_template = style.create_template([(" ", OUT_LEN - 2)], [])
  if prev != "script":
    print(divider)
  for test in tests:
    for setting in settings(params):
      name = "{0} {1}".format(test, flags(setting))
      print(result_template.format(name))
      (rc, out) = run(test, flags(setting))
      if rc != 0:
        for line in out.split("\n"):
          while line:
            print(result_template.format(line[:OUT_LEN - 2]))
            line = line[OUT_LEN - 2:]
  prev = "script"

def check(tests, params):
  global divider
  global prev
  global passed
  global failed
  global not_run
  global style

  assert tests, "Error: no tests to check."

  if prev not in {"feed", "check"}:
    print(divider)
  
  result_template = style.create_template([(" ", OUT_LEN - 7)], [(" ", 4)])
  divider = style.create_divider([(" ", OUT_LEN - 7)], [(" ", 4)])
  error_template = style.create_template([(" ", OUT_LEN - 2)], [])
  if prev != "check":
    print(divider)
    print(result_template.format("NAME", "RES"))
    print(divider)

  for test in tests:
    for setting in settings(params):
      name = run(test, "{0} -p".format(flags(setting)))[1].replace("\n", "")
      #assert len(name) < OUT_LEN - 7, "Error: check name too long."
      (rc, out) = run(test, flags(setting))
      if rc == 0:
        print(result_template.format(name, "PASS"))
        passed += 1
      elif rc == 125:
        print(result_template.format(name, "N/A"))
        not_run += 1
      else:
        print(result_template.format(name, "FAIL"))
        print(error_template.format("error {}:".format(rc)))
        failed += 1
      if rc != 0:
        for line in out.split("\n"):
          while line:
            print(error_template.format(line[:OUT_LEN - 2]))
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
  parser.add_argument('-o', '--format', default="term", choices=styles.keys(), help='output format')

  args = parser.parse_args()
  style = styles[args.format]
  for name in args.suites:
    suite(name)
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
