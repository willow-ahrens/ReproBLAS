#Harness.py
#the intended test harness for all tests here.
#By peter ahrens
#comments/code niceness to come


import argparse
import os, subprocess, sys
import re
import decimal
import texttable.texttable as texttable
import itertools
import copy

class Harness(object)
  def __init__(self, name):
    self.name = name
    parser = argparse.ArgumentParser(description = name)
    parser.add_argument('-f', '--format', default="term", choices=["term", "csv"], help='output format')
    args = parser.parse_args()
    if args.format == "term":
      self.table = textable.Texttable(max_width = 80)
      self.table.set_precision(9)
    if args.format == "csv":
      self.table = texttable.Texttable(max_width = 0)
      self.table.set_chars(["", ", ", "", ""]
      self.table.set_deco(texttable.Texttable.VLINES)
    self.suites = []

  def add_suite(self, suite):
    self.suites.append(suite)

  def add_suites(self, suites):
    self.suites += suites

  def run(self):
    command_list = []
    for suite in self.suites:
      suite.setup()
      command_list += suite.get_command_list()
    output_list = terminal.run(command_list)
    for suite in self.suites:
      suite.parse_output_list(output_list[:len(suite.get_command_list())])
      output_list = output_list[len(suite.get_command_list()):]
    for suite in self.suites:
      self.table.set_header(suite.get_header())
      self.table.set_cols_align(suite.get_align())
      self.table.set_cols_dtype(suite.get_dtype())
      if suite.get_cols_width():
        if args.format == "term":
          self.table.set_cols_width(suite.get_cols_width(80))
        if args.format == "csv":
      self.table.add_rows(suite.get_rows())
      self.table.draw()
      self.table.reset()

class Suite(object):
  def setup(self):
    """
    call necessary commands to setup suite (i.e. build the executables)
    """
    raise(NotImplementedError())
  def get_command_list(self):
    """
    return a list of commands that constitute the suite to be run on the
    target architecture
    """
    raise(NotImplementedError())
  def parse_output_list(self, output_list):
    """
    parse the output of the command set. The output will be given as a list of
    (return code, output)
    """
    raise(NotImplementedError())
  def get_header(self):
    """
    return the header names for our output as a list.
    see texttable.py for details.
    """
    raise(NotImplementedError())
  def get_align(self):
    """
    return the alignment for our output fields as a list of "l", "c", or "r".
    see texttable.py for details.
    """
    raise(NotImplementedError())
  def get_dtype(self):
    """
    return the datatype for our output fields as a list of "a", "t", "f", "e",
    or "r".
    see texttable.py for details.
    """
    raise(NotImplementedError())
  def get_cols_width(self, max_width):
    """
    return either a list of integer column widths or None if you want autosize.
    """
    raise(NotImplementedError())
  def get_rows(self):
    """
    return a list of lists of each datatype in the rows. parse_output_list will
    have been called when this is called.
    """
    raise(NotImplementedError())

class CheckSuite(Suite):

  def __init__(self):
    self.checks = []
    self.check_rows = []
    self.args = []
    self.params = []

  def add_checks(self, checks, params, ranges):
    for args in itertools.product(*ranges):
      check_row = [copy.deepcopy(checks) for check in checks]
      self.check_rows.append(check_row)
      self.checks += check_row
      self.params.append(params)
      self.args.append(args)

  def setup(self, **kwargs):
    for check_row, params, args in zip(self.check_rows, self.params, self.args):
      for check in check_row:
        check.setup(flags = terminal.flags(params, args), **kwargs)

  def get_command_list(self):
    command_list = []
    for check in self.checks:
      command_list += check.get_command_list():
    return command_list

  def parse_output_list(self, output_list):
    for check in self.checks:
      check.parse_output_list(output_list[:len(check.get_command_list)])
      output_list = output_list[len(check.get_command_list):]

  def get_header(self):
    return ["Check", "Res"]

  def get_align(self):
    return ["l", "c"]

  def get_dtype(self):
    return ["t", "t"]

  def set_cols_width(self, max_width):
    return [max_width - 1 - 1 - 4 - 1, 4]

  def get_rows(self):
    passed = 0
    failed = 0
    na = 0
    rows = []
    for check in self.checks:
      if check.get_result() == 0:
        rows.append([check.get_name(), "Pass"])
        passed += 1
      elif check.get_result() == 125:
        rows.append([check.get_output(), "N/A"])
        na += 1
      else:
        rows.append([check.get_output(), "Fail"])
        failed += 1
    emoticon = ":("
    if passed == len(self.checks):
      emoticon = ":D"
    rows.append(["Passed: {} Failed: {} N/A: {}".format(passed, failed, na), emoticon])

class MetricSuite(Suite):

  def __init__(self, metrics, params, ranges):
    self.params = params
    self.ranges = ranges
    self.metric_rows = []
    self.args = []
    self.metrics = []
    for args in itertools.product(*ranges):
      self.args.append(args)
      row = [copy.deepcopy(metric) for metric in metrics]
      for metric in row:
        self.metrics += metric
      self.metric_rows.append(row)

  def setup(self, **kwargs):
    for (metric_row, args) in zip(self.metric_rows, self.args):
      for metric in metric_row:
        metric.setup(flags = terminal.flags(params, args),**kwargs)

  def get_command_list(self):
    command_list = []
    for metric in self.metrics:
      command_list += metric.get_command_list():
    return command_list

  def parse_output_list(self, output_list):
    for metric in self.metrics:
      metric.parse_output_list(output_list[:len(metric.get_command_list())])
      output_list = output_list[len(metric.get_command_list()):]

  def get_header(self):
    return self.params + [metric.get_name() for metric in self.metric_rows[0]]

  def get_align(self):
    return ["l" for _ in self.get_header()]

  def get_dtype(self):
    return ["a" for _ in self.params] + ["f" for _ in self.metric_rows[0]]

  def get_cols_width(self, max_width):
    return None

  def get_rows(self):
    rows = []
    for metric_row in self.metric_rows:
      row = copy.deepcopy(self.params)
      for metric in metric_row:
        row.append(metric.result())
    return rows

class Test(object):

  def get_name(self):
    """
    return the name of the test
    """
    raise NotImplementedError()

  def setup(self, **kwargs):
    """
    call necessary commands to setup test (i.e. build the executables)
    """
    raise(NotImplementedError())

  def get_command_list(self):
    """
    return a list of commands that constitute the test to be run on the
    target architecture
    """
    raise(NotImplementedError())

  def parse_output_list(self, output_list):
    """
    parse the output of the command set. The output will be given as a list of
    (return code, output)
    """
    raise(NotImplementedError())

  def get_output(self):
    """
    return all relevant output (mostly for debugging)
    """
    raise(NotImplementedError())

  def get_output(self):
    """
    return test result
    """
    raise(NotImplementedError())

class ExecutableTest(Test):

  def setup(self, **kwargs):
    make(self.executable, **kwargs)

  def (self):
    rai

      try:
        name = re.findall("\[(.*)\]", run(test, "{0} -p")[1])[0]
      except IndexError:
        assert False, "Error: benchmark name does not match format."
      assert len(name) < BENCH_LEN, "Error: test name too long."
    
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
