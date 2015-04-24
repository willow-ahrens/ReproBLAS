#Harness.py
#the intended test harness for all tests here.
#By peter ahrens
#comments/code niceness to come


import argparse
import copy
import itertools
import os
import re

import config
import scripts.terminal as terminal
import scripts.texttable.texttable as texttable

class Harness(object):
  def __init__(self, name):
    self.name = name
    parser = argparse.ArgumentParser(description = name)
    parser.add_argument('-f', '--format', default="term", choices=["term", "csv"], help='output format')
    self.args = parser.parse_args()
    if self.args.format == "term":
      self.table = texttable.Texttable(max_width = 80)
      self.table.set_precision(4)
    if self.args.format == "csv":
      self.table = texttable.Texttable(max_width = 0)
      self.table.set_chars(["", ", ", "", ""])
      self.table.set_deco(texttable.Texttable.VLINES)
    self.suites = []

  def add_suite(self, suite):
    self.suites.append(suite)

  def add_suites(self, suites):
    self.suites += suites

  def run(self, **kwargs):
    command_list = []
    for suite in self.suites:
      suite.setup(**kwargs)
      command_list += suite.get_command_list()
    output_list = config.run(command_list)
    for suite in self.suites:
      suite.parse_output_list(output_list[:len(suite.get_command_list())])
      output_list = output_list[len(suite.get_command_list()):]
    for suite in self.suites:
      tablecopy = copy.deepcopy(self.table)
      tablecopy.set_cols_align(suite.get_align())
      tablecopy.set_cols_dtype(suite.get_dtype())
      if self.args.format == "term":
        if suite.get_cols_width(80):
          tablecopy.set_cols_width(suite.get_cols_width(80))
      tablecopy.add_rows([suite.get_header()] + suite.get_rows())
      print(tablecopy.draw())

class Suite(object):

  def setup(self):
    """
    call necessary commands to setup suite (i.e. build the executables)
    """
    raise NotImplementedError()

  def get_command_list(self):
    """
    return a list of commands that constitute the suite to be run on the
    target architecture
    """
    raise NotImplementedError()

  def parse_output_list(self, output_list):
    """
    parse the output of the command set. The output will be given as a list of
    (return code, output)
    """
    raise NotImplementedError()

  def get_header(self):
    """
    return the header names for our output as a list.
    see texttable.py for details.
    """
    raise NotImplementedError()

  def get_align(self):
    """
    return the alignment for our output fields as a list of "l", "c", or "r".
    see texttable.py for details.
    """
    raise NotImplementedError()

  def get_dtype(self):
    """
    return the datatype for our output fields as a list of "a", "t", "f", "e",
    or "r".
    see texttable.py for details.
    """
    raise NotImplementedError()

  def get_cols_width(self, max_width):
    """
    return the width of each column as a list of integers or None for auto
    see texttable.py for details.
    """
    raise NotImplementedError()

  def get_rows(self):
    """
    return a list of lists of each datatype in the rows. parse_output_list will
    have been called when this is called.
    """
    raise(NotImplementedError())

class MetricSuite(Suite):

  def __init__(self, metrics, params, ranges):
    self.params = params
    self.ranges = ranges
    self.metric_rows = []
    self.argss = []
    self.metrics = []
    for args in itertools.product(*ranges):
      self.argss.append(args)
      row = [copy.deepcopy(metric) for metric in metrics]
      for metric in row:
        self.metrics.append(metric)
      self.metric_rows.append(row)

  def setup(self, **kwargs):
    for (metric_row, args) in zip(self.metric_rows, self.argss):
      for metric in metric_row:
        metric.setup(flags = terminal.flags(self.params, args),**kwargs)

  def get_command_list(self):
    command_list = []
    for metric in self.metrics:
      command_list += metric.get_command_list()
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
    return ["a" for _ in self.params] + ["e" for _ in self.metric_rows[0]]

  def get_cols_width(self, max_width):
    return None

  def get_rows(self):
    rows = []
    for (metric_row, args) in zip(self.metric_rows, self.argss):
      row = list(args)
      for metric in metric_row:
        row.append(metric.get_result())
      rows.append(row)
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

  def get_result(self):
    """
    return test result
    """
    raise(NotImplementedError())

class ExecutableTest(Test):
  base_flags = ""

  def setup(self, flags="", **kwargs):
    self.flags = flags
    self.executable_output = terminal.make(self.executable, **kwargs)

  def get_command_list(self):
    """
    return a list of commands that constitute the test to be run on the
    target architecture
    """
    return ["{} {} {}".format(self.executable_output, self.base_flags, self.flags)]

class MetricTest(ExecutableTest):

  def get_name(self):
    """
    return the name of the test
    """
    return self.name

  def parse_output_list(self, output_list):
    """
    parse the output of the command set. The output will be given as a list of
    (return code, output)
    """
    assert len(output_list) == 1, "ReproBLAS error: unexpected test output"
    self.output = output_list[0][1]
    self.result = float(re.findall("([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)", output_list[0][1])[0][0])

  def get_output(self):
    """
    return all relevant output (mostly for debugging)
    """
    return self.output

  def get_result(self):
    """
    return test result
    """
    return self.result

