#Harness.py
#the intended test harness for all tests here.
#By peter ahrens
#comments/code niceness to come

import argparse
import copy
import itertools
import os
import re
import json

import config
import scripts.terminal as terminal
import scripts.texttable.texttable as texttable

import timeit

class Harness(object):
  def __init__(self, name):
    self.name = name
    parser = argparse.ArgumentParser(description = name)
    parser.add_argument('-f', '--format', default="term", choices=["term", "csv"], help='output format')
    parser.add_argument('-r', '--runmode', default="sequential", choices=["sequential", "parallel"], help='run mode')
    parser.add_argument('-v', '--verbose', default="false", type=str, nargs="?", help='verbose if "true"')
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
      suite.setup(verbose = self.args.verbose, **kwargs)
      command_list += suite.get_command_list()
    if self.args.runmode == "sequential":
      output_list = config.run(command_list, verbose = self.args.verbose)
    else:
      output_list = config.run_parallel(command_list, verbose=self.args.verbose)
    i = 0
    for suite in self.suites:
      suite.parse_output_list(output_list[i:i + len(suite.get_command_list())])
      i += len(suite.get_command_list())
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

  def __init__(self, metrics, params, ranges, attribute, silent_flags=""):
    self.params = params
    self.ranges = ranges
    self.attribute = attribute
    self.metric_rows = []
    self.argss = []
    self.metrics = []
    self.silent_flags = silent_flags
    for args in itertools.product(*ranges):
      self.argss.append(args)
      row = [copy.deepcopy(metric) for metric in metrics]
      for metric in row:
        self.metrics.append(metric)
      self.metric_rows.append(row)

  def setup(self, **kwargs):
    for (metric_row, args) in zip(self.metric_rows, self.argss):
      for metric in metric_row:
        metric.setup(attribute = self.attribute, flagss = [terminal.flags(self.params, args)], silent_flags = self.silent_flags, **kwargs)

  def get_command_list(self):
    command_list = []
    for metric in self.metrics:
      command_list += metric.get_command_list()
    return command_list

  def parse_output_list(self, output_list):
    i = 0
    for metric in self.metrics:
      metric.parse_output_list(output_list[i:i + metric.get_num_commands()])
      i += metric.get_num_commands()

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

  def get_num_commands(self):
    """
    return the number of commands in the command_list
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

  def setup(self, flagss=[""], silent_flags="", **kwargs):
    self.silent_flags = silent_flags
    self.flagss = flagss
    self.executable_output = terminal.make(self.executable, **kwargs)

  def get_command_list(self):
    return ["{} {} {} {}".format(self.executable_output, self.base_flags, self.silent_flags, flags) for flags in self.flagss]

  def get_num_commands(self):
    return 1

class MetricTest(ExecutableTest):

  def get_name(self):
    return self.name

  def setup(self, attribute="", **kwargs):
    self.attribute = attribute
    super(MetricTest, self).setup(**kwargs)

  def parse_output_list(self, output_list):
    assert len(output_list) == len(self.flagss), "ReproBLAS error: unexpected test output"

    self.result = 0
    self.output = []
    for output_item in output_list:
      self.output.append(json.loads(output_item[1]))
      self.result += self.parse_output(self.output[-1])

  def parse_output(self, output):
    raise(NotImplementedError())

  def get_output(self):
    return self.output

  def get_result(self):
    return self.result
