################################################################################
# generate.py                                                                  #
#                                                                              #
#     These functions are called when one actually wants to generate a file.   #
#                                                                              #
#                                                            Peter Ahrens 2014 #
################################################################################

import os,json
import cog
from utils import *
from vectorizations import *

def deserialize_arguments(arguments_file_name):
  arguments_file = open(arguments_file_name, "r")
  arguments = json.load(arguments_file)
  arguments_file.close()
  assert type(arguments) == dict, "ReproBLAS error: invalid argument file format"
  for key in arguments.keys():
    assert type(key) == str, "ReproBLAS error: invalid argument file format"
  return arguments

def serialize_arguments(arguments, arguments_file_name):
  assert type(arguments) == dict, "ReproBLAS error: invalid argument file format"
  for key in arguments.keys():
    assert type(key) == str, "ReproBLAS error: invalid argument file format"
  arguments_file = open(arguments_file_name, "w")
  json.dump(arguments_file, arguments)
  arguments_file.close()

def serialize_parameter_space(parameter_space, parameter_space_file_name):
  parameter_space_file = open(parameter_space_file_name, "w")
  json.dump(parameter_space_file, parameter_space)
  parameter_space_file.close()

def deserialize_parameter_space(parameter_space, parameter_space_file_name):
  parameter_space_file = open(parameter_space_file_name, "r")
  parameter_space = ParameterSpace.decode(json.dump(parameter_space_file, parameter_space))
  parameter_space_file.close()

class Parameter:
  def __init__(self, name):
    self.name = name

  def parse_value(self, value):
    raise NotImplementedError

  def encode(self):
    return {"name":self.name, "flavor":self.flavor}

  @classmethod
  def decode(data):
    assert type(data) == dict, "ReproBLAS error: invalid parameter file format"
    assert "name" in data, "ReproBLAS error: invalid parameter file format"
    assert "flavor" in data, "ReproBLAS error: invalid parameter file format"
    assert "default" in data, "ReproBLAS error: invalid parameter file format"
    if data["flavor"] == "integer":
      assert "minimum" in data, "ReproBLAS error: invalid parameter file format"
      assert "maximum" in data, "ReproBLAS error: invalid parameter file format"
      assert "step" in data, "ReproBLAS error: invalid parameter file format"
      return IntegerParameter(data["name"], data["minimum"], data["maximum"], data["step"], data["default"])

class IntegerParameter(Parameter):
  def __init__(self, name, minimum, maximum, step, default):
    super(IntegerParemeter, self).__init__(self, name)
    self.flavor = "integer"
    self.minimum = minimum
    self.maximum = maximum
    self.step = step
    assert minimum % step == 0, "ReproBLAS error: integer parameter minimum must be multiple of step"
    assert maximum % step == 0, "ReproBLAS error: integer parameter maximum must be multiple of step"
    self.default = self.parse_value(default)

  def parse_value(self, value):
    value = int(value)
    assert value >= minimum, "ReproBLAS error: integer parameter must be >= min"
    assert value < maximum, "ReproBLAS error: integer parameter must be < max"
    assert value % step == 0, "ReproBLAS error: integer parameter value must be multiple of step"

  def encode(self):
    data = super(IntegerParameter, self).encode
    data["minimum"] = self.minimum
    data["maximum"] = self.maximum
    data["step"] = self.step
    data["default"] = self.default

class ParameterSpace:
  def __init__(self):
    self.forward_dependencies = {} #file_name > arguments
    self.backward_dependencies = {} #argument > file_names
    self.forward_metrics = {} #metric > arguments
    self.backward_metrics = {} #argument > metrics
    self.parameters = {}

  def encode(self):
    return {"forward_dependencies":{file_name:list(arguments) for (file_name, arguments) in self.forward_dependencies.items()},\
            "backward_dependencies":{argument:list(file_names) for (argument, file_names) in self.backward_dependencies.items()},\
            "forward_metrics":{metric:list(arguments) for (metric, arguments) in self.forward_metrics.items()},\
            "backward_metrics":{argument:list(metrics) for (argument, metrics) in self.backward_metrics.items()},\
            "parameters":{parameter_name:parameter.encode() for (parameter_name, parameter) in self.parameters.items()]}

  @classmethod
  def decode(data):
    parameter_space = ParameterSpace()
    assert type(data) == dict, "ReproBLAS error: invalid parameter file format"
    assert "forward_dependencies" in data, "ReproBLAS error: invalid parameter file format"
    assert "backward_dependencies" in data, "ReproBLAS error: invalid parameter file format"
    assert "forward_metrics" in data, "ReproBLAS error: invalid parameter file format"
    assert "backward_metrics" in data, "ReproBLAS error: invalid parameter file format"
    assert "parameters" in data, "ReproBLAS error: invalid parameter file format"
    parameter_space.forward_dependencies = {file_name:set(arguments) for (file_name, arguments) in data["forward_dependencies"].items()}
    parameter_space.backward_dependencies = {argument:set(file_names) for (argument, file_names) in data["backward_dependencies"].items()}
    parameter_space.forward_metrics = {metric:set(arguments) for (metric, arguments) in data["forward_metrics"].items()}
    parameter_space.backward_metrics = {argument:set(metrics) for (argument, metrics) in data["backward_metrics"].items()}
    parameter_space.parameters = {parameter_name:Parameter.decode(parameter) for (parameter_name, parameter) in data["parameters"].items()]

  def add_target(self, target):
    if target.file_name not in self.forward_dependencies:
      self.forward_dependencies[target.file_name] = set()
    for argument in target.get_arguments():
      self.forward_dependencies[target.file_name].add(argument)

    for argument in target.get_arguments():
      if argument not in self.backward_dependencies:
        self.backward_dependencies[argument] = set()
      self.backward_dependencies[argument].add(target.file_name)

    for argument in target.get_arguments():
      if argument not in self.backward_metrics:
        self.backward_metrics[arguments] = set()
      for metric in target.get_metrics():
        self.backward_metrics[argument].add(metric)

    for metric in target.get_metrics():
      if metric not in self.forward_metrics:
        self.forward_metrics[metric] = set()
      for argument in target.get_arguments():
        self.forward_metrics[metric].add(argument)

    for parameter in target.get_parameters():
      if parameter.name in self.parameters:
        assert False, "ReproBLAS error: duplicate parameter"
      self.parameters[parameter.name] = parameter

  def get_value(self, argument, arguments):
    assert argument in arguments, "ReproBLAS error: missing argument data"
    assert argument in self.parameters, "ReproBLAS error: missing parameter data"
    return self.parameters[argument].parse_value(arguments[argument])

class Target(object):
    """
    A Target is a target for code generation. Override "get_parameters",
    "get_arguments", "get_metrics", and "write" to use.
    """

  def __init__(self):
    self.file_name = cog.inFile

  def get_parameters(self):
    """
    Return a list of Parameter objects that this target is responsible for.
    """
    raise(NotImplementedError())

  def get_metrics(self):
    """
    Return a list of metric names target affects.
    """
    raise(NotImplementedError())

  def get_arguments(self):
    """
    Return a list of argument names target will use during it's execution.
    """
    raise(NotImplementedError())

  def set_arguments(self, arguments, parameter_space):
    """
    Given the current argument list and ParameterSpace object, extract the necessary
    information.
    """
    self.arguments = {}
    for argument in self.get_arguments():
      self.arguments[name] = parameter_space.get_value(argument, arguments)

  def write(self, code_block):
    """
    Write the generated code into the given code_block. At this point, the
    arguments field should have been initialized.
    """
    raise(NotImplementedError())

def generate(target, args, params, mode):
  """
  Given the target, run the generator and return the output. This function
  should be called on a target <target> as follows:
  cog.outl(generate.generate(<target>, args, params, mode))
  """
  code_block = CodeBlock()

  if mode == "generate":
    target.set_arguments(deserialize_arguments(args), deserialize_parameter_space(params))
    target.write(code_block)
  if mode == "params":
    if os.path.isfile(params):
      parameter_space = deserialize_parameter_space(params)
    else:
      parameter_space = ParameterSpace()
    parameter_space.add_target(target)
    serialize_parameter_space(parameter_space, params)
  return str(code_block)
