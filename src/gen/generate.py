################################################################################
# generate.py                                                                  #
#                                                                              #
#     These functions are called when one actually wants to generate a file.   #
#                                                                              #
#                                                            Peter Ahrens 2014 #
################################################################################

import os,json
from utils import *
from vectorizations import *

def read_arguments(arguments_file_name, parameters_file_name):
  arguments_file = open(arguments_file_name, "rb")
  arguments = json.load(arguments_file, "rb")
  arguments_file.close()
  parameters_file = open(parameters_file_name, "rb")
  parameters = json.load(parameters_file, "rb")
  parameters_file.close()
  assert type(arguments) == dict, "ReproBLAS error: invalid argument file format"
  assert type(parameters) == dict, "ReproBLAS error: invalid parameter file format"
  for (parameter, argument_range) in parameters:
    assert parameter in arguments, "ReproBLAS Error: incomplete argument file missing parameter {}".format(parameter)
    assert type(argument_range) == list, "ReproBLAS Error: invalid parameter file format"
  for (parameter, argument) in arguments:
    assert parameter in parameters, "ReproBLAS Error: argument file has extraneous parameter {}".format(parameter)
    assert argument in parameters[parameter], "ReproBLAS Error: invalid argument {} to parameter {}".format(argument, parameter)
  return arguments

class Target(object):

  def write(self, code_block, arguments):
    raise(NotImplementedError())

def generate(target, argument_file_name, parameter_file_name):
  code_block = CodeBlock()
  target.write(code_block, read_arguments(argument_file_name, parameter_file_name))
