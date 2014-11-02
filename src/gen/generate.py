################################################################################
# generate.py                                                                  #
#                                                                              #
#     These functions are called when one actually wants to generate a file.   #
#                                                                              #
#                                                            Peter Ahrens 2014 #
################################################################################

import os
from utils import *
from vectorizations import *

def get_settings(file_name):
  settings_name = os.path.splitext(os.path.abspath(file_name))[0] + ".set"
  assert os.path.isfile(settings_name), "Error: settings file does not exist."
  f = open(settings_name, 'r')
  settings = None
  try:
    return eval(f.readlines()[0])
  except (ValueError, SyntaxError):
    assert False, "Error: empty or corrupt settings file."

class Target(object):

  def write(self, code_block, vectorization, settings):
    raise(NotImplementedError())

class Function(Target):
  name = ""

  def __init__(self, data_type_class):
    self.data_type_class = data_type_class

  def write(self, code_block, vectorization, settings):
    self.data_type = self.data_type_class(code_block)
    self.write_declaration(code_block, settings)
    code_block.indent()
    self.data_type = self.data_type_class(code_block)
    body_block = code_block.sub_block()
    self.vec = vectorization(body_block, self.data_type_class)
    self.write_body(body_block, settings)
    code_block.dedent()
    code_block.write("}")

  def write_declaration(self, code_block, settings):
    raise(NotImplementedError())

  def write_body(self, code_block, settings):
    raise(NotImplementedError())

def generate(target, file_name, settings):
  output_name = os.path.splitext(os.path.abspath(file_name))[0] + ".c"
  output_file = SrcFile(output_name)
  output_block = output_file.sub_block()

  output_block.new_line()

  for (i, (vec_name, vec_settings)) in enumerate(settings):
    vectorization = vectorization_lookup[vec_name]
    if len(settings) > 1:
      if i == 0:
        output_block.write("#if defined( " + vectorization.defined_macro + " )")
      elif i < len(settings) - 1:
        output_block.write("#elif defined( " + vectorization.defined_macro + " )")
      else:
        output_block.write("#else")
      output_block.indent()
    target.write(output_block, vectorization, vec_settings)
    if len(settings) > 1:
      output_block.dedent()
  if len(settings) > 1:
    output_block.write("#endif")

  output_file.dump()
