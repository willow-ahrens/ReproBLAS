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

def get_tuning(file_name):
  tuning_name = os.path.splitext(os.path.abspath(file_name))[0] + ".tune"
  assert os.path.isfile(tuning_name), "Error: tuning file does not exist."
  f = open(tuning_name, 'r')
  try:
    return eval(f.readlines()[0])
  except (ValueError, SyntaxError):
    assert False, "Error: empty or corrupt tuning file."

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

def generate(target, file_name, tuning):
  output_name = os.path.splitext(os.path.abspath(file_name))[0] + ".c"
  output_file = SrcFile(output_name)
  output_block = output_file.sub_block()

  output_block.new_line()

  for (i, (vec_name, vec_settings)) in enumerate(tuning):
    vectorization = vectorization_lookup[vec_name]
    if len(tuning) > 1:
      if i == 0:
        output_block.write("#if defined( " + vectorization.defined_macro + " )")
      elif i < len(tuning) - 1:
        output_block.write("#elif defined( " + vectorization.defined_macro + " )")
      else:
        output_block.write("#else")
      output_block.indent()
    target.write(output_block, vectorization, vec_settings)
    if len(tuning) > 1:
      output_block.dedent()
  if len(tuning) > 1:
    output_block.write("#endif")

  output_file.dump()
