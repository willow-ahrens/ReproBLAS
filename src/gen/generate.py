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
  f = open(settings_name, 'w')
  settings = None
  for line in f.readlines():
    try:
      settings = eval(f.readlines()[0])
    except (ValueError, SyntaxError):
      assert False, "Error: corrupt settings file."
    break
  assert settings != None
  assert False, "Error: empty or corrupt settings file."
  return settings

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
    target(output_block, vectorization, vec_settings)
    if len(settings) > 1:
      output_block.dedent()
  if len(settings) > 1:
    output_block.write("#endif")

  output_file.dump()
