import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
from utils import *
from dataTypes import *
from vectorizations import *
from generate import *
import itertools

class AMaxP(Function):
  standard_incs = ["incv"]

  def __init__(self, data_type_class):
    super(AMaxP, self).__init__(data_type_class)

  def write_declaration(self, code_block, settings):
    code_block.srcFile.include("#include <stdio.h>")
    code_block.srcFile.include("#include <stdlib.h>")
    code_block.srcFile.include("#include <float.h>")
    code_block.srcFile.include("#include <math.h>")
    code_block.srcFile.include("#include \"config.h\"")
    code_block.srcFile.include("#include \"Common/Common.h\"")
    #code_block.srcFile.include("#include \"IndexedFP/" + self.data_type.base_type.name_char + "Indexed.h\"")
    #code_block.srcFile.include("#include \"rblas1.h\"")

  def write_body(self, code_block, settings = [1]):
      max_unroll = settings[0]
      max_process_width = self.compute_process_width(max_unroll)
      code_block.write("int i;")
      code_block.define_vars(self.data_type.name, ["max"])
      code_block.new_line()
      self.define_load_ptrs(code_block, max_process_width)
      self.define_load_vars(code_block, max_process_width)
      self.m_vars = ["m_" + str(i) for i in range(self.vec.suf_width)]
      code_block.define_vars(self.vec.type_name, self.m_vars)
      code_block.set_equal(self.m_vars, itertools.repeat(self.vec.zero));

      code_block.new_line()

      code_block.write("if(" + " && ".join([inc + " == 1" for inc in self.standard_incs]) + "){")
      code_block.indent()
      self.write_core(code_block, max_process_width, max_unroll, [1 for inc in self.standard_incs])
      code_block.dedent()
      code_block.write("}else{")
      code_block.indent()
      self.write_core(code_block, max_process_width, max_unroll, self.standard_incs)
      code_block.dedent()
      code_block.write("}")
      self.vec.max_into("(&max)", 0, 1, self.m_vars)
      code_block.write("return max;")

  def write_core(self, code_block, max_process_width, max_unroll, incs):
    code_block.new_line()
    def body(unroll, align = False):
      if type(unroll) == str:
        process_width = self.compute_process_width(self.vec.type_size)
        self.preprocess(code_block, self.vec.type_size, incs, partial=unroll)
        self.process(code_block, process_width)
      else:
        process_width = self.compute_process_width(unroll)
        self.preprocess(code_block, unroll, incs, align=align)
        self.process(code_block, process_width)
    self.vec.iterate_unrolled_aligned("i", "n", self.load_ptrs, incs, max_unroll, 1, body)

  def process(self, code_block, process_width):
    code_block.set_equal(itertools.cycle(self.m_vars), self.vec.max(itertools.cycle(self.m_vars), self.load_vars[0][:process_width]))

class AMax(AMaxP):
  standard_incs = ["incv"]

  def __init__(self, data_type_class):
    super(AMax, self).__init__(data_type_class)

  def write_declaration(self, code_block, settings):
    super(AMax, self).write_declaration(code_block, settings)
    code_block.write("{0} {1}amax(int n, {0}* v, int incv){{".format(self.data_type.name, self.data_type.name_char))

  def define_load_vars(self, code_block, process_width):
    self.load_vars = [["v_" + str(i) for i in range(process_width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])

  def define_load_ptrs(self, code_block, process_width):
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* v_base = (" + self.data_type.base_type.name + "*) v;")
      self.load_ptrs = ["v_base"]
    else:
      self.load_ptrs = ["v"]

  def compute_process_width(self, unroll):
    return (unroll * self.data_type.base_size)//self.vec.base_size

  def preprocess(self, code_block, unroll, incs, partial="", align = False):
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load(self.load_ptrs[0], 0, incs[0], unroll, align)))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial)))

class AMaxM(AMaxP):
  standard_incs = ["incv", "incy"]

  def __init__(self, data_type_class):
    super(AMaxM, self).__init__(data_type_class)

  def write_declaration(self, code_block, settings):
    super(AMaxM, self).write_declaration(code_block, settings)
    code_block.write("{0} {1}amaxm(int n, {0}* v, int incv, {0}* y, int incy){{".format(self.data_type.name, self.data_type.name_char))

  def define_load_ptrs(self, code_block, process_width):
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* v_base = (" + self.data_type.base_type.name + "*) v;")
      code_block.write(self.data_type.base_type.name + "* y_base = (" + self.data_type.base_type.name + "*) y;")
      self.load_ptrs = ["v_base", "y_base"]
    else:
      self.load_ptrs = ["v", "y"]

  def define_load_vars(self, code_block, process_width):
    if self.data_type.is_complex:
      self.load_vars = [["v_" + str(i) for i in range(process_width)], ["y_" + str(i) for i in range(process_width//2)]]
    else:
      self.load_vars = [["v_" + str(i) for i in range(process_width)], ["y_" + str(i) for i in range(process_width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])
    code_block.define_vars(self.vec.type_name, self.load_vars[1])

  def compute_process_width(self, unroll):
    return (unroll * self.data_type.base_size)//self.vec.base_size * self.data_type.base_size

  def preprocess(self, code_block, unroll, incs, partial="", align = False):
    process_width = self.compute_process_width(unroll)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll, align))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], unroll))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    if self.data_type.is_complex:
      code_block.set_equal(self.load_vars[0][process_width//2:process_width], self.vec.abs(self.vec.mul(self.vec.swap_pairwise(self.load_vars[0]), self.vec.rep_odds(self.load_vars[1]))))
      code_block.set_equal(self.load_vars[0][:process_width//2], self.vec.abs(self.vec.mul(self.load_vars[0], self.vec.rep_evens(self.load_vars[1]))))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.mul(self.load_vars[0], self.load_vars[1][:process_width])))
