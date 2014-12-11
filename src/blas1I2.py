import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
from utils import *
from dataTypes import *
from vectorizations import *
from generate import *

class OneDimensionalAccumulation(Function):

  def __init__(self, data_type_class):
    super(OneDimensionalAccumulation, self).__init__(data_type_class)

  def write_declaration(self, code_block, settings):
    code_block.srcFile.include("#include <stdio.h>")
    code_block.srcFile.include("#include <stdlib.h>")
    code_block.srcFile.include("#include <float.h>")
    code_block.srcFile.include("#include <math.h>")
    code_block.srcFile.include("#include \"config.h\"")
    code_block.srcFile.include("#include \"Common/Common.h\"")
    #code_block.srcFile.include("#include \"IndexedFP/" + self.data_type.base_type.name_char + "Indexed.h\"")
    #code_block.srcFile.include("#include \"rblas1.h\"")

  def write_body(self, code_block, settings = [(-1, 1, 8)]):
    code_block.write("SET_DAZ_FLAG;")
    if len(settings) == 1:
      self.write_cores(code_block, settings[0][0], settings[0][1], settings[0][2])
    else:
      code_block.write("switch(fold){")
      code_block.indent()
      for (fold, max_stride, max_unroll) in settings:
        if fold == -1:
          code_block.write("default:{")
        else:
          code_block.write("case " + str(fold) + ":{")
        code_block.indent()
        self.write_cores(code_block, fold, max_stride, max_unroll)
        #code_block.write("break;")
        code_block.write("RESET_DAZ_FLAG")
        code_block.write("return;")
        code_block.dedent()
        code_block.write("}")
      code_block.dedent()
      code_block.write("}")

  def write_cores(self, code_block, fold, max_stride, max_unroll):
#    if self.vec.name == "AVX":
#      self.code_block.write('printf("Hi im avx {0}\\n");'.format(self.name))
#    elif self.vec.name == "SSE":
#      self.code_block.write('printf("Hi im sse {0}\\n");'.format(self.name))
#    else:
#      self.code_block.write('printf("Hi im sisd {0}\\n");'.format(self.name))

# max_unroll_width is the maximum number of elements of the input vector that are processed in each iteration.
# unroll_width is the number of elements of the input vector that are processed in each iteration.
# process_width is the number of variables that are processed in each iteration.
    max_width = self.compute_width(max_stride)
    if fold == -1:
      code_block.write("int i, j;")
    else:
      code_block.write("int i;")
    code_block.new_line()
    sum_ptr = "sum"
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* sum_base = (" + self.data_type.base_type.name + "*) sum;")
      sum_ptr = "sum_base"
    self.define_load_ptrs(code_block, max_width * max_unroll)
    self.define_load_vars(code_block, max_width * max_unroll)
    if fold == -1:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(max_width)]
      code_block.define_vars(self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_" + str(i) for i in range(max_width)]]
      code_block.define_vars(self.vec.type_name, self.s_vars[0])
    else:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(max_width)]
      code_block.define_vars(self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_{0}_{1}".format(j, i) for i in range(max_width)] for j in range(fold)]
      for j in range(fold):
        code_block.define_vars(self.vec.type_name, self.s_vars[j])

    if fold == -1:
      code_block.write("{0} s_buffer[{1}];".format(self.vec.type_name, mix("*", max_width, "MAX_FOLD")))
      self.buffer_vars = ["s_buffer[{0}]".format(mix("+", mix("*", "j", max_width), i)) for i in range(max_width)]

    code_block.new_line()

    #propagate sum to buffer
    if fold == -1:
      code_block.write("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.propagate_into(self.buffer_vars, sum_ptr, "j", 1)
      code_block.dedent()
      code_block.write("}")
    else:
      for j in range(fold):
        self.vec.propagate_into(self.s_vars[j], sum_ptr, j, 1)

    code_block.write("if(" + " && ".join([inc + " == 1" for inc in self.standard_incs]) + "){")
    code_block.indent()
    self.write_core(code_block, fold, max_stride, max_unroll, [1 for inc in self.standard_incs])
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_stride, max_unroll, self.standard_incs)
    code_block.dedent()
    code_block.write("}")
    #consolidate
    if fold == -1:
      code_block.write("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.consolidate_into("sum", "j", 1, self.buffer_vars, sum_ptr, "j", 1, self.q_vars[0])
      code_block.dedent()
      code_block.write("}")
    else:
      for j in range(fold):
        self.vec.consolidate_into("sum", j, 1, self.s_vars[j], sum_ptr, j, 1, self.q_vars[0])

  def write_core(self, code_block, fold, max_stride, max_unroll, incs):
    max_width = self.compute_width(max_stride);
    code_block.new_line()
    def body(n, align = False):
      if type(n) == str:
        width = self.compute_width(self.vec.type_size)
        self.preprocess(code_block, self.vec.type_size, incs, partial=n, align=align)
        self.process(code_block, fold, width, 1)
      else:
        width = self.compute_width(min(n, max_stride))
        self.preprocess(code_block, n, incs, align=align)
        self.process(code_block, fold, width, n // max_stride)
    #self.vec.iterate_unrolled_aligned("i", "n", self.load_ptrs, incs, max_stride * max_unroll, 1, body)
    self.vec.iterate_unrolled("i", "n", self.load_ptrs, incs, max_stride * max_unroll, 1, body)

  def define_load_vars(self, code_block, width):
    raise(NotImplementedError())

  def define_load_ptrs(self, code_block, width):
    raise(NotImplementedError())

  def preprocess(self, code_block, n, incs, partial="", align = False):
    raise(NotImplementedError())

  def process(self, code_block, fold, width, unroll):
    if(fold == -1):
      code_block.write("for(j = 0; j < fold - 1; j++){")
      code_block.indent()
      for i in range(max(unroll, 1)):
        code_block.set_equal(self.s_vars[0], self.buffer_vars[:width])
        self.vec.add_BLP_into(self.q_vars, self.s_vars[0], self.load_vars[0][i * width:], width)
        code_block.set_equal(self.buffer_vars, self.q_vars[:width])
        code_block.set_equal(self.q_vars, self.vec.sub(self.s_vars[0], self.q_vars[:width]))
        code_block.set_equal(self.load_vars[0][i * width:], self.vec.add(self.load_vars[0][i * width:], self.q_vars[:width]))
      code_block.dedent()
      code_block.write("}")
      self.vec.add_BLP_into(self.buffer_vars, self.buffer_vars, self.load_vars[0], width)
    else:
      for i in range(max(unroll, 1)):
        for j in range(fold - 1):
          code_block.set_equal(self.q_vars, self.s_vars[j])
          self.vec.add_BLP_into(self.s_vars[j], self.s_vars[j], self.load_vars[0][i * width:], width)
          code_block.set_equal(self.q_vars, self.vec.sub(self.q_vars, self.s_vars[j][:width]))
          code_block.set_equal(self.load_vars[0][i * width:], self.vec.add(self.load_vars[0][i * width:], self.q_vars[:width]))
        self.vec.add_BLP_into(self.s_vars[fold - 1], self.s_vars[fold - 1], self.load_vars[0], width)

  def compute_width(self, n):
    raise(NotImplementedError())

class NonDotOneDimensionalAccumulation(OneDimensionalAccumulation):
  standard_incs = ["incv"]

  def __init__(self, data_type_class):
    super(NonDotOneDimensionalAccumulation, self).__init__(data_type_class)

  def define_load_vars(self, code_block, width):
    self.load_vars = [["v_" + str(i) for i in range(width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])

  def define_load_ptrs(self, code_block, width):
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* v_base = (" + self.data_type.base_type.name + "*) v;")
      self.load_ptrs = ["v_base"]
    else:
      self.load_ptrs = ["v"]

  def compute_width(self, n):
    return (n * self.data_type.base_size)//self.vec.base_size

class DotOneDimensionalAccumulation(OneDimensionalAccumulation):
  standard_incs = ["incv", "incy"]

  def __init__(self, data_type_class):
    super(DotOneDimensionalAccumulation, self).__init__(data_type_class)

  def define_load_ptrs(self, code_block, width):
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* v_base = (" + self.data_type.base_type.name + "*) v;")
      code_block.write(self.data_type.base_type.name + "* y_base = (" + self.data_type.base_type.name + "*) y;")
      self.load_ptrs = ["v_base", "y_base"]
    else:
      self.load_ptrs = ["v", "y"]

  def define_load_vars(self, code_block, width):
    if self.data_type.is_complex:
      self.load_vars = [["v_" + str(i) for i in range(width)], ["y_" + str(i) for i in range(width//2)]]
    else:
      self.load_vars = [["v_" + str(i) for i in range(width)], ["y_" + str(i) for i in range(width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])
    code_block.define_vars(self.vec.type_name, self.load_vars[1])

  def compute_width(self, n):
    return (n * self.data_type.base_size)//self.vec.base_size * self.data_type.base_size

