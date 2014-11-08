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

  def write_body(self, code_block, settings = [(-1, 1)]):
    code_block.write("SET_DAZ_FLAG;")
    if len(settings) == 1:
      self.write_cores(code_block,settings[0][0], settings[0][1])
    else:
      code_block.write("switch(fold){")
      code_block.indent()
      for (fold, max_unroll_width) in settings:
        if fold == -1:
          code_block.write("default:{")
        else:
          code_block.write("case " + str(fold) + ":{")
        code_block.indent()
        self.write_cores(code_block, fold, max_unroll_width)
        #code_block.write("break;")
        code_block.write("RESET_DAZ_FLAG")
        code_block.write("return;")
        code_block.dedent()
        code_block.write("}")
      code_block.dedent()
      code_block.write("}")

  def write_cores(self, code_block, fold, max_unroll):
#    if self.vec.name == "AVX":
#      self.code_block.write('printf("Hi im avx {0}\\n");'.format(self.name))
#    elif self.vec.name == "SSE":
#      self.code_block.write('printf("Hi im sse {0}\\n");'.format(self.name))
#    else:
#      self.code_block.write('printf("Hi im sisd {0}\\n");'.format(self.name))

# max_unroll_width is the maximum number of elements of the input vector that are processed in each iteration.
# unroll_width is the number of elements of the input vector that are processed in each iteration.
# process_width is the number of variables that are processed in each iteration.
    max_process_width = self.compute_process_width(max_unroll)
    if fold == -1:
      code_block.write("int i, j;")
    else:
      code_block.write("int i;")
    code_block.new_line()
    sum_ptr = "sum"
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* sum_base = (" + self.data_type.base_type.name + "*) sum;")
      sum_ptr = "sum_base"
    self.define_load_ptrs(code_block, max_process_width)
    self.define_load_vars(code_block, max_process_width)
    if fold == -1:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(max_process_width)]
      code_block.define_vars(self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_" + str(i) for i in range(max_process_width)]]
      code_block.define_vars(self.vec.type_name, self.s_vars[0])
    else:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(self.vec.suf_width)]
      code_block.define_vars(self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_{0}_{1}".format(j, i) for i in range(self.vec.suf_width)] for j in range(fold)]
      for j in range(fold):
        code_block.define_vars(self.vec.type_name, self.s_vars[j])

    if fold == -1:
      code_block.write("{0} s_buffer[{1}];".format(self.vec.type_name, mix("*", max_process_width, "MAX_FOLD")))
      self.buffer_vars = ["s_buffer[{0}]".format(mix("+", mix("*", "j", max_process_width), i)) for i in range(max_process_width)]

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
    self.write_core(code_block, fold, max_process_width, max_unroll, [1 for inc in self.standard_incs])
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_process_width, max_unroll, self.standard_incs)
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

  def write_core(self, code_block, fold, max_process_width, max_unroll, incs):
    code_block.new_line()
    def body(unroll, align = False):
      if type(unroll) == str:
        process_width = self.compute_process_width(self.vec.type_size)
        self.preprocess(code_block, self.vec.type_size, incs, partial=unroll)
        self.process(code_block, fold, process_width)
      else:
        process_width = self.compute_process_width(unroll)
        self.preprocess(code_block, unroll, incs, align=align)
        self.process(code_block, fold, process_width)
    self.vec.iterate_unrolled("i", "n", self.load_ptrs, incs, max_unroll, 1, body)

  def define_load_vars(self, code_block, process_width):
    raise(NotImplementedError())

  def define_load_ptrs(self, code_block, process_width):
    raise(NotImplementedError())

  def preprocess(self, code_block, unroll, incs, partial="", align = False):
    raise(NotImplementedError())

  def process(self, code_block, fold, process_width):
    if(fold == -1):
      code_block.write("for(j = 0; j < fold - 1; j++){")
      code_block.indent()
      code_block.set_equal(self.s_vars[0], self.buffer_vars[:process_width])
      self.vec.add_BLP_into(self.q_vars, self.s_vars[0], self.load_vars[0], process_width)
      code_block.set_equal(self.buffer_vars, self.q_vars[:process_width])
      code_block.set_equal(self.q_vars, self.vec.sub(self.s_vars[0], self.q_vars[:process_width]))
      code_block.set_equal(self.load_vars[0], self.vec.add(self.load_vars[0], self.q_vars[:process_width]))
      code_block.dedent()
      code_block.write("}")
      self.vec.add_BLP_into(self.buffer_vars, self.buffer_vars, self.load_vars[0], process_width)
    else:
      for i in range(process_width // self.vec.suf_width):
        for j in range(fold - 1):
          code_block.set_equal(self.q_vars, self.s_vars[j])
          self.vec.add_BLP_into(self.s_vars[j], self.s_vars[j], self.load_vars[0][i * self.vec.suf_width:], self.vec.suf_width)
          code_block.set_equal(self.q_vars, self.vec.sub(self.q_vars, self.s_vars[j]))
          code_block.set_equal(self.load_vars[0][i * self.vec.suf_width:], self.vec.add(self.load_vars[0][i * self.vec.suf_width:], self.q_vars))      
        self.vec.add_BLP_into(self.s_vars[fold - 1], self.s_vars[fold - 1], self.load_vars[0][i * self.vec.suf_width:], self.vec.suf_width)
    
  def compute_process_width(self, unroll):
    raise(NotImplementedError())

class NonDotOneDimensionalAccumulation(OneDimensionalAccumulation):
  standard_incs = ["incv"]

  def __init__(self, data_type_class):
    super(NonDotOneDimensionalAccumulation, self).__init__(data_type_class)

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

class DotOneDimensionalAccumulation(OneDimensionalAccumulation):
  standard_incs = ["incv", "incy"]

  def __init__(self, data_type_class):
    super(DotOneDimensionalAccumulation, self).__init__(data_type_class)

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
