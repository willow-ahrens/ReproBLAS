import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
from utils import *
from dataTypes import *
from vectorizations import *
from generate import *

class OneDimensionalAccumulation(Target):

  def __init__(self, data_type_class, vec_class):
    self.data_type_class = data_type_class
    self.vec_class = vec_class


  #SUM_WIDTH = number of indexed sums used at once
  #PIPE_WIDTH = number of independently loaded input elements processed per indexed sum
  #REG_WIDTH = number of variables needed to hold the independently loaded elements
  #UNROLL_WIDTH = number of times PIPE_WIDTH elements are to be processed in the inner loop
  def write(self, code_block, arguments):
    self.data_type = self.data_type_class(code_block)
    self.vec = self.vec_class(code_block, self.data_type_class)
    code_block.write("SET_DAZ_FLAG;")
    expanded_folds = []
    for i in range(1, 9):
      if arguments["{}_expand_fold_{}".format(self.name, i)]:
        expanded_folds.append(i)
    expanded_folds.append(0)
    if len(expanded_folds) == 1:
      self.write_fold(code_block, 0, arguments["{}_max_pipe_width_{}_fold_{}".format(self.name, self.vec.name, 0)], arguments["{}_max_unroll_width_{}_fold_{}".format(self.name, self.vec.name, 0)])
    else:
      code_block.write("switch(fold){")
      code_block.indent()
      for fold in expanded_folds:
        if fold == 0:
          code_block.write("default:{")
        else:
          code_block.write("case " + str(fold) + ":{")
        code_block.indent()
        self.write_fold(code_block, fold, arguments["{}_max_pipe_width_{}_fold_{}".format(self.name, self.vec.name, fold)], arguments["{}_max_unroll_width_{}_fold_{}".format(self.name, self.vec.name, fold)])
        #code_block.write("break;")
        code_block.write("RESET_DAZ_FLAG")
        code_block.write("return;")
        code_block.dedent()
        code_block.write("}")
      code_block.dedent()
      code_block.write("}")

  #TODO this should probably be two methods, one for generic fold and one for specific
  def write_fold(self, code_block, fold, max_pipe_width, max_unroll_width):
#    if self.vec.name == "AVX":
#      self.code_block.write('printf("Hi im avx {0}\\n");'.format(self.name))
#    elif self.vec.name == "SSE":
#      self.code_block.write('printf("Hi im sse {0}\\n");'.format(self.name))
#    else:
#      self.code_block.write('printf("Hi im sisd {0}\\n");'.format(self.name))

    max_reg_width = self.compute_reg_width(max_pipe_width)
    if fold == 0:
      code_block.write("int i, j;")
    else:
      code_block.write("int i;")
    code_block.new_line()
    sum_ptr = "sum"
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* sum_base = (" + self.data_type.base_type.name + "*) sum;")
      sum_ptr = "sum_base"
    self.define_load_ptrs(code_block, max_reg_width * max_unroll_width)
    self.define_load_vars(code_block, max_reg_width * max_unroll_width)
    if fold == 0:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(max_reg_width)]
      code_block.define_vars(self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_" + str(i) for i in range(max_reg_width)]]
      code_block.define_vars(self.vec.type_name, self.s_vars[0])
    else:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(max_reg_width)]
      code_block.define_vars(self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_{0}_{1}".format(j, i) for i in range(max_reg_width)] for j in range(fold)]
      for j in range(fold):
        code_block.define_vars(self.vec.type_name, self.s_vars[j])

    if fold == 0:
      code_block.write("{0} s_buffer[{1}];".format(self.vec.type_name, mix("*", max_reg_width, "MAX_FOLD")))
      self.buffer_vars = ["s_buffer[{0}]".format(mix("+", mix("*", "j", max_reg_width), i)) for i in range(max_reg_width)]

    code_block.new_line()

    #propagate sum to buffer
    if fold == 0:
      code_block.write("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.propagate_into(self.buffer_vars, sum_ptr, "j", 1)
      code_block.dedent()
      code_block.write("}")
    else:
      for j in range(fold):
        self.vec.propagate_into(self.s_vars[j], sum_ptr, j, 1)

    self.write_cores(code_block, fold, max_pipe_width, max_unroll_width)

    #consolidate
    if fold == 0:
      code_block.write("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.consolidate_into("sum", "j", 1, self.buffer_vars, sum_ptr, "j", 1, self.q_vars[0])
      code_block.dedent()
      code_block.write("}")
    else:
      for j in range(fold):
        #TODO consolidate_into should be a subtract and a reduce
        self.vec.consolidate_into("sum", j, 1, self.s_vars[j], sum_ptr, j, 1, self.q_vars[0])

  def write_core(self, code_block, fold, max_pipe_width, max_unroll_width, incs):
    max_reg_width = self.compute_reg_width(max_pipe_width);
    code_block.new_line()
    def body(n, align = False):
      if type(n) == str:
        reg_width = self.compute_reg_width(self.vec.type_size)
        self.preprocess(code_block, self.vec.type_size, incs, partial=n, align=align)
        self.process(code_block, fold, reg_width, 1)
      else:
        reg_width = self.compute_reg_width(min(n, max_pipe_width))
        self.preprocess(code_block, n, incs, align=align)
        self.process(code_block, fold, reg_width, n // max_pipe_width)
    #self.vec.iterate_unrolled_aligned("i", "n", self.load_ptrs, incs, max_pipe_width * max_unroll_width, 1, body)
    self.vec.iterate_unrolled("i", "n", self.load_ptrs, incs, max_pipe_width * max_unroll_width, 1, body)

  def define_load_vars(self, code_block, width):
    raise(NotImplementedError())

  def define_load_ptrs(self, code_block, width):
    raise(NotImplementedError())

  def preprocess(self, code_block, n, incs, partial="", align = False):
    raise(NotImplementedError())

  def process(self, code_block, fold, reg_width, unroll_width):
    if(fold == 0):
      code_block.write("for(j = 0; j < fold - 1; j++){")
      code_block.indent()
      for i in range(max(unroll_width, 1)):
        code_block.set_equal(self.s_vars[0], self.buffer_vars[:reg_width])
        self.vec.add_BLP_into(self.q_vars, self.s_vars[0], self.load_vars[0][i * reg_width:], reg_width)
        code_block.set_equal(self.buffer_vars, self.q_vars[:reg_width])
        code_block.set_equal(self.q_vars, self.vec.sub(self.s_vars[0], self.q_vars[:reg_width]))
        code_block.set_equal(self.load_vars[0][i * reg_width:], self.vec.add(self.load_vars[0][i * reg_width:], self.q_vars[:reg_width]))
      code_block.dedent()
      code_block.write("}")
      self.vec.add_BLP_into(self.buffer_vars, self.buffer_vars, self.load_vars[0][i * reg_width:], reg_width)
    else:
      for i in range(max(unroll_width, 1)):
        for j in range(fold - 1):
          code_block.set_equal(self.q_vars, self.s_vars[j][:reg_width])
          self.vec.add_BLP_into(self.s_vars[j], self.s_vars[j], self.load_vars[0][i * reg_width:], reg_width)
          code_block.set_equal(self.q_vars, self.vec.sub(self.q_vars, self.s_vars[j][:reg_width]))
          code_block.set_equal(self.load_vars[0][i * reg_width:], self.vec.add(self.load_vars[0][i * reg_width:], self.q_vars[:reg_width]))
        self.vec.add_BLP_into(self.s_vars[fold - 1], self.s_vars[fold - 1], self.load_vars[0][i * reg_width:], reg_width)

  def compute_reg_width(self, pipe_width):
    raise(NotImplementedError())

class NonDotOneDimensionalAccumulation(OneDimensionalAccumulation):
  standard_incs = ["incv"]

  def __init__(self, data_type_class, vec_class):
    super(NonDotOneDimensionalAccumulation, self).__init__(data_type_class, vec_class)

  def write_cores(self, code_block, fold, max_pipe_width, max_unroll_width):
    code_block.write("if(incv == 1){")
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, [1])
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, self.standard_incs)
    code_block.dedent()
    code_block.write("}")

  def define_load_vars(self, code_block, width):
    self.load_vars = [["v_" + str(i) for i in range(width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])

  def define_load_ptrs(self, code_block, width):
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* v_base = (" + self.data_type.base_type.name + "*) v;")
      self.load_ptrs = ["v_base"]
    else:
      self.load_ptrs = ["v"]

  def compute_reg_width(self, pipe_width):
    return (pipe_width * self.data_type.base_size)//self.vec.base_size

class DotOneDimensionalAccumulation(OneDimensionalAccumulation):
  standard_incs = ["incv", "incy"]

  def __init__(self, data_type_class, vec_class):
    super(DotOneDimensionalAccumulation, self).__init__(data_type_class, vec_class)

  def write_cores(self, code_block, fold, max_pipe_width, max_unroll_width):
    code_block.write("if(incv == 1){")
    code_block.indent()
    code_block.write("if(incy == 1){")
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, [1, 1])
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, [1, "incy"])
    code_block.dedent()
    code_block.write("}")
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    code_block.write("if(incy == 1){")
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, ["incv", 1])
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, ["incv", "incy"])
    code_block.dedent()
    code_block.write("}")
    code_block.dedent()
    code_block.write("}")

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

  def compute_reg_width(self, pipe_width):
    return (pipe_width * self.data_type.base_size)//self.vec.base_size * self.data_type.base_size

