import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
from utils import *
from dataTypes import *
from vectorizations import *
from generate import *
from scripts import terminal
import itertools

  #SUM_WIDTH = number of indexed sums used at once
  #PIPE_WIDTH = number of independently loaded input elements processed per indexed sum
  #REG_WIDTH = number of variables needed to hold the independently loaded elements
  #UNROLL_WIDTH = number of times PIPE_WIDTH elements per indexed sum are to be processed in the inner loop
class Deposit(Target):
  def __init__(self, data_type_class, N, X, incX, manY, incmanY):
    super(Deposit, self).__init__()
    self.default_fold = terminal.get_default_fold()
    self.max_fold = terminal.get_max_fold()
    self.max_expand_fold = min(self.max_fold, 6)
    self.data_type_class = data_type_class
    self.N = N
    self.X = X
    self.incX = incX
    self.manY = manY
    self.incmanY = incmanY

  def get_arguments(self):
    arguments = []
    for i in range(self.max_expand_fold + 1):
      for vectorization in vectorization_lookup.values():
        if(i != 0):
          arguments.append("{}_expand_{}_fold_{}".format(self.name, vectorization.name, i))
        arguments.append("{}_max_pipe_width_{}_fold_{}".format(self.name, vectorization.name, i))
        arguments.append("{}_max_unroll_width_{}_fold_{}".format(self.name, vectorization.name, i))
    return arguments

  def get_metrics(self):
    metrics = {}
    for i in range(self.max_expand_fold + 1):
      for vectorization in vectorization_lookup.values():
        if(i != 0):
          metrics["{}_expand_{}_fold_{}".format(self.name, vectorization.name, i)] = ["bench_{}_fold_{}".format(self.metric_name, i)]
        metrics["{}_max_pipe_width_{}_fold_{}".format(self.name, vectorization.name, i)] = ["bench_{}_fold_{}".format(self.metric_name, i)]
        metrics["{}_max_unroll_width_{}_fold_{}".format(self.name, vectorization.name, i)] = ["bench_{}_fold_{}".format(self.metric_name, i)]
    return metrics

  def get_parameters(self):
    parameters = []
    for i in range(self.max_expand_fold + 1):
      for vectorization in vectorization_lookup.values():
        vec = vectorization(CodeBlock(), self.data_type_class)
        if(i != 0):
          parameters.append(BooleanParameter("{}_expand_{}_fold_{}".format(self.name, vec.name, i), {"vectorization":vec.name}, i == self.default_fold))
        parameters.append(IntegerParameter("{}_max_unroll_width_{}_fold_{}".format(self.name, vec.name, i), {"vectorization":vec.name}, 1, 16, 1, 1))
        name = "{}_max_pipe_width_{}_fold_{}".format(self.name, vec.name, i)
        minimum = max(1, vec.type_size)
        maximum = minimum * 16
        default = minimum
        parameters.append(PowerOfTwoParameter(name, {"vectorization":vec.name}, minimum, maximum, default))
    return parameters

  def write(self, code_block):
    code_block.write("{")
    code_block.indent()
    self.code_block = code_block.sub_block()
    iterate_all_vectorizations(self.write_vec, self.code_block)
    code_block.indent()
    code_block.write("}")

  def write_vec(self, vec_class, code_block):
    self.data_type = self.data_type_class(code_block)
    self.vec = vec_class(code_block, self.data_type_class)
    self.vec.set_SIMD_daz_ftz()
    code_block.new_line()
    expanded_folds = []
    for i in range(1, self.max_expand_fold + 1):
      if self.arguments["{}_expand_{}_fold_{}".format(self.name, self.vec.name, i)]:
        expanded_folds.append(i)
    expanded_folds.append(0)
    if len(expanded_folds) == 1:
      self.write_fold(code_block, 0, self.arguments["{}_max_pipe_width_{}_fold_{}".format(self.name, self.vec.name, 0)], self.arguments["{}_max_unroll_width_{}_fold_{}".format(self.name, self.vec.name, 0)])
    else:
      code_block.write("switch(fold){")
      code_block.indent()
      for fold in expanded_folds:
        if fold == 0:
          code_block.write("default:")
        else:
          code_block.write("case " + str(fold) + ":")
        code_block.indent()
        code_block.write("{")
        code_block.indent()
        self.write_fold(code_block, fold, self.arguments["{}_max_pipe_width_{}_fold_{}".format(self.name, self.vec.name, fold)], self.arguments["{}_max_unroll_width_{}_fold_{}".format(self.name, self.vec.name, fold)])
        code_block.new_line()
        self.vec.reset_SIMD_daz_ftz()
        code_block.dedent()
        code_block.write("}")
        code_block.write("break;")
        code_block.dedent()
      code_block.dedent()
      code_block.write("}")

  def write_fold(self, code_block, fold, max_pipe_width, max_unroll_width):
    max_reg_width = self.compute_reg_width(max_pipe_width)
    if fold == 0:
      code_block.write("int i, j;")
    else:
      code_block.write("int i;")
    self.define_load_ptrs(code_block)
    self.define_load_vars(code_block, max_reg_width * max_unroll_width)
    self.compression_vars = ["compression_" + str(i) for i in range(self.vec.suf_width)]
    self.expansion_vars = ["expansion_" + str(i) for i in range(self.vec.suf_width)]
    code_block.define_vars(self.vec.type_name, self.compression_vars)
    code_block.define_vars(self.vec.type_name, self.expansion_vars)
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
      self.buffer0_vars = ["s_buffer[{0}]".format(mix("+", mix("*", "0", max_reg_width), i)) for i in range(max_reg_width)]

    code_block.new_line()

    #propagate sum to buffer
    if fold == 0:
      code_block.write("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.propagate_into(self.buffer_vars, self.manY, "j", self.incmanY)
      code_block.dedent()
      code_block.write("}")
    else:
      for j in range(fold):
        self.vec.propagate_into(self.s_vars[j], self.manY, j, self.incmanY)

    code_block.new_line()

    self.write_increments(code_block, fold, max_pipe_width, max_unroll_width)

    code_block.new_line()

    #consolidate
    if fold == 0:
      code_block.write("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.consolidate_into(self.manY, "j", self.incmanY, self.buffer_vars, self.manY, "j", self.incmanY)
      code_block.dedent()
      code_block.write("}")
    else:
      for j in range(fold):
        self.vec.consolidate_into(self.manY, j, self.incmanY, self.s_vars[j], self.manY, j, self.incmanY)

  def write_increments(self, code_block, fold, max_pipe_width, max_unroll_width):
    code_block.write("if({} == 1){{".format(self.incX))
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, [1])
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, [self.incX])
    code_block.dedent()
    code_block.write("}")

  def write_core(self, code_block, fold, max_pipe_width, max_unroll_width, incs):
    max_reg_width = self.compute_reg_width(max_pipe_width);

    def body(n, align = False):
      if type(n) == str:
        reg_width = self.compute_reg_width(self.vec.type_size)
        self.preprocess(code_block, self.vec.type_size, incs, partial=n, align=align)
        code_block.new_line()
        self.process(code_block, fold, reg_width, 1)
      else:
        reg_width = self.compute_reg_width(min(n, max_pipe_width))
        self.preprocess(code_block, n, incs, align=align)
        code_block.new_line()
        self.process(code_block, fold, reg_width, n // max_pipe_width)

    def body0(n, align = False):
      if type(n) == str:
        reg_width = self.compute_reg_width(self.vec.type_size)
        self.preprocess(code_block, self.vec.type_size, incs, partial=n, align=align)
        code_block.new_line()
        self.process0(code_block, fold, reg_width, 1)
      else:
        reg_width = self.compute_reg_width(min(n, max_pipe_width))
        self.preprocess(code_block, n, incs, align=align)
        code_block.new_line()
        self.process0(code_block, fold, reg_width, n // max_pipe_width)

    if self.data_type.is_complex:
      code_block.write("if({0}mindex0({1}) || {0}mindex0({1} + 1)){{".format(self.data_type.base_type.name_char, self.manY))
      code_block.indent()

      code_block.write("if({0}mindex0({1})){{".format(self.data_type.base_type.name_char, self.manY))
      code_block.indent()
      code_block.write("if({0}mindex0({1} + 1)){{".format(self.data_type.base_type.name_char, self.manY))
      code_block.indent()
      code_block.set_equal(self.compression_vars, self.vec.set("{0}MCOMPRESSION".format(self.data_type.base_type.name_char.upper())))
      code_block.set_equal(self.expansion_vars, self.vec.set("{0}MEXPANSION".format(self.data_type.base_type.name_char.upper())))
      code_block.dedent()
      code_block.write("}else{")
      code_block.indent()
      code_block.set_equal(self.compression_vars, self.vec.set_real_imag("{0}MCOMPRESSION".format(self.data_type.base_type.name_char.upper()), "1.0"))
      code_block.set_equal(self.expansion_vars, self.vec.set_real_imag("{0}MEXPANSION".format(self.data_type.base_type.name_char.upper()), "1.0"))
      code_block.dedent()
      code_block.write("}")
      code_block.dedent()
      code_block.write("}else{")
      code_block.indent()
      code_block.set_equal(self.compression_vars, self.vec.set_real_imag("1.0", "{0}MCOMPRESSION".format(self.data_type.base_type.name_char.upper())))
      code_block.set_equal(self.expansion_vars, self.vec.set_real_imag("1.0", "{0}MEXPANSION".format(self.data_type.base_type.name_char.upper())))
      code_block.dedent()
      code_block.write("}")
      self.vec.iterate_unrolled("i", self.N, self.load_ptrs, incs, max_pipe_width * max_unroll_width, 1, body0)
      code_block.dedent()
      code_block.write("}else{")
      code_block.indent()
      self.vec.iterate_unrolled("i", self.N, self.load_ptrs, incs, max_pipe_width * max_unroll_width, 1, body)
      code_block.dedent()
      code_block.write("}")
    else:
      code_block.write("if({0}mindex0({1})){{".format(self.data_type.base_type.name_char, self.manY))
      code_block.indent()
      code_block.set_equal(self.compression_vars, self.vec.set("{0}MCOMPRESSION".format(self.data_type.base_type.name_char.upper())))
      code_block.set_equal(self.expansion_vars, self.vec.set("{0}MEXPANSION".format(self.data_type.base_type.name_char.upper())))
      self.vec.iterate_unrolled("i", self.N, self.load_ptrs, incs, max_pipe_width * max_unroll_width, 1, body0)
      code_block.dedent()
      code_block.write("}else{")
      code_block.indent()
      self.vec.iterate_unrolled("i", self.N, self.load_ptrs, incs, max_pipe_width * max_unroll_width, 1, body)
      code_block.dedent()
      code_block.write("}")

  def preprocess(self, code_block, n, incs, partial="", align = False):
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], n, align))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))

  def process(self, code_block, fold, reg_width, unroll_width):
    if(fold == 0):
      for i in range(max(unroll_width, 1)):
        code_block.write("for(j = 0; j < fold - 1; j++){")
        code_block.indent()
        code_block.set_equal(self.s_vars[0], self.buffer_vars[:reg_width])
        self.vec.add_blp_into(self.q_vars, self.s_vars[0], self.load_vars[0][i * reg_width:], reg_width)
        code_block.set_equal(self.buffer_vars, self.q_vars[:reg_width])
        code_block.set_equal(self.q_vars, self.vec.sub(self.s_vars[0], self.q_vars[:reg_width]))
        code_block.set_equal(self.load_vars[0][i * reg_width:], self.vec.add(self.load_vars[0][i * reg_width:], self.q_vars[:reg_width]))
        code_block.dedent()
        code_block.write("}")
        self.vec.add_blp_into(self.buffer_vars, self.buffer_vars, self.load_vars[0][i * reg_width:], reg_width)
    else:
      for i in range(max(unroll_width, 1)):
        for j in range(fold - 1):
          code_block.set_equal(self.q_vars, self.s_vars[j][:reg_width])
          self.vec.add_blp_into(self.s_vars[j], self.s_vars[j], self.load_vars[0][i * reg_width:], reg_width)
          code_block.set_equal(self.q_vars, self.vec.sub(self.q_vars, self.s_vars[j][:reg_width]))
          code_block.set_equal(self.load_vars[0][i * reg_width:], self.vec.add(self.load_vars[0][i * reg_width:], self.q_vars[:reg_width]))
        self.vec.add_blp_into(self.s_vars[fold - 1], self.s_vars[fold - 1], self.load_vars[0][i * reg_width:], reg_width)

  def process0(self, code_block, fold, reg_width, unroll_width):
    if(fold == 0):
      for i in range(max(unroll_width, 1)):
        code_block.set_equal(self.s_vars[0], self.buffer0_vars[:reg_width])
        self.vec.add_blp_into(self.q_vars, self.s_vars[0], self.vec.mul(self.load_vars[0][i * reg_width:], itertools.cycle(self.compression_vars)), reg_width)
        code_block.set_equal(self.buffer0_vars, self.q_vars[:reg_width])
        code_block.write("if(fold > 1){") #TODO we know this at compile time
        code_block.indent()
        code_block.set_equal(self.q_vars, self.vec.sub(self.s_vars[0], self.q_vars[:reg_width]))
        code_block.set_equal(self.load_vars[0][i * reg_width:], self.vec.add(self.load_vars[0][i * reg_width:], self.vec.mul(self.q_vars[:reg_width], itertools.cycle(self.expansion_vars))))
        code_block.write("for(j = 1; j < fold - 1; j++){")
        code_block.indent()
        code_block.set_equal(self.s_vars[0], self.buffer_vars[:reg_width])
        self.vec.add_blp_into(self.q_vars, self.s_vars[0], self.load_vars[0][i * reg_width:], reg_width)
        code_block.set_equal(self.buffer_vars, self.q_vars[:reg_width])
        code_block.set_equal(self.q_vars, self.vec.sub(self.s_vars[0], self.q_vars[:reg_width]))
        code_block.set_equal(self.load_vars[0][i * reg_width:], self.vec.add(self.load_vars[0][i * reg_width:], self.q_vars[:reg_width]))
        code_block.dedent()
        code_block.write("}")
        self.vec.add_blp_into(self.buffer_vars, self.buffer_vars, self.load_vars[0][i * reg_width:], reg_width)
        code_block.dedent()
        code_block.write("}")
    else:
      for i in range(max(unroll_width, 1)):
        if(fold == 1):
          self.vec.add_blp_into(self.s_vars[0], self.s_vars[0], self.vec.mul(self.load_vars[0][i * reg_width:], itertools.cycle(self.compression_vars)), reg_width)
        else:
          code_block.set_equal(self.q_vars, self.s_vars[0][:reg_width])
          self.vec.add_blp_into(self.s_vars[0], self.s_vars[0], self.vec.mul(self.load_vars[0][i * reg_width:], itertools.cycle(self.compression_vars)), reg_width)
          code_block.set_equal(self.q_vars, self.vec.sub(self.q_vars, self.s_vars[0][:reg_width]))
          code_block.set_equal(self.load_vars[0][i * reg_width:], self.vec.add(self.load_vars[0][i * reg_width:], self.vec.mul(self.q_vars[:reg_width], itertools.cycle(self.expansion_vars))))
          for j in range(1, fold - 1):
            code_block.set_equal(self.q_vars, self.s_vars[j][:reg_width])
            self.vec.add_blp_into(self.s_vars[j], self.s_vars[j], self.load_vars[0][i * reg_width:], reg_width)
            code_block.set_equal(self.q_vars, self.vec.sub(self.q_vars, self.s_vars[j][:reg_width]))
            code_block.set_equal(self.load_vars[0][i * reg_width:], self.vec.add(self.load_vars[0][i * reg_width:], self.q_vars[:reg_width]))
          self.vec.add_blp_into(self.s_vars[fold - 1], self.s_vars[fold - 1], self.load_vars[0][i * reg_width:], reg_width)

  def define_load_vars(self, code_block, reg_width):
    self.load_vars = [["{}_{}".format(self.X, i) for i in range(reg_width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])

  def define_load_ptrs(self, code_block):
    self.load_ptrs = [self.X]

  def compute_reg_width(self, pipe_width):
    return (pipe_width * self.data_type.base_size)//self.vec.base_size
