import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
from utils import *
from dataTypes import *
from vectorizations import *
from generate import *
import itertools

class AMax(Target):
  name = "amax"

  def __init__(self, data_type_class, N, X, incX, amax):
    super(AMax, self).__init__()
    self.data_type_class = data_type_class
    self.N = N
    self.X = X
    self.incX = incX
    self.amax = amax
    self.standard_incs = [incX]

  def get_arguments(self):
    return ["{}{}_max_unroll_width_{}".format(self.data_type_class.name_char, self.name, vectorization.name) for vectorization in vectorization_lookup.values()]

  def get_parameters(self):
    parameters = []
    for vectorization in vectorization_lookup.values():
      vec = vectorization(CodeBlock(), self.data_type_class)
      name = "{}{}_max_unroll_width_{}".format(self.data_type_class.name_char, self.name, vec.name)
      step = max(1, vec.type_size)
      minimum = step
      maximum = step * 8
      default = step
      parameters.append(IntegerParameter(name, minimum, maximum, step, default))
    return parameters

  def get_metrics(self):
    return {argument: ["bench_{}{}".format(self.data_type_class.name_char, self.name)] for argument in self.get_arguments()}

  def write(self, code_block):
    iterate_all_vectorizations(self.write_vec, code_block)

  def write_vec(self, vec_class, code_block):
      self.data_type = self.data_type_class(code_block)
      self.vec = vec_class(code_block, self.data_type_class)
      max_unroll_width = self.arguments["{}{}_max_unroll_width_{}".format(self.data_type.name_char, self.name, self.vec.name)]
      max_reg_width = self.compute_reg_width(max_unroll_width)
      code_block.write("int i;")
      code_block.new_line()
      self.define_load_ptrs(code_block, max_reg_width)
      self.define_load_vars(code_block, max_reg_width)
      self.m_vars = ["m_" + str(i) for i in range(self.vec.suf_width)]
      code_block.define_vars(self.vec.type_name, self.m_vars)
      code_block.set_equal(self.m_vars, itertools.repeat(self.vec.zero));

      code_block.new_line()

      code_block.write("if(" + " && ".join([inc + " == 1" for inc in self.standard_incs]) + "){")
      code_block.indent()
      self.write_core(code_block, max_reg_width, max_unroll_width, [1 for inc in self.standard_incs])
      code_block.dedent()
      code_block.write("}else{")
      code_block.indent()
      self.write_core(code_block, max_reg_width, max_unroll_width, self.standard_incs)
      code_block.dedent()
      code_block.write("}")
      self.vec.max_into(self.amax, 0, 1, self.m_vars)

  def write_core(self, code_block, max_reg_width, max_unroll_width, incs):
    code_block.new_line()
    def body(unroll_width, align = False):
      if type(unroll_width) == str:
        reg_width = self.compute_reg_width(self.vec.type_size)
        self.preprocess(code_block, self.vec.type_size, incs, partial=unroll_width)
        self.process(code_block, reg_width)
      else:
        reg_width = self.compute_reg_width(unroll_width)
        self.preprocess(code_block, unroll_width, incs, align=align)
        self.process(code_block, reg_width)
    self.vec.iterate_unrolled("i", self.N, self.load_ptrs, incs, max_unroll_width, 1, body)

  def preprocess(self, code_block, unroll_width, incs, partial="", align = False):
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load(self.load_ptrs[0], 0, incs[0], unroll_width, align)))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial)))

  def process(self, code_block, reg_width):
    code_block.set_equal(itertools.cycle(self.m_vars), self.vec.max(itertools.cycle(self.m_vars), self.load_vars[0][:reg_width]))

  def define_load_vars(self, code_block, reg_width):
    self.load_vars = [["{}_{}".format(self.X, i) for i in range(reg_width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])

  def define_load_ptrs(self, code_block, reg_width):
    self.load_ptrs = [self.X]

  def compute_reg_width(self, unroll_width):
    return (unroll_width * self.data_type.base_size)//self.vec.base_size
