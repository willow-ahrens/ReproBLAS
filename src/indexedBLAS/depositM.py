import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
from utils import *
from dataTypes import *
from vectorizations import *
from generate import *
from src.indexed import deposit

class DepositM(deposit.Deposit):
  def __init__(self, data_type_class, N, X, incX, manY, incmanY, Z, incZ):
    super(DepositM, self).__init__(data_type_class, N, X, incX, manY, incmanY)
    self.Z = Z
    self.incZ = incZ

  def write_increments(self, code_block, fold, max_pipe_width, max_unroll_width):
    code_block.write("if({} == 1){{".format(self.incX))
    code_block.indent()
    code_block.write("if({} == 1){{".format(self.incZ))
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, [1, 1])
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, [1, self.incZ])
    code_block.dedent()
    code_block.write("}")
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    code_block.write("if({} == 1){{".format(self.incZ))
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, [self.incX, 1])
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_pipe_width, max_unroll_width, [self.incX, self.incZ])
    code_block.dedent()
    code_block.write("}")
    code_block.dedent()
    code_block.write("}")

  def define_load_ptrs(self, code_block, width):
      self.load_ptrs = [self.X, self.Z]

  def define_load_vars(self, code_block, width):
    if self.data_type.is_complex:
      self.load_vars = [["{}_{}".format(self.X, i) for i in range(width)], ["{}_{}".format(self.Z, i) for i in range(width//2)]]
    else:
      self.load_vars = [["{}_{}".format(self.X, i) for i in range(width)], ["{}_{}".format(self.Z, i) for i in range(width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])
    code_block.define_vars(self.vec.type_name, self.load_vars[1])

  def compute_reg_width(self, pipe_width):
    return (pipe_width * self.data_type.base_size)//self.vec.base_size * self.data_type.base_size

