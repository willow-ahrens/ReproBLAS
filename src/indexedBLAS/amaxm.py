import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
from utils import *
from dataTypes import *
from vectorizations import *
from generate import *
import itertools
import amax

class AMaxM(amax.AMax):
  name = "amaxm"

  def __init__(self, data_type_class, N, X, incX, Y, incY, amaxm):
    super(AMaxM, self).__init__(data_type_class, N, X, incX, amaxm)
    self.Y = Y
    self.incY = incY
    self.standard_incs = [incX, incY]

  def define_load_ptrs(self, code_block, reg_width):
    self.load_ptrs = [self.X, self.Y]

  def define_load_vars(self, code_block, reg_width):
    if self.data_type.is_complex:
      self.load_vars = [["{}_{}".format(self.X, i) for i in range(reg_width)], ["{}_{}".format(self.Y, i) for i in range(reg_width//2)]]
    else:
      self.load_vars = [["{}_{}".format(self.X, i) for i in range(reg_width)], ["{}_{}".format(self.Y, i) for i in range(reg_width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])
    code_block.define_vars(self.vec.type_name, self.load_vars[1])

  def compute_reg_width(self, unroll_width):
    return (unroll_width * self.data_type.base_size)//self.vec.base_size * self.data_type.base_size

  def preprocess(self, code_block, unroll_width, incs, partial="", align = False):
    reg_width = self.compute_reg_width(unroll_width)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll_width, align))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], unroll_width))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    if self.data_type.is_complex:
      code_block.set_equal(self.load_vars[0][reg_width//2:reg_width], self.vec.abs(self.vec.mul(self.vec.swap_pairwise(self.load_vars[0]), self.vec.rep_odds(self.load_vars[1]))))
      code_block.set_equal(self.load_vars[0][:reg_width//2], self.vec.abs(self.vec.mul(self.load_vars[0], self.vec.rep_evens(self.load_vars[1]))))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.mul(self.load_vars[0], self.load_vars[1][:reg_width])))
