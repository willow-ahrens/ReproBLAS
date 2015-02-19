import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import blas1I2

class DotCI2(blas1I2.DotOneDimensionalAccumulation):
  def __init__(self, data_type_class):
    assert data_type_class.is_complex, "dotc is only for complex types"
    super(DotCI2, self).__init__(data_type_class)
    self.name = "{0}dotcI2".format(self.data_type.name_char)

  def preprocess(self, code_block, n, incs, partial="", align = False):
    width = self.compute_width(n)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], n, align))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], n))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    code_block.set_equal(self.load_vars[0][width//2:width], self.vec.conj(self.vec.mul(self.vec.swap_pairwise(self.load_vars[0]), self.vec.rep_odds(self.load_vars[1]))))
    code_block.set_equal(self.load_vars[0][:width//2], self.vec.mul(self.load_vars[0], self.vec.rep_evens(self.load_vars[1])))
