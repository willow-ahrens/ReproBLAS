import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import blas1I2

class DotUI2(blas1I2.DotOneDimensionalAccumulation):
  def __init__(self, data_type_class):
    assert data_type_class.is_complex, "dotu is only for complex types"
    super(DotUI2, self).__init__(data_type_class)

  def write_declaration(self, code_block, settings):
    super(DotUI2, self).write_declaration(code_block, settings)
    code_block.write("void {0}dotuI2(int n, {1}* v, int incv, {1}* y, int incy, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    process_width = self.compute_process_width(unroll)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], unroll))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    code_block.set_equal(self.load_vars[0][process_width//2:process_width], self.vec.nconj(self.vec.mul(self.vec.swap_pairwise(self.load_vars[0]), self.vec.rep_odds(self.load_vars[1]))))
    code_block.set_equal(self.load_vars[0][:process_width//2], self.vec.mul(self.load_vars[0], self.vec.rep_evens(self.load_vars[1])))
