import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import blas1I2

class DotI2(blas1I2.DotOneDimensionalAccumulation):
  def __init__(self, data_type):
    assert not data_type.is_complex, "dot is only for real types"
    super(DotI2, self).__init__(data_type)

  def write_declaration(self, code_block, settings):
    super(DotI2, self).write_declaration(code_block, settings)
    code_block.write("void {0}dotI2(int n, {1}* v, int incv, {1}* y, int incy, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def define_preprocess_vars(self):
    return

  def preprocess(self, code_block, n, incs, partial="", align = False):
    width = self.compute_width(n)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], n, align))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], n))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    
    code_block.set_equal(self.load_vars[0], self.vec.mul(self.load_vars[0], self.load_vars[1][:width]))
