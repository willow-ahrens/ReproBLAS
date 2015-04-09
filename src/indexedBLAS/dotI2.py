import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import blas1I2

class DotI2(blas1I2.DotOneDimensionalAccumulation):
  def __init__(self, data_type_class):
    assert not data_type_class.is_complex, "dot is only for real types"
    super(DotI2, self).__init__(data_type_class)
    self.name = "{0}dotI2".format(self.data_type_class.name_char, self.data_type_class.name)
    self.metric_name = "r{0}dot".format(self.data_type_class.name_char, self.data_type_class.name)

  def define_preprocess_vars(self):
    return

  #TODO vec.load should just handle the partial case separately.
  def preprocess(self, code_block, n, incs, partial="", align = False):
    reg_width = self.compute_reg_width(n)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], n, align))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], n))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    code_block.set_equal(self.load_vars[0], self.vec.mul(self.load_vars[0], self.load_vars[1][:reg_width]))
