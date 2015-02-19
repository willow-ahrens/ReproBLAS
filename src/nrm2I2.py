import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import blas1I2

class Nrm2I2(blas1I2.NonDotOneDimensionalAccumulation):
  def __init__(self, data_type_class):
    super(Nrm2I2, self).__init__(data_type_class)
    redundant_char = ""
    if self.data_type.is_complex:
      redundant_char = self.data_type.base_type.name_char
    self.name = "{0}{1}nrm2I2".format(redundant_char, self.data_type.name_char)

  def preprocess(self, code_block, n, incs, partial="", align=False):
    code_block.include("{0} scale_mask = {1};".format(self.vec.type_name, self.vec.set("scale")[0]))
    width = self.compute_width(n)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.mul(self.vec.load(self.load_ptrs[0], 0, incs[0], n, align), ["scale_mask"] * width))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.mul(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial), ["scale_mask"] * width))
    code_block.set_equal(self.load_vars[0], self.vec.mul(self.load_vars[0], self.load_vars[0][:width]))
