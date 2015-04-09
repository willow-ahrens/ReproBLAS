import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import blas1I2

class ASumI2(blas1I2.NonDotOneDimensionalAccumulation):
  def __init__(self, data_type_class):
    super(ASumI2, self).__init__(data_type_class)
    redundant_char = ""
    if self.data_type_class.is_complex:
      redundant_char = self.data_type_class.base_type.name_char
    self.name = "{0}{1}asumI2".format(redundant_char, self.data_type_class.name_char)

  def get_metrics(self):
    return ["bench_{}".format(self.name[:-1])]

  def preprocess(self, code_block, n, incs, partial="", align = False):
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load(self.load_ptrs[0], 0, incs[0], n, align)))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial)))
