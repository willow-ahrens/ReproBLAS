import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import blas1I2

class SumI2(blas1I2.NonDotOneDimensionalAccumulation):
  name = "sumI"

  def __init__(self, data_type_class):
    super(SumI2, self).__init__(data_type_class)

  def write_declaration(self, code_block, settings):
    super(SumI2, self).write_declaration(code_block, settings)
    code_block.write("void {0}sumI2(int n, {1}* v, int incv, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
