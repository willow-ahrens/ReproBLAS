import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import blas1I2

class ASumI2(blas1I2.NonDotOneDimensionalAccumulation):
  def __init__(self, data_type_class):
    super(ASumI2, self).__init__(data_type_class)

  def write_declaration(self, code_block, settings):
    super(ASumI2, self).write_declaration(code_block, settings)
    redundant_char = ""
    if self.data_type.is_complex:
      redundant_char = self.data_type.base_type.name_char
    code_block.write("void {0}{1}asumI2(int n, {2}* v, int incv, int fold, {2}* sum){{".format(redundant_char, self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial="", align = False):
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load(self.load_ptrs[0], 0, incs[0], unroll, align)))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial)))
