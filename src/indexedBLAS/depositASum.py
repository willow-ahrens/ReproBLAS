import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
from src.idxd import deposit

class DepositASum(deposit.Deposit):
  def __init__(self, data_type_class, fold_name, N_name, X_name, incX_name, manY_name, incmanY_name):
    super(DepositASum, self).__init__(data_type_class, fold_name, N_name, X_name, incX_name, manY_name, incmanY_name)
    redundant_char = ""
    if self.data_type_class.is_complex:
      redundant_char = self.data_type_class.base_type.name_char
    self.name = "{0}{1}depositASum".format(redundant_char, self.data_type_class.name_char)
    self.metric_name = "r{0}{1}asum".format(redundant_char, self.data_type_class.name_char)

  def preprocess(self, code_block, n, incs, partial="", align = False):
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load(self.load_ptrs[0], 0, incs[0], n, align)))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial)))
