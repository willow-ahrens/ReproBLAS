import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import depositM

class DepositDotC(depositM.DepositM):
  def __init__(self, data_type_class, fold_name, N_name, X_name, incX_name, manY_name, incmanY_name, Z_name, incZ_name):
    assert data_type_class.is_complex, "dotc is only for complex types"
    super(DepositDotC, self).__init__(data_type_class, fold_name, N_name, X_name, incX_name, manY_name, incmanY_name, Z_name, incZ_name)
    self.name = "{0}depositDotC".format(self.data_type_class.name_char)
    self.metric_name = "r{0}dotc".format(self.data_type_class.name_char)

  def preprocess(self, code_block, n, incs, partial="", align = False):
    reg_width = self.compute_reg_width(n)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], n, align))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], n))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    code_block.set_equal(self.load_vars[0][reg_width//2:reg_width], self.vec.mul(self.vec.swap_pairwise(self.load_vars[0]), self.vec.rep_odds(self.load_vars[1])))
    code_block.set_equal(self.load_vars[0][:reg_width//2], self.vec.conj(self.vec.mul(self.load_vars[0], self.vec.rep_evens(self.load_vars[1]))))
