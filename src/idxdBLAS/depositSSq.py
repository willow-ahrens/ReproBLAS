import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
from src.idxd import deposit

class DepositSSq(deposit.Deposit):
  def __init__(self, data_type_class, fold_name, N_name, X_name, incX_name, manY_name, incmanY_name, scale_name):
    super(DepositSSq, self).__init__(data_type_class, fold_name, N_name, X_name, incX_name, manY_name, incmanY_name)
    redundant_char = ""
    if self.data_type_class.is_complex:
      redundant_char = self.data_type_class.base_type.name_char
    self.name = "{0}{1}depositSSq".format(redundant_char, self.data_type_class.name_char)
    self.metric_name = "r{0}{1}nrm2".format(redundant_char, self.data_type_class.name_char)
    self.scale_name = scale_name

  def preprocess(self, code_block, n, incs, partial="", align=False):
    code_block.include("{0} scale_mask = {1};".format(self.vec.type_name, self.vec.set(self.scale_name)[0]))
    reg_width = self.compute_reg_width(n)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.div(self.vec.load(self.load_ptrs[0], 0, incs[0], n, align), ["scale_mask"] * reg_width))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.div(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial), ["scale_mask"] * reg_width))
    code_block.set_equal(self.load_vars[0], self.vec.mul(self.load_vars[0], self.load_vars[0][:reg_width]))

  def set_daz_ftz(self, code_block):
    code_block.write("if(!{}mdenorm({}, {}) && {} >= {}_MIN){{".format(self.data_type.name_char, self.fold_name, self.priY_name, self.scale_name, "DBL" if self.data_type.base_type.name == "double" else "FLT"))
    code_block.indent()
    self.vec.set_SIMD_daz_ftz()
    code_block.dedent()
    code_block.write("}")
