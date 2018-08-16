import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
from src.binned import deposit

class DepositSum(deposit.Deposit):
  def __init__(self, data_type_class, fold_name, N_name, X_name, incX_name, manY_name, incmanY_name):
    super(DepositSum, self).__init__(data_type_class, fold_name, N_name, X_name, incX_name, manY_name, incmanY_name)
    self.name = "{0}depositSum".format(self.data_type_class.name_char)
    self.metric_name = "r{0}sum".format(self.data_type_class.name_char)

