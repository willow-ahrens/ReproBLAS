import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
from src.indexed import deposit

class DepositSum(deposit.Deposit):
  def __init__(self, data_type_class, N, X, incX, manY, incmanY):
    super(DepositSum, self).__init__(data_type_class, N, X, incX, manY, incmanY)
    self.name = "{0}depositSum".format(self.data_type_class.name_char)
    self.metric_name = "r{0}sum".format(self.data_type_class.name_char)

