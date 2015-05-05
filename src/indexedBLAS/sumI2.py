import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import blas1I2

class SumI2(blas1I2.Deposit):
  name = "sumI"

  def __init__(self, data_type_class, N, X, incX, manY, incmanY):
    super(SumI2, self).__init__(data_type_class, N, X, incX, manY, incmanY)
    self.name = "{0}sumI2".format(self.data_type_class.name_char)
    self.metric_name = "r{0}sum".format(self.data_type_class.name_char)

