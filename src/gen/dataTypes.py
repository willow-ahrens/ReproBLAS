################################################################################
# dataTypes.py                                                                 #
#                                                                              #
#     A set of classes used to represent the four main data types supported by #
# ReproBLAS.                                                                   #
#                                                                              #
#                                                            Peter Ahrens 2014 #
################################################################################

from utils import *

class DataType(object):
  name_char = ''
  name = "" #name of type
  int_float_name = ""
  int_char = ""
  float_char = ""
  base_type = None #base type (base type of double complex is double)
  byte_size = 0
  bit_size = 0
  base_size = 0
  is_complex = False

  def index(self, offset, inc, i, base_index = True):
    if(self.is_complex and base_index):
      return mix("+", mix("*", self.base_size, inc, mix("+", offset, mix("//", i, 2))), mix("%", i, 2))
    else:
      return mix("*", inc, mix("+", offset, i))

  def data_increment(self, src_ptrs, src_incs, i):
    return ", ".join(["{0} += {1}".format(src_ptr, mix("*", src_inc, self.base_size, i)) for (src_ptr, src_inc) in zip(src_ptrs, src_incs)])
 # consider here the possible optimization where we multiply incv once by two if the data type is complex.

class Double(DataType):
  name_char = 'd'
  name = "double"
  int_float_name = "l_double"
  int_char = "l"
  float_char = "d"
  byte_size = 8
  bit_size = 64
  base_size = 1

  def __init__(self, code):
    code.srcFile.include("#include <float.h>")
Double.base_type = Double

class DoubleComplex(DataType):
  name_char = 'z'
  name = "double complex"
  byte_size = Double.byte_size * 2
  bit_size = Double.bit_size * 2
  base_size = 2
  is_complex = True
  base_type = Double

  def __init__(self, code):
    code.srcFile.include("#include <complex.h>")

class Float(DataType):
  name_char = 's'
  name = "float"
  int_float_name = "i_float"
  int_char = "i"
  float_char = "f"
  byte_size = 4
  bit_size = 32
  base_size = 1

  def __init__(self, code):
    code.srcFile.include("#include <float.h>")
Float.base_type = Float


class FloatComplex(DataType):
  name_char = 'c'
  name = "float complex"
  byte_size = Float.byte_size * 2
  bit_size = Float.bit_size * 2
  base_size = 2
  is_complex = True
  base_type = Float

  def __init__(self, code):
    code.srcFile.include("#include <complex.h>")
