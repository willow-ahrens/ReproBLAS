#Generator.py
#a file used to generate the reproducible blas routines
#by Peter Ahrens
#comments and reformatting to come.

def mix(op, *args, paren = True):
  int_args = [str(arg) for arg in args if str(arg).isdigit()]
  str_args = [str(arg) for arg in args if not str(arg).isdigit()]
  if int_args:
    int_result = eval(" {0} ".format(op).join(int_args))
    if not str_args:
      return int_result
    identities = {"+" : {0}, "*" : {1}, "%" : {}, "//" : {}}
    if int_result not in identities[op]:
      str_args += [str(int_result)]
    zeros = {"+" : {}, "*" : {0}, "%" : {}, "//" : {}}
    if int_result in zeros[op]:
      return 0
  str_result = " {0} ".format(op).join(str_args)
  if paren and len(str_args) > 1:
    str_result = "({0})".format(str_result)
  return str_result

class DataType:
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
  
#todo naming conventions! make l_double into long_double
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
    code.srcFile.include_line("#include <float.h>")
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
    code.srcFile.include_line("#include <complex.h>")

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
    code.srcFile.include_line("#include <float.h>")
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
    code.srcFile.include_line("#include <complex.h>")

class Vectorization:
  name = ""
  defined_macro = ""

  def __init__(self, code_block, data_type_class):
    self.code_block = code_block
    self.data_type = data_type_class(code_block)

  def define_unroll_step(self, fold, unroll):
    code_block.write_line("#define " + self.unroll_step_string(fold) + " " + str(unroll))

  def define_max_unroll_step(self, fold, max_unroll_width):
    code_block.write_line("#define MAX_" + self.unroll_step_string(fold) + " " + str(max_unroll_width))

  def unroll_step_string(self, fold):
    if fold == -1:
      return "UNROLL_STEP_" + self.name
    else:
      return "UNROLL_STEP_" + self.name + "_" + str(fold) + "-FOLD"

  def consolidate_into(self, dst_ptr, offset, inc, src_vars, common_summand_ptr, common_summand_offset, common_summand_inc): 
    raise(NotImplementedError())

  #propagates the value pointed to by pointer to the variables listed in variables
  def propagate_into(self, dst_vars, src_ptr, offset, inc):
    raise(NotImplementedError())

  def add_BLP_into(dst, src, blp, width):
    raise(NotImplementedError())

  #loads n values from pointer into the variables listed in variables (n must be a multiple of type_size)
  def load(self, src_ptr, offset, inc, n):
    raise(NotImplementedError())

  #loads less than type_size values from pointer into the variables listed in variables
  def load_partial(self, src_ptr, offset, inc, n):
    raise(NotImplementedError())

  def sub(self, src_vars, amt_vars):
    raise(NotImplementedError())

  def add(self, src_vars, amt_vars):
    raise(NotImplementedError())

  def mul(self, src_vars, amt_vars):
    raise(NotImplementedError())

  def abs(self, src_vars):
    raise(NotImplementedError())

  def conj(self, src_vars):
    raise(NotImplementedError())

  def set(self, src_var):
    raise(NotImplementedError())

  def rep_evens(self, src_vars):
    raise(NotImplementedError())

  def rep_odds(self, src_vars):
    raise(NotImplementedError())

  def swap_pairwise(self, src_vars):
    raise(NotImplementedError())

class SISD(Vectorization):
  name = "SISD"

  def __init__(self, code_block, data_type_class):
    super().__init__(code_block, data_type_class)
    self.code_block.srcFile.include_line("#include <emmintrin.h>") #needed for DAZ_FLAG. Unclear if the flag is even necessary in SISD.
    self.bit_size = self.data_type.base_type.bit_size
    self.byte_size = self.data_type.base_type.byte_size
    self.type_name = self.data_type.base_type.name
    self.base_size = self.bit_size//self.data_type.base_type.bit_size
    self.type_size = self.bit_size//self.data_type.bit_size #notice that type_size is 0 in SISD vectorization if the data type is complex. It is probably better not to use type_size for SISD purposes.
    self.zero = 0
    self.suf_width = self.data_type.base_size


  def consolidate_into(self, dst_ptr, offset, inc, src_vars, common_summand_ptr, common_summand_offset, common_summand_inc, common_summand_var): 
    if self.data_type.is_complex:
      if(len(src_vars) > 2):
        self.code_block.write_line("{0} = (({1}*){2})[{3}];".format(common_summand_var, self.type_name, dst_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
      for src_var in src_vars[2::2]:
        self.code_block.write_line("{0} = {0} + ({1} - {2});".format(src_vars[0], src_var, common_summand_var))
      if(len(src_vars) > 2):
        self.code_block.write_line("{0} = (({1}*){2})[{3}];".format(common_summand_var, self.type_name, dst_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 1)))
      for src_var in src_vars[3::2]:
        self.code_block.write_line("{0} = {0} + ({1} - {2});".format(src_vars[1], src_var, common_summand_var))
      self.code_block.write_line("(({0}*){1})[{2}] = {3};".format(self.type_name, dst_ptr, self.data_type.index(offset, inc, 0), src_vars[0]))
      self.code_block.write_line("(({0}*){1})[{2}] = {3};".format(self.type_name, dst_ptr, self.data_type.index(offset, inc, 1), src_vars[1]))
    else:
      if(len(src_vars) > 1):
        self.code_block.write_line("{0} = {1}[{2}];".format(common_summand_var, dst_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
      for src_var in src_vars[1:]:
        self.code_block.write_line("{0} = {0} + ({1} - {2});".format(src_vars[0], src_var, common_summand_var))
      self.code_block.write_line("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0), src_vars[0]))

  def propagate_into(self, dst_vars, src_ptr, offset, inc):
    if self.data_type.is_complex:
      assert len(dst_vars) % 2 == 0, "cannot propagate complex value to odd number of base type dst_vars"
      self.code_block.write_line(" = ".join(dst_vars[0::2]) + " = {0}[{1}];".format(src_ptr, self.data_type.index(offset, inc, 0))) 
      self.code_block.write_line(" = ".join(dst_vars[1::2]) + " = {0}[{1}];".format(src_ptr, self.data_type.index(offset, inc, 1))) 
    else:
      self.code_block.write_line(" = ".join(dst_vars) + " = {0}[{1}];".format(src_ptr, self.data_type.index(offset, inc, 0)))

  def include_BLP_vars(self):
    self.code_block.include_line("{0} tmp_BLP;".format(self.data_type.base_type.int_float_name))

  def add_BLP_into(self, dst, src, blp, width):
    assert len(dst) >= width
    assert len(src) >= width
    assert len(blp) >= width

    self.include_BLP_vars()
    for i in range(width):
      self.code_block.write_line("tmp_BLP.{0} = {1};".format(self.data_type.base_type.float_char, blp[i]))
      self.code_block.write_line("tmp_BLP.{0} |= 1;".format(self.data_type.base_type.int_char))
      self.code_block.write_line("{0} = {1} + tmp_BLP.{2};".format(dst[i], src[i], self.data_type.base_type.float_char))

  def load(self, src_ptr, offset, inc, n):
    assert n > 0, "n must be nonzero"
    assert inc != 0, "inc must be nonzero"

    return ["{0}[{1}]".format(src_ptr, self.data_type.index(offset, inc, i)) for i in range(n * self.data_type.base_size)]

  def load_partial(self, src_ptr, offset, inc, n):
    if(isinstance(n, int)):
      assert n > 0, "n must be nonzero"
      assert n < self.type_size, "n must be less than the number of types that fit in a vector"
    assert False, "there is no way to partially fill SISD vectors"

  def set(self, src_var):
    return [src_var]

  def sub(self, src_vars, amt_vars):
    return ["{0} - {1}".format(src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def add(self, src_vars, amt_vars):
    return ["{0} + {1}".format(src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def abs(self, src_vars):
    return ["fabs({0})".format(src_var) for src_var in src_vars]

  def mul(self, src_vars, amt_vars):
    return ["{0} * {1}".format(src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def rep_evens(self, src_vars):
    return [src_var for src_var in src_vars[::2] for _ in (0, 1)]

  def rep_odds(self, src_vars):
    return [src_var for src_var in src_vars[1::2] for _ in (0, 1)]

  def swap_pairwise(self, src_vars):
    assert len(src_vars) % 2 == 0, "If you want to swap pairwise, you need an even number of things to swap pairwise."
    return [src_vars[((i + 1) % 2) + (2 * (i // 2))] for i in range(len(src_vars))]

  def conj(self, src_vars):
    if self.data_type.is_complex:
      return [src_var if i % 2 == 0 else self.mul([src_var], ["-1"])[0] for (i, src_var) in enumerate(src_vars)]
    else:
      return src_vars

  def nconj(self, src_vars):
    if self.data_type.is_complex:
      return [src_var if i % 2 == 1 else self.mul([src_var], ["-1"])[0] for (i, src_var) in enumerate(src_vars)]
    else:
      return src_vars

class SIMD(Vectorization):

  def __init__(self, code_block, data_type_class):
    super().__init__(code_block, data_type_class)
    self.suf_width = 1

  def include_consolidation_vars(self):
    self.code_block.include_line("{0} tmp[{1}] __attribute__((aligned({2})));".format(self.data_type.name, self.type_size, self.byte_size))

  def include_ABS_vars(self):
    self.code_block.include_line("{0} mask_ABS; {1}_ABS_MASK{2}(mask_ABS);".format(self.type_name, self.name, self.data_type.base_type.name_char.upper()))

  def include_BLP_vars(self):
    self.code_block.include_line("{0} mask_BLP; {1}_BLP_MASK{2}(mask_BLP);".format(self.type_name, self.name, self.data_type.base_type.name_char.upper()))

  def include_CONJ_vars(self):
    self.code_block.include_line("{0} mask_CONJ; {1}_CONJ_MASK{2}(mask_CONJ);".format(self.type_name, self.name, self.data_type.base_type.name_char.upper())) 

  def include_NCONJ_vars(self):
    self.code_block.include_line("{0} mask_NCONJ; {1}_NCONJ_MASK{2}(mask_NCONJ);".format(self.type_name, self.name, self.data_type.base_type.name_char.upper()))


class SSE(SIMD):
  name = "SSE"
  defined_macro = "__SSE2__"

  def __init__(self, code_block, data_type_class):
    super().__init__(code_block, data_type_class)
    self.code_block.srcFile.include_line("#include <emmintrin.h>")
    self.bit_size = 128
    self.byte_size = 16
    self.type_name = {"float": "__m128", "double": "__m128d"}[self.data_type.base_type.name]
    self.base_size = self.bit_size//self.data_type.base_type.bit_size
    self.type_size = self.bit_size//self.data_type.bit_size
    self.zero = "_mm_setzero_p{0}()".format(self.data_type.base_type.name_char)

  def consolidate_into(self, dst_ptr, offset, inc, src_vars, common_summand_ptr, common_summand_offset, common_summand_inc, common_summand_var): 
    self.include_consolidation_vars()
    if self.data_type.name == "float":
      self.code_block.write_line("{0} = _mm_sub_ps({0}, _mm_set_ps({1}[{2}], {1}[{2}], {1}[{2}], 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
    elif self.data_type.name == "double":
      self.code_block.write_line("{0} = _mm_sub_pd({0}, _mm_set_pd({1}[{2}], 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
    elif self.data_type.name == "float complex":
      self.code_block.write_line("{0} = _mm_sub_ps({0}, _mm_set_ps({1}[{3}], {1}[{2}], 0, 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0), self.data_type.index(common_summand_offset, common_summand_inc, 1)))
    if len(src_vars) > 1:
      self.propagate_into([common_summand_var], common_summand_ptr, common_summand_offset, common_summand_inc)
      for src_var in src_vars[1:]:
        self.code_block.write_line("{0} = _mm_add_p{1}({0}, _mm_sub_p{1}({2}, {3}));".format(src_vars[0], self.data_type.base_type.name_char, src_var, common_summand_var))
    if self.data_type.name == "double complex":
      self.code_block.write_line("_mm_store_pd((double*){0}, {1});".format(mix("+", dst_ptr, self.data_type.index(offset, inc, 0), paren=False), src_vars[0]))
    else:
      if self.data_type.name == "float complex":
        self.code_block.write_line("_mm_store_ps((float*)tmp, {0});".format(src_vars[0]))
      else:
        self.code_block.write_line("_mm_store_p{0}(tmp, {1});".format(self.data_type.name_char, src_vars[0]))
      self.code_block.write_line("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0, False), " + ".join(["tmp[{0}]".format(i) for i in range(self.type_size)])))

  def propagate_into(self, dst_vars, src_ptr, offset, inc):
    if self.data_type.is_complex:
      broadcast = {"float": "(__m128)_mm_load1_pd((double *)({0}));", "double": "_mm_loadu_pd({0});"}[self.data_type.base_type.name]
    else:
      broadcast = "_mm_load1_p{0}({{0}});".format(self.data_type.name_char)
    self.code_block.write_line(" = ".join(dst_vars) + " = " + broadcast.format(mix("+", src_ptr, self.data_type.index(offset, inc, 0), paren=False)))

  def add_BLP_into(self, dst, src, blp, width):
    assert len(dst) >= width
    assert len(src) >= width
    assert len(blp) >= width

    self.include_BLP_vars()
    for i in range(width):
      self.code_block.write_line("{0} = _mm_add_p{1}({2}, _mm_or_p{1}({3}, mask_BLP));".format(dst[i], self.data_type.base_type.name_char, src[i], blp[i]))

  def load(self, src_ptr, offset, inc, n):
    assert n > 0, "n must be nonzero"
    assert n % self.type_size == 0, "n must be a multiple of the number of types that fit in a vector"
    result = []
    for i in range(n//self.type_size):
      if inc == 1 or self.type_size == 1:
        result += ["_mm_loadu_p{0}({1})".format(self.data_type.base_type.name_char, mix("+", src_ptr, self.data_type.index(offset, inc, i * self.base_size), paren=False))]
      else:
        if self.data_type.base_type.name == "float":
          result += ["_mm_set_ps({0}[{4}], {0}[{3}], {0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, i * self.base_size), self.data_type.index(offset, inc, i * self.base_size + 1), self.data_type.index(offset, inc, i * self.base_size + 2), self.data_type.index(offset, inc, i * self.base_size + 3))]
        if self.data_type.base_type.name == "double":
          result += ["_mm_set_pd({0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, i * self.base_size), self.data_type.index(offset, inc, i * self.base_size + 1))]
    return result

  def load_partial(self, src_ptr, offset, inc, n):
    if(isinstance(n, int)):
      assert n > 0, "n must be nonzero"
      assert n < self.type_size, "n must be less than the number of types that fit in a vector"
    if self.data_type.name == "float complex":
      return ["_mm_set_ps(0, 0, {0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, 0), self.data_type.index(offset, inc, 1))]
    elif self.data_type.name == "double":
      return ["_mm_set_pd(0, {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, 0))]
    elif self.data_type.name == "float":
      return ["_mm_set_ps(0, {1}>2?{0}[{4}]:0, {1}>1?{0}[{3}]:0, {0}[{2}])".format(src_ptr, n, self.data_type.index(offset, inc, 0), self.data_type.index(offset, inc, 1), self.data_type.index(offset, inc, 2))]

  def sub(self, src_vars, amt_vars):
    return ["_mm_sub_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def add(self, src_vars, amt_vars):
    return ["_mm_add_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def abs(self, src_vars):
    self.include_ABS_vars()
    return ["_mm_and_p{0}({1}, mask_ABS)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]

  def mul(self, src_vars, amt_vars):
    return ["_mm_mul_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def conj(self, src_vars):
    if self.data_type.is_complex:
      self.include_CONJ_vars()
      return ["_mm_xor_p{0}({1}, mask_CONJ)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]
    else:
      return src_vars

  def nconj(self, src_vars):
    if self.data_type.is_complex:
      self.include_NCONJ_vars()
      return ["_mm_xor_p{0}({1}, mask_NCONJ)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]
    else:
      return src_vars

  def set(self, src_var):
    return ["_mm_set1_p{0}({1})".format(self.data_type.base_type.name_char, src_var)]

  def rep_evens(self, src_vars):
    if self.data_type.base_type.name == "double":
      return ["_mm_shuffle_pd({0}, {0}, 0b00)".format(src_var) for src_var in src_vars]
    elif self.data_type.base_type.name == "float":
      return ["_mm_shuffle_ps({0}, {0}, 0b10100000)".format(src_var) for src_var in src_vars]

  def rep_odds(self, src_vars):
    if self.data_type.base_type.name == "double":
      return ["_mm_shuffle_pd({0}, {0}, 0b11)".format(src_var) for src_var in src_vars]
    elif self.data_type.base_type.name == "float":
      return ["_mm_shuffle_ps({0}, {0}, 0b11110101)".format(src_var) for src_var in src_vars]

  def swap_pairwise(self, src_vars):
    if self.data_type.base_type.name == "double":
      return ["_mm_shuffle_pd({0}, {0}, 0b01)".format(src_var) for src_var in src_vars]
    elif self.data_type.base_type.name == "float":
      return ["_mm_shuffle_ps({0}, {0}, 0b10110001)".format(src_var) for src_var in src_vars]
    
class AVX(SIMD):
  name = "AVX"
  defined_macro = "__AVX__"

  def __init__(self, code_block, data_type_class):
    super().__init__(code_block, data_type_class)
    self.code_block.srcFile.include_line("#include <immintrin.h>")
    self.bit_size = 256
    self.byte_size = 32 
    self.type_name = {"float": "__m256", "double": "__m256d"}[self.data_type.base_type.name]
    self.base_size = self.bit_size//self.data_type.base_type.bit_size
    self.type_size = self.bit_size//self.data_type.bit_size
    self.zero = "_mm256_setzero_p{0}()".format(self.data_type.base_type.name_char)

  def consolidate_into(self, dst_ptr, offset, inc, src_vars, common_summand_ptr, common_summand_offset, common_summand_inc, common_summand_var): 
    self.include_consolidation_vars()
    if self.data_type.name == "float":
      self.code_block.write_line("{0} = _mm256_sub_ps({0}, _mm256_set_ps({1}[{2}], {1}[{2}], {1}[{2}], {1}[{2}], {1}[{2}], {1}[{2}], {1}[{2}], 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
    elif self.data_type.name == "double":
      self.code_block.write_line("{0} = _mm256_sub_pd({0}, _mm256_set_pd({1}[{2}], {1}[{2}], {1}[{2}], 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
    elif self.data_type.name == "float complex":
      self.code_block.write_line("{0} = _mm256_sub_ps({0}, _mm256_set_ps({1}[{3}], {1}[{2}], {1}[{3}], {1}[{2}], {1}[{3}], {1}[{2}], 0, 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0), self.data_type.index(common_summand_offset, common_summand_inc, 1)))
    elif self.data_type.name == "double complex":
      self.code_block.write_line("{0} = _mm256_sub_pd({0}, _mm256_set_pd({1}[{3}], {1}[{2}], 0, 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0), self.data_type.index(common_summand_offset, common_summand_inc, 1)))
    if len(src_vars) > 1:
      self.propagate_into([common_summand_var], common_summand_ptr, common_summand_offset, common_summand_inc)
      for src_var in src_vars[1:]:
        self.code_block.write_line("{0} = _mm256_add_p{1}({0}, _mm256_sub_p{1}({2}, {3}));".format(src_vars[0], self.data_type.base_type.name_char, src_var, common_summand_var))
    if self.data_type.is_complex:
      self.code_block.write_line("_mm256_store_p{0}(({1}*)tmp, {2});".format(self.data_type.base_type.name_char, self.data_type.base_type.name, src_vars[0]))
    else:
      self.code_block.write_line("_mm256_store_p{0}(tmp, {1});".format(self.data_type.name_char, src_vars[0]))
    self.code_block.write_line("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0, False), " + ".join(["tmp[{0}]".format(i) for i in range(self.type_size)])))

  def propagate_into(self, dst_vars, src_ptr, offset, inc):
    if self.data_type.is_complex:
      broadcast = {"float": "(__m256)_mm256_broadcast_sd((double *)({0}));", "double": "_mm256_broadcast_pd((__m128d *)({0}));"}[self.data_type.base_type.name]
    else:
      broadcast = "_mm256_broadcast_s{0}({{0}});".format(self.data_type.name_char)
    self.code_block.write_line(" = ".join(dst_vars) + " = " + broadcast.format(mix("+", src_ptr, self.data_type.index(offset, inc, 0), paren=False)))

  def add_BLP_into(self, dst, src, blp, width):
    assert len(dst) >= width
    assert len(src) >= width
    assert len(blp) >= width
    self.include_BLP_vars()
    for i in range(width):
      self.code_block.write_line("{0} = _mm256_add_p{1}({2}, _mm256_or_p{1}({3}, mask_BLP));".format(dst[i], self.data_type.base_type.name_char, src[i], blp[i]))

  def load(self, src_ptr, offset, inc, n):
    assert n > 0, "n must be nonzero"
    assert n % self.type_size == 0, "n must be a multiple of the number of types that fit in a vector"
    result = []
    for i in range(n//self.type_size):
      if inc == 1 or self.type_size == 1:
        result += ["_mm256_loadu_p{0}({1})".format(self.data_type.base_type.name_char, mix("+", src_ptr, self.data_type.index(offset, inc, i * self.base_size), paren=False))]
      else:
        if self.data_type.base_type.name == "float":
          result += ["_mm256_set_ps({0}[{8}], {0}[{7}], {0}[{6}], {0}[{5}], {0}[{4}], {0}[{3}], {0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, i * self.base_size), self.data_type.index(offset, inc, i * self.base_size + 1), self.data_type.index(offset, inc, i * self.base_size + 2), self.data_type.index(offset, inc, i * self.base_size + 3), self.data_type.index(offset, inc, i * self.base_size + 4), self.data_type.index(offset, inc, i * self.base_size + 5), self.data_type.index(offset, inc, i * self.base_size + 6), self.data_type.index(offset, inc, i * self.base_size + 7))]
        elif self.data_type.base_type.name == "double":
          result += ["_mm256_set_pd({0}[{4}], {0}[{3}], {0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, i * self.base_size), self.data_type.index(offset, inc, i * self.base_size + 1), self.data_type.index(offset, inc, i * self.base_size + 2), self.data_type.index(offset, inc, i * self.base_size + 3))]
    return result

  def load_partial(self, src_ptr, offset, inc, n):
    if(isinstance(n, int)):
      assert n > 0, "n must be nonzero"
      assert n < self.type_size, "n must be less than the number of types that fit in a vector"
    if self.data_type.name == "double complex":
      return ["_mm256_set_pd(0, 0, {0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, 0), self.data_type.index(offset, inc, 1))]
    elif self.data_type.name == "float complex":
      return ["(__m256)_mm256_set_pd(0, {1}>2?((double*){0})[{4}]:0, {1}>1?((double*){0})[{3}]:0, ((double*){0})[{2}])".format(src_ptr, n, self.data_type.index(offset, inc, 0, False), self.data_type.index(offset, inc, 1, False), self.data_type.index(offset, inc, 2, False))]
    elif self.data_type.name == "double":
      return ["_mm256_set_pd(0, {1}>2?{0}[{4}]:0, {1}>1?{0}[{3}]:0, {0}[{2}])".format(src_ptr, n, self.data_type.index(offset, inc, 0), self.data_type.index(offset, inc, 1), self.data_type.index(offset, inc, 2))]
    elif self.data_type.name == "float":
      return ["_mm256_set_ps(0, {1}>6?{0}[{8}]:0, {1}>5?{0}[{7}]:0, {1}>4?{0}[{6}]:0, {1}>3?{0}[{5}]:0, {1}>2?{0}[{4}]:0, {1}>1?{0}[{3}]:0, {0}[{2}])".format(src_ptr, n, self.data_type.index(offset, inc, 0), self.data_type.index(offset, inc, 1), self.data_type.index(offset, inc, 2), self.data_type.index(offset, inc, 3), self.data_type.index(offset, inc, 4), self.data_type.index(offset, inc, 5), self.data_type.index(offset, inc, 6))]

  def sub(self, src_vars, amt_vars):
    return ["_mm256_sub_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def add(self, src_vars, amt_vars):
    return ["_mm256_add_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def abs(self, src_vars):
    self.include_ABS_vars()
    return ["_mm256_and_p{0}({1}, mask_ABS)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]

  def mul(self, src_vars, amt_vars):
    return ["_mm256_mul_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def conj(self, src_vars):
    if self.data_type.is_complex:
      self.include_CONJ_vars()
      return ["_mm256_xor_p{0}({1}, mask_CONJ)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]
    else:
      return src_vars

  def nconj(self, src_vars):
    if self.data_type.is_complex:
      self.include_NCONJ_vars()
      return ["_mm256_xor_p{0}({1}, mask_NCONJ)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]
    else:
      return src_vars
    
  def set(self, src_var):
    return ["_mm256_set1_p{0}({1})".format(self.data_type.base_type.name_char, src_var)]

  def rep_evens(self, src_vars):
    if self.data_type.base_type.name == "double":
      return ["_mm256_permute_pd({0}, 0b0000)".format(src_var) for src_var in src_vars]
    elif self.data_type.base_type.name == "float":
      return ["_mm256_permute_ps({0}, 0b10100000)".format(src_var) for src_var in src_vars]

  def rep_odds(self, src_vars):
    if self.data_type.base_type.name == "double":
      return ["_mm256_permute_pd({0}, 0b1111)".format(src_var) for src_var in src_vars]
    elif self.data_type.base_type.name == "float":
      return ["_mm256_permute_ps({0}, 0b11110101)".format(src_var) for src_var in src_vars]

  def swap_pairwise(self, src_vars):
    if self.data_type.base_type.name == "double":
      return ["_mm256_permute_pd({0}, 0b0101)".format(src_var) for src_var in src_vars]
    elif self.data_type.base_type.name == "float":
      return ["_mm256_permute_ps({0}, 0b10110001)".format(src_var) for src_var in src_vars]

class Function:
  name = ""

  def __init__(self, data_type_class, vectorization_class):
    self.data_type_class = data_type_class
    self.vectorization_class = vectorization_class

  def write(self, code_block):
    self.data_type = self.data_type_class(code_block)
    self.vec = self.vectorization_class(code_block, self.data_type_class)
    self.write_declaration(code_block)
    code_block.indent()
    self.data_type = self.data_type_class(code_block)
    body_block = code_block.sub_block()
    self.vec = self.vectorization_class(body_block, self.data_type_class)
    self.write_body(body_block)
    code_block.dedent()
    code_block.write_line("}")

  def write_declaration(self, code_block):
    raise(NotImplementedError())

  def write_body(self, code_block):
    raise(NotImplementedError())

class TargetFunction(Function):

  def __init__(self, data_type_class, vectorization_class):
    super().__init__(data_type_class, vectorization_class)

  @classmethod
  def file_name(cls, data_type_class):
    raise(NotImplementedError())

class OneDimensionalAccumulation(TargetFunction):

  def __init__(self, data_type_class, vectorization_class, settings = [(-1, 1)]):
    super().__init__(data_type_class, vectorization_class)
    self.set_settings(settings)

  def set_settings(self, settings):
    #put settings validation here
    self.settings = settings

  def write_declaration(self, code_block):
    code_block.srcFile.include_line("#include <stdio.h>")
    code_block.srcFile.include_line("#include <stdlib.h>")
    code_block.srcFile.include_line("#include <float.h>")
    code_block.srcFile.include_line("#include <math.h>")
    code_block.srcFile.include_line("#include \"config.h\"")
    code_block.srcFile.include_line("#include \"Common/Common.h\"")
    #code_block.srcFile.include_line("#include \"IndexedFP/" + self.data_type.base_type.name_char + "Indexed.h\"")
    #code_block.srcFile.include_line("#include \"rblas1.h\"")

  def write_body(self, code_block):
    code_block.write_line("SET_DAZ_FLAG;")
    if len(self.settings) == 1:
      self.write_cores(-1, self.settings[0][1])
      code_block.write_line("break;")
    else:
      code_block.write_line("switch(fold){")
      code_block.indent()
      for (fold, max_unroll_width) in self.settings:
        if fold == -1:
          code_block.write_line("default:{")
        else:
          code_block.write_line("case " + str(fold) + ":{")
        code_block.indent()
        self.write_cores(code_block, fold, max_unroll_width)
        #code_block.write_line("break;")
        code_block.write_line("RESET_DAZ_FLAG")
        code_block.write_line("return;")
        code_block.dedent()
        code_block.write_line("}")
      code_block.dedent()
      code_block.write_line("}")

  def write_cores(self, code_block, fold, max_unroll):
#    if self.vec.name == "AVX":
#      self.code_block.write_line('printf("Hi im avx {0}\\n");'.format(self.name))
#    elif self.vec.name == "SSE":
#      self.code_block.write_line('printf("Hi im sse {0}\\n");'.format(self.name))
#    else:
#      self.code_block.write_line('printf("Hi im sisd {0}\\n");'.format(self.name))

# max_unroll_width is the maximum number of elements of the input vector that are processed in each iteration.
# unroll_width is the number of elements of the input vector that are processed in each iteration.
# process_width is the number of variables that are processed in each iteration.
    max_process_width = self.compute_process_width(max_unroll)
    if fold == -1:
      code_block.write_line("int i, j;")
    else:
      code_block.write_line("int i;")
    code_block.new_line()
    sum_ptr = "sum"
    if self.data_type.is_complex:
      code_block.write_line(self.data_type.base_type.name + "* sum_base = (" + self.data_type.base_type.name + "*) sum;")
      sum_ptr = "sum_base"
    self.define_load_ptrs(code_block, max_process_width)
    self.define_load_vars(code_block, max_process_width)
    if fold == -1:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(max_process_width)]
      self.define_vars(code_block, self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_" + str(i) for i in range(max_process_width)]]
      self.define_vars(code_block, self.vec.type_name, self.s_vars[0])
    else:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(self.vec.suf_width)]
      self.define_vars(code_block, self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_{0}_{1}".format(j, i) for i in range(self.vec.suf_width)] for j in range(fold)]
      for j in range(fold):
        self.define_vars(code_block, self.vec.type_name, self.s_vars[j])

    if fold == -1:
      code_block.write_line("{0} s_buffer[{1}];".format(self.vec.type_name, mix("*", max_process_width, "MAX_FOLD")))
      self.buffer_vars = ["s_buffer[{0}]".format(mix("+", mix("*", "j", max_process_width), i)) for i in range(max_process_width)]

    code_block.new_line()

    #propagate sum to buffer
    if fold == -1:
      code_block.write_line("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.propagate_into(self.buffer_vars, sum_ptr, "j", 1)
      code_block.dedent()
      code_block.write_line("}")
    else:
      for j in range(fold):
        self.vec.propagate_into(self.s_vars[j], sum_ptr, j, 1)

    code_block.write_line("if(" + " && ".join([inc + " == 1" for inc in self.standard_incs]) + "){")
    code_block.indent()
    self.write_core(code_block, fold, max_process_width, max_unroll, [1 for inc in self.standard_incs])
    code_block.dedent()
    code_block.write_line("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_process_width, max_unroll, self.standard_incs)
    code_block.dedent()
    code_block.write_line("}")
    #consolidate
    if fold == -1:
      code_block.write_line("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.consolidate_into("sum", "j", 1, self.buffer_vars, sum_ptr, "j", 1, self.q_vars[0])
      code_block.dedent()
      code_block.write_line("}")
    else:
      for j in range(fold):
        self.vec.consolidate_into("sum", j, 1, self.s_vars[j], sum_ptr, j, 1, self.q_vars[0])

  def write_core(self, code_block, fold, max_process_width, max_unroll, incs):
    code_block.new_line()
    code_block.write_line("for(i = 0; i + {0} <= n; i += {0}, {1}){{".format(max_unroll, self.data_increment(max_unroll, incs)))
    code_block.indent()
    self.preprocess(code_block, max_unroll, incs)
    self.process(code_block, fold, max_process_width)
    code_block.dedent()
    code_block.write_line("}")
    unroll = max_unroll // 2
    while(unroll >= self.vec.type_size and unroll > 0):
      process_width = self.compute_process_width(unroll)
      code_block.write_line("if(i + {0} <= n){{".format(unroll))
      code_block.indent()
      self.preprocess(code_block, unroll, incs)
      self.process(code_block, fold, process_width)
      code_block.write_line("i += {0}, {1};".format(unroll, self.data_increment(unroll, incs)))
      code_block.dedent()
      code_block.write_line("}")
      unroll //= 2
    if(unroll > 0):
      unroll = self.vec.type_size
      process_width = self.compute_process_width(unroll)
      code_block.write_line("if(i < n){")
      code_block.indent()
      self.preprocess(code_block, unroll, incs, "(n - i)")
      self.process(code_block, fold, process_width)
      code_block.dedent()
      code_block.write_line("}")

  def define_vars(self, code_block, type_name, variables):
    code_block.write_line(type_name + " " + ", ".join(variables) + ";")

  def set_equal(self, code_block, a_vars, b_vars):
    for (a_var, b_var) in zip(a_vars, b_vars):
      code_block.write_line("{0} = {1};".format(a_var, b_var))

  def define_load_vars(self, code_block, process_width):
    raise(NotImplementedError())

  def define_load_ptrs(self, code_block, process_width):
    raise(NotImplementedError())

  def data_increment(self, n, incs):
    return ", ".join(["{0} += {1}".format(load_ptr, mix("*", inc, self.data_type.base_size, n)) for (load_ptr, inc) in zip(self.load_ptrs, incs)])
 # consider here the possible optimization where we multiply incv once by two if the data type is complex.
  
  def preprocess(self, code_block, unroll, incs, partial=""):
    raise(NotImplementedError())

  def process(self, code_block, fold, process_width):
    if(fold == -1):
      code_block.write_line("for(j = 0; j < fold - 1; j++){")
      code_block.indent()
      self.set_equal(code_block, self.s_vars[0], self.buffer_vars[:process_width])
      self.vec.add_BLP_into(self.q_vars, self.s_vars[0], self.load_vars[0], process_width)
      self.set_equal(code_block, self.buffer_vars, self.q_vars[:process_width])
      self.set_equal(code_block, self.q_vars, self.vec.sub(self.s_vars[0], self.q_vars[:process_width]))
      self.set_equal(code_block, self.load_vars[0], self.vec.add(self.load_vars[0], self.q_vars[:process_width]))
      code_block.dedent()
      code_block.write_line("}")
      self.vec.add_BLP_into(self.buffer_vars, self.buffer_vars, self.load_vars[0], process_width)
    else:
      for i in range(process_width // self.vec.suf_width):
        for j in range(fold - 1):
          self.set_equal(code_block, self.q_vars, self.s_vars[j])
          self.vec.add_BLP_into(self.s_vars[j], self.s_vars[j], self.load_vars[0][i * self.vec.suf_width:], self.vec.suf_width)
          self.set_equal(code_block, self.q_vars, self.vec.sub(self.q_vars, self.s_vars[j]))
          self.set_equal(code_block, self.load_vars[0][i * self.vec.suf_width:], self.vec.add(self.load_vars[0][i * self.vec.suf_width:], self.q_vars))      
        self.vec.add_BLP_into(self.s_vars[fold - 1], self.s_vars[fold - 1], self.load_vars[0][i * self.vec.suf_width:], self.vec.suf_width)
    
  def compute_process_width(self, unroll):
    raise(NotImplementedError())

class NonDotOneDimensionalAccumulation(OneDimensionalAccumulation):
  standard_incs = ["incv"]

  def __init__(self, data_type_class, vectorization_class):
    super().__init__(data_type_class, vectorization_class)

  def define_load_vars(self, code_block, process_width):
    self.load_vars = [["v_" + str(i) for i in range(process_width)]]
    self.define_vars(code_block, self.vec.type_name, self.load_vars[0])

  def define_load_ptrs(self, code_block, process_width):
    if self.data_type.is_complex:
      code_block.write_line(self.data_type.base_type.name + "* v_base = (" + self.data_type.base_type.name + "*) v;")
      self.load_ptrs = ["v_base"]
    else:
      self.load_ptrs = ["v"]

  def compute_process_width(self, unroll):
    return (unroll * self.data_type.base_size)//self.vec.base_size

class DotOneDimensionalAccumulation(OneDimensionalAccumulation):
  standard_incs = ["incv", "incy"]

  def __init__(self, data_type_class, vectorization_class):
    super().__init__(data_type_class, vectorization_class)

  def define_load_ptrs(self, code_block, process_width):
    if self.data_type.is_complex:
      code_block.write_line(self.data_type.base_type.name + "* v_base = (" + self.data_type.base_type.name + "*) v;")
      code_block.write_line(self.data_type.base_type.name + "* y_base = (" + self.data_type.base_type.name + "*) y;")
      self.load_ptrs = ["v_base", "y_base"]
    else:
      self.load_ptrs = ["v", "y"]

  def define_load_vars(self, code_block, process_width):
    if self.data_type.is_complex:
      self.load_vars = [["v_" + str(i) for i in range(process_width)], ["y_" + str(i) for i in range(process_width//2)]]
    else:
      self.load_vars = [["v_" + str(i) for i in range(process_width)], ["y_" + str(i) for i in range(process_width)]]
    self.define_vars(code_block, self.vec.type_name, self.load_vars[0])
    self.define_vars(code_block, self.vec.type_name, self.load_vars[1])

  def compute_process_width(self, unroll):
    return (unroll * self.data_type.base_size)//self.vec.base_size * self.data_type.base_size

class SumI(NonDotOneDimensionalAccumulation):
  name = "sumI"

  def __init__(self, data_type_class, vectorization_class):
    super().__init__(data_type_class, vectorization_class)

  @classmethod
  def file_name(cls, data_type_class):
    return "{0}sumI2.c".format(data_type_class.name_char)

  def write_declaration(self, code_block):
    super().write_declaration(code_block)
    code_block.write_line("void {0}sumI2(int n, {1}* v, int incv, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    if partial == "":
      self.set_equal(code_block, self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
    else:
      self.set_equal(code_block, self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))

class ASumI(NonDotOneDimensionalAccumulation):
  name = "asumI"

  def __init__(self, data_type_class, vectorization_class):
    super().__init__(data_type_class, vectorization_class)

  @classmethod
  def file_name(cls, data_type_class):
    redundant_char = ""
    if data_type_class.is_complex:
      redundant_char = data_type_class.base_type.name_char
    return "{0}{1}asumI2.c".format(redundant_char, data_type_class.name_char)

  def write_declaration(self, code_block):
    super().write_declaration(code_block)
    redundant_char = ""
    if self.data_type.is_complex:
      redundant_char = self.data_type.base_type.name_char
    code_block.write_line("void {0}{1}asumI2(int n, {2}* v, int incv, int fold, {2}* sum){{".format(redundant_char, self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    if partial == "":
      self.set_equal(code_block, self.load_vars[0], self.vec.abs(self.vec.load(self.load_ptrs[0], 0, incs[0], unroll)))
    else:
      self.set_equal(code_block, self.load_vars[0], self.vec.abs(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial)))

class Nrm2I(NonDotOneDimensionalAccumulation):
  name = "nrm2I"

  def __init__(self, data_type_class, vectorization_class):
    super().__init__(data_type_class, vectorization_class)

  @classmethod
  def file_name(cls, data_type_class):
    redundant_char = ""
    if data_type_class.is_complex:
      redundant_char = data_type_class.base_type.name_char
    return "{0}{1}nrm2I2.c".format(redundant_char, data_type_class.name_char)

  def write_declaration(self, code_block):
    super().write_declaration(code_block)
    redundant_char = ""
    if self.data_type.is_complex:
      redundant_char = self.data_type.base_type.name_char
    code_block.write_line("void {0}{1}nrm2I2(int n, {2}* v, int incv, {3} scale, int fold, {2}* sum){{".format(redundant_char, self.data_type.name_char, self.data_type.name, self.data_type.base_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    code_block.include_line("{0} scale_mask = {1};".format(self.vec.type_name, self.vec.set("scale")[0]))
    process_width = self.compute_process_width(unroll)
    if partial == "":
      self.set_equal(code_block, self.load_vars[0], self.vec.mul(self.vec.load(self.load_ptrs[0], 0, incs[0], unroll), ["scale_mask"] * process_width))
    else:
      self.set_equal(code_block, self.load_vars[0], self.vec.mul(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial), ["scale_mask"] * process_width))
    self.set_equal(code_block, self.load_vars[0], self.vec.mul(self.load_vars[0], self.load_vars[0][:process_width]))

class DotI(DotOneDimensionalAccumulation):
  def __init__(self, data_type, vectorization):
    assert not data_type.is_complex, "dot is only for real types"
    super().__init__(data_type, vectorization)

  @classmethod
  def file_name(cls, data_type_class):
    return "{0}dotI2.c".format(data_type_class.name_char)

  def write_declaration(self, code_block):
    super().write_declaration(code_block)
    code_block.write_line("void {0}dotI2(int n, {1}* v, int incv, {1}* y, int incy, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def define_preprocess_vars(self):
    return

  def preprocess(self, code_block, unroll, incs, partial=""):
    process_width = self.compute_process_width(unroll)
    if partial == "":
      self.set_equal(code_block, self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
      self.set_equal(code_block, self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], unroll))
    else:
      self.set_equal(code_block, self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      self.set_equal(code_block, self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    
    self.set_equal(code_block, self.load_vars[0], self.vec.mul(self.load_vars[0], self.load_vars[1][:process_width]))

class DotUI(DotOneDimensionalAccumulation):
  def __init__(self, data_type_class, vectorization_class):
    assert data_type_class.is_complex, "dotu is only for complex types"
    super().__init__(data_type_class, vectorization_class)

  @classmethod
  def file_name(cls, data_type_class):
    return "{0}dotuI2.c".format(data_type_class.name_char)

  def write_declaration(self, code_block):
    super().write_declaration(code_block)
    code_block.write_line("void {0}dotuI2(int n, {1}* v, int incv, {1}* y, int incy, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    process_width = self.compute_process_width(unroll)
    if partial == "":
      self.set_equal(code_block, self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
      self.set_equal(code_block, self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], unroll))
    else:
      self.set_equal(code_block, self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      self.set_equal(code_block, self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    self.set_equal(code_block, self.load_vars[0][process_width//2:process_width], self.vec.nconj(self.vec.mul(self.vec.swap_pairwise(self.load_vars[0]), self.vec.rep_odds(self.load_vars[1]))))
    self.set_equal(code_block, self.load_vars[0][:process_width//2], self.vec.mul(self.load_vars[0], self.vec.rep_evens(self.load_vars[1])))

class DotCI(DotOneDimensionalAccumulation):
  def __init__(self, data_type_class, vectorization_class):
    assert data_type_class.is_complex, "dotc is only for complex types"
    super().__init__(data_type_class, vectorization_class)

  @classmethod
  def file_name(self, data_type_class):
    return "{0}dotcI2.c".format(data_type_class.name_char)

  def write_declaration(self, code_block):
    super().write_declaration(code_block)
    code_block.write_line("void {0}dotcI2(int n, {1}* v, int incv, {1}* y, int incy, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    process_width = self.compute_process_width(unroll)
    if partial == "":
      self.set_equal(code_block, self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
      self.set_equal(code_block, self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], unroll))
    else:
      self.set_equal(code_block, self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      self.set_equal(code_block, self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    self.set_equal(code_block, self.load_vars[0][process_width//2:process_width], self.vec.conj(self.vec.mul(self.vec.swap_pairwise(self.load_vars[0]), self.vec.rep_odds(self.load_vars[1]))))
    self.set_equal(code_block, self.load_vars[0][:process_width//2], self.vec.mul(self.load_vars[0], self.vec.rep_evens(self.load_vars[1])))

class CodeBlock:
  def __init__(self, srcFile, base_indent_level = 0):
    self.base_indent_level = base_indent_level
    self.indent_level = base_indent_level
    self.blocks = []
    self.srcFile = srcFile
    self.included = set()
    self.includes = []

  def indent(self):
    self.indent_level += 1

  def dedent(self):
    assert self.indent_level > self.base_indent_level, "attempting to dedent beyond writable area"
    self.indent_level -= 1

  def write_line(self, line):
    self.blocks += ["  " * self.indent_level + line]

  def include_line(self, line):
    if line not in self.included:
      self.includes += ["  " * self.base_indent_level + line]
      self.included.add(line)

  def new_line(self):
    self.blocks += [""]

  def sub_block(self):
    block = CodeBlock(self.srcFile, base_indent_level = self.indent_level)
    self.blocks += [block]
    return block

  def __str__(self):
    return "\n".join([str(block) for block in (self.includes + self.blocks)])

class SrcFile:
  def __init__(self, name, prelude = []):
    self.name = name
    self.code = CodeBlock(self)
    self.prelude = self.code.sub_block()
    for line in prelude:
      self.prelude.write_line(line)

  def write_line(self, line):
    self.code.write_line(line)

  def include_line(self, line):
    self.code.include_line(line)

  def sub_block(self):
    return self.code.sub_block()

  def write(self):
    f = open(self.name, 'w')
    f.write(str(self.code))
    f.close()
  
def IBLASFile(target_function, data_type_class, implementations):
  output_file = SrcFile(target_function.file_name(data_type_class))
  output_block = output_file.sub_block()

  output_block.new_line()

  for (i, (vectorization, unroll_settings)) in enumerate(implementations):
    func = target_function(data_type_class, vectorization)
    if len(implementations) > 1:
      if i == 0:
        output_block.write_line("#if defined( " + vectorization.defined_macro + " )")
      elif i < len(implementations) - 1:
        output_block.write_line("#elif defined( " + vectorization.defined_macro + " )")
      else:
        output_block.write_line("#else")
      output_block.indent()
    func.set_settings(unroll_settings)
    func.write(output_block)
    if len(implementations) > 1:
      output_block.dedent()
  if len(implementations) > 1:
    output_block.write_line("#endif")

  output_file.write()
  return output_block


vec = 0
        
IBLASFile(SumI, Float, [(AVX, [(3, 32), (-1, 8)]), (SSE, [(3, 8), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(SumI, Double, [(AVX, [(3, 16), (-1, 8)]), (SSE, [(3, 8), (-1, 4)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(SumI, FloatComplex, [(AVX, [(3, 16), (-1, 8)]), (SSE, [(3, 4), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(SumI, DoubleComplex, [(AVX, [(3, 8), (-1, 8)]), (SSE, [(3, 4), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(ASumI, Float, [(AVX, [(3, 32), (-1, 8)]), (SSE, [(3, 8), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(ASumI, Double, [(AVX, [(3, 16), (-1, 8)]), (SSE, [(3, 4), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(ASumI, FloatComplex, [(AVX, [(3, 16), (-1, 8)]), (SSE, [(3, 4), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(ASumI, DoubleComplex, [(AVX, [(3, 8), (-1, 8)]), (SSE, [(3, 2), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(Nrm2I, Float, [(AVX, [(3, 32), (-1, 8)]), (SSE, [(3, 8), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(Nrm2I, Double, [(AVX, [(3, 16), (-1, 8)]), (SSE, [(3, 4), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(Nrm2I, FloatComplex, [(AVX, [(3, 16), (-1, 8)]), (SSE, [(3, 4), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(Nrm2I, DoubleComplex, [(AVX, [(3, 8), (-1, 8)]), (SSE, [(3, 2), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotI, Float, [(AVX, [(3, 32), (-1, 8)]), (SSE, [(3, 8), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotI, Double, [(AVX, [(3, 16), (-1, 8)]), (SSE, [(3, 4), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotUI, FloatComplex, [(AVX, [(3, 16), (-1, 8)]), (SSE, [(3, 4), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotUI, DoubleComplex, [(AVX, [(3, 8), (-1, 8)]), (SSE, [(3, 2), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotCI, FloatComplex, [(AVX, [(3, 16), (-1, 8)]), (SSE, [(3, 4), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotCI, DoubleComplex, [(AVX, [(3, 8), (-1, 8)]), (SSE, [(3, 2), (-1, 8)]), (SISD, [(3,     2), (-1, 2)])][vec:]) 
