################################################################################
# vectorizations.py                                                            #
#                                                                              #
#     A set of classes used to generate generically vectorized code. Currently #
# supports SISD, Intel SSE, Intel AVX.                                         #
#                                                                              #
#                                                            Peter Ahrens 2014 #
################################################################################

import math
from utils import *
from dataTypes import *

class Vectorization(object):
  name = ""
  defined_macro = ""

  def __init__(self, code_block, data_type_class):
    self.code_block = code_block
    self.data_type = data_type_class(code_block)

  def define_unroll_step(self, fold, unroll):
    code_block.write("#define " + self.unroll_step_string(fold) + " " + str(unroll))

  def define_max_unroll_step(self, fold, max_unroll_width):
    code_block.write("#define MAX_" + self.unroll_step_string(fold) + " " + str(max_unroll_width))

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
  def load(self, src_ptr, offset, inc, n, align=False):
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

  def max(self, src1_vars, src2_vars):
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

  def iterate_unrolled(self, i_var, n_var, src_ptrs, src_incs, max_unroll, min_unroll, body):
    i = 0
    unroll = (max_unroll // max(self.type_size, 1)) * max(self.type_size, 1)
    while(unroll >= min_unroll):
      if i == 0:
        self.code_block.write("for({0} = 0; {0} + {1} <= {2}; {0} += {1}, {3}){{".format(i_var, max_unroll, n_var, self.data_type.data_increment(src_ptrs, src_incs, max_unroll)))
        self.code_block.indent()
        body(unroll)
        self.code_block.dedent()
        self.code_block.write("}")
      elif unroll < self.type_size:
        self.code_block.write("if({0} < {1}){{".format(i_var, n_var))
        self.code_block.indent()
        body("({0} - {1})".format(n_var, i_var))
        self.code_block.dedent()
        self.code_block.write("}")
        break
      else:
        self.code_block.write("if({0} + {1} <= {2}){{".format(i_var, unroll, n_var))
        self.code_block.indent()
        body(unroll)
        self.code_block.write("{0} += {1}, {2};".format(i_var, unroll, self.data_type.data_increment(src_ptrs, src_incs, unroll)))
        self.code_block.dedent()
        self.code_block.write("}")
      if math.log(unroll, 2) % 1 != 0:
        unroll = 2**int(math.log(unroll, 2))
      else:
        unroll //= 2
      i += 1;

"""
  ##TODO This is currently broken on wierd unroll sizes. fix that.
  def iterate_unrolled_aligned(self, i_var, n_var, src_ptrs, src_incs, max_unroll, min_unroll, body):
    align = False
    if self.type_size > 1 and src_incs[0] == 1:
      self.code_block.set_equal(i_var, ["((uintptr_t){0} & {1}) >> {2}".format(src_ptrs[0], self.byte_size - 1, int(math.floor(math.log(self.data_type.byte_size, 2))))])
      self.code_block.write("if({0} != 0 && {0} < {1}){{".format(i_var, n_var))
      self.code_block.indent()
#      self.code_block.write("printf(\"v = %p\\n\", {0});".format(src_ptrs[0]))
      body("{0}".format(i_var))
      self.code_block.write("{0};".format(self.data_type.data_increment(src_ptrs, src_incs, i_var)))
#      self.code_block.write("printf(\"v = %p, i = %d\\n\", {0}, {1});".format(src_ptrs[0], i_var))
      self.code_block.dedent()
      self.code_block.write("}else{")
      self.code_block.indent()
      self.code_block.set_equal(i_var, ["0"])
      self.code_block.dedent()
      self.code_block.write("}")
      self.code_block.write("for(; {0} + {1} <= {2}; {0} += {1}, {3}){{".format(i_var, max_unroll, n_var, self.data_type.data_increment(src_ptrs, src_incs, max_unroll)))
      align = True
    else:
      self.code_block.write("for({0} = 0; {0} + {1} <= {2}; {0} += {1}, {3}){{".format(i_var, max_unroll, n_var, self.data_type.data_increment(src_ptrs, src_incs, max_unroll)))
    self.code_block.indent()
    body(max_unroll, align)
    self.code_block.dedent()
    self.code_block.write("}")
    unroll = max_unroll // 2
    while(unroll >= self.type_size and unroll >= min_unroll):
      self.code_block.write("if({0} + {1} <= {2}){{".format(i_var, unroll, n_var))
      self.code_block.indent()
      body(unroll, align)
      self.code_block.write("{0} += {1}, {2};".format(i_var, unroll, self.data_type.data_increment(src_ptrs, src_incs, unroll)))
      self.code_block.dedent()
      self.code_block.write("}")
      unroll //=2
    if(unroll >= min_unroll):
      self.code_block.write("if({0} < {1}){{".format(i_var, n_var))
      self.code_block.indent()
      body("({0} - {1})".format(n_var, i_var), align)
      self.code_block.dedent()
      self.code_block.write("}")
"""


class SISD(Vectorization):
  name = "SISD"

  def __init__(self, code_block, data_type_class):
    super(SISD, self).__init__(code_block, data_type_class)
    self.bit_size = self.data_type.base_type.bit_size
    self.byte_size = self.data_type.base_type.byte_size
    self.type_name = self.data_type.base_type.name
    self.base_size = self.bit_size//self.data_type.base_type.bit_size
    self.type_size = self.bit_size//self.data_type.bit_size #notice that type_size is 0 in SISD vectorization if the data type is complex. It is probably better not to use type_size for SISD purposes.
    self.zero = 0
    self.suf_width = self.data_type.base_size

  def consolidate_into(self, dst_ptr, offset, inc, src_vars, common_summand_ptr, common_summand_offset, common_summand_inc):
    self.code_block.include("{} cons_tmp;".format(self.type_name))

    if self.data_type.is_complex:
      dst_ptr = "(({0}*){1})".format(self.data_type.base_type.name, dst_ptr)
      common_summand_ptr = "(({0}*){1})".format(self.data_type.base_type.name, common_summand_ptr)

    if self.data_type.is_complex:
      if(len(src_vars) > 2):
        self.code_block.write("cons_tmp = {0}[{1}];".format(common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
      for src_var in src_vars[2::2]:
        self.code_block.write("{0} = {0} + ({1} - cons_tmp);".format(src_vars[0], src_var))
      if(len(src_vars) > 2):
        self.code_block.write("cons_tmp = {0}[{1}];".format(common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 1)))
      for src_var in src_vars[3::2]:
        self.code_block.write("{0} = {0} + ({1} - cons_tmp);".format(src_vars[1], src_var))
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0), src_vars[0]))
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 1), src_vars[1]))
    else:
      if(len(src_vars) > 1):
        self.code_block.write("cons_tmp = {0}[{1}];".format(common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
      for src_var in src_vars[1:]:
        self.code_block.write("{0} = {0} + ({1} - cons_tmp);".format(src_vars[0], src_var))
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0), src_vars[0]))

  def max_into(self, dst_ptr, offset, inc, src_vars):
    if self.data_type.is_complex:
      dst_ptr = "(({0}*){1})".format(self.data_type.base_type.name, dst_ptr)

    for i in range(self.data_type.base_size, self.base_size):
      self.code_block.set_equal(src_vars[i % self.data_type.base_size], self.max(src_vars[i % self.data_type.base_size], src_vars[i]))
    if self.data_type.is_complex:
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0), src_vars[0]))
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 1), src_vars[1]))
    else:
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0), src_vars[0]))

  def propagate_into(self, dst_vars, src_ptr, offset, inc):
    if self.data_type.is_complex:
      src_ptr = "(({0}*){1})".format(self.data_type.base_type.name, src_ptr)

    if self.data_type.is_complex:
      assert len(dst_vars) % 2 == 0, "cannot propagate complex value to odd number of base type dst_vars"
      self.code_block.write(" = ".join(dst_vars[0::2]) + " = {0}[{1}];".format(src_ptr, self.data_type.index(offset, inc, 0)))
      self.code_block.write(" = ".join(dst_vars[1::2]) + " = {0}[{1}];".format(src_ptr, self.data_type.index(offset, inc, 1)))
    else:
      self.code_block.write(" = ".join(dst_vars) + " = {0}[{1}];".format(src_ptr, self.data_type.index(offset, inc, 0)))

  def include_BLP_vars(self):
    self.code_block.include("{0} tmp_BLP;".format(self.data_type.base_type.int_float_name))

  def add_BLP_into(self, dst, src, blp, width):
    assert len(dst) >= width
    assert len(src) >= width
    assert len(blp) >= width

    self.include_BLP_vars()
    for i in range(width):
      self.code_block.write("tmp_BLP.{0} = {1};".format(self.data_type.base_type.float_char, blp[i]))
      self.code_block.write("tmp_BLP.{0} |= 1;".format(self.data_type.base_type.int_char))
      self.code_block.write("{0} = {1} + tmp_BLP.{2};".format(dst[i], src[i], self.data_type.base_type.float_char))

  def load(self, src_ptr, offset, inc, n, align=False):
    assert n > 0, "n must be nonzero"
    assert inc != 0, "inc must be nonzero"

    if self.data_type.is_complex:
      src_ptr = "(({0}*){1})".format(self.data_type.base_type.name, src_ptr)

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

  def max(self, src1_vars, src2_vars):
    return ["({0} > {1}? {0}: {1})".format(src1_var, src2_var) for (src1_var, src2_var) in zip(src1_vars, src2_vars)]

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
    super(SIMD, self).__init__(code_block, data_type_class)
    self.suf_width = 1

  def include_max_vars(self):
    self.code_block.include("{0} max_tmp[{1}] __attribute__((aligned({2})));".format(self.data_type.base_type.name, self.base_size, self.byte_size))

  def include_consolidation_vars(self):
    self.code_block.include("{0} cons_tmp;".format(self.type_name))
    self.code_block.include("{0} cons_buffer[{1}] __attribute__((aligned({2})));".format(self.data_type.base_type.name, self.base_size, self.byte_size))

  #TODO masks should be defined here and not necessarily in COMMON?
  def include_ABS_vars(self):
    self.code_block.include("{0} abs_mask; {1}_ABS_MASK{2}(abs_mask);".format(self.type_name, self.name, self.data_type.base_type.name_char.upper()))

  def include_BLP_vars(self):
    self.code_block.include("{0} blp_mask; {1}_BLP_MASK{2}(blp_mask);".format(self.type_name, self.name, self.data_type.base_type.name_char.upper()))

  def include_CONJ_vars(self):
    self.code_block.include("{0} conj_mask; {1}_CONJ_MASK{2}(conj_mask);".format(self.type_name, self.name, self.data_type.base_type.name_char.upper())) 

  def include_NCONJ_vars(self):
    self.code_block.include("{0} nconj_mask; {1}_NCONJ_MASK{2}(nconj_mask);".format(self.type_name, self.name, self.data_type.base_type.name_char.upper()))


class SSE(SIMD):
  name = "SSE"
  defined_macro = "__SSE2__"

  def __init__(self, code_block, data_type_class):
    super(SSE, self).__init__(code_block, data_type_class)
    self.bit_size = 128
    self.byte_size = 16
    self.type_name = {"float": "__m128", "double": "__m128d"}[self.data_type.base_type.name]
    self.base_size = self.bit_size//self.data_type.base_type.bit_size
    self.type_size = self.bit_size//self.data_type.bit_size
    self.zero = "_mm_setzero_p{0}()".format(self.data_type.base_type.name_char)

  def consolidate_into(self, dst_ptr, offset, inc, src_vars, common_summand_ptr, common_summand_offset, common_summand_inc):
    self.include_consolidation_vars()

    if self.data_type.is_complex:
      dst_ptr = "(({0}*){1})".format(self.data_type.base_type.name, dst_ptr)
      common_summand_ptr = "(({0}*){1})".format(self.data_type.base_type.name, common_summand_ptr)

    if self.data_type.name == "float":
      self.code_block.write("{0} = _mm_sub_ps({0}, _mm_set_ps({1}[{2}], {1}[{2}], {1}[{2}], 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
    elif self.data_type.name == "double":
      self.code_block.write("{0} = _mm_sub_pd({0}, _mm_set_pd({1}[{2}], 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
    elif self.data_type.name == "float complex":
      self.code_block.write("{0} = _mm_sub_ps({0}, _mm_set_ps({1}[{3}], {1}[{2}], 0, 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0), self.data_type.index(common_summand_offset, common_summand_inc, 1)))
    if len(src_vars) > 1:
      self.propagate_into(["cons_tmp"], common_summand_ptr, common_summand_offset, common_summand_inc)
      for src_var in src_vars[1:]:
        self.code_block.write("{0} = _mm_add_p{1}({0}, _mm_sub_p{1}({2}, cons_tmp));".format(src_vars[0], self.data_type.base_type.name_char, src_var))
    self.code_block.write("_mm_store_p{0}(cons_buffer, {1});".format(self.data_type.base_type.name_char, src_vars[0]))
    if self.data_type.is_complex:
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0), " + ".join(["cons_buffer[{0}]".format(self.data_type.index(i, 1, 0)) for i in range(self.type_size)])))
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 1), " + ".join(["cons_buffer[{0}]".format(self.data_type.index(i, 1, 1)) for i in range(self.type_size)])))
    else:
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0), " + ".join(["cons_buffer[{0}]".format(i) for i in range(self.base_size)])))

  def max_into(self, dst_ptr, offset, inc, src_vars):
    self.include_max_vars()

    if self.data_type.is_complex:
      dst_ptr = "(({0}*){1})".format(self.data_type.base_type.name, dst_ptr)

    for src_var in src_vars[1:]:
      self.set_equal(self.max(src_var[0], src_var));
    self.code_block.write("_mm_store_p{0}(max_tmp, {1});".format(self.data_type.base_type.name_char, src_vars[0]))
    for i in range(self.data_type.base_size, self.base_size):
      self.code_block.write("max_tmp[{0}] = (max_tmp[{0}] > max_tmp[{1}] ? max_tmp[{0}]: max_tmp[{1}]);".format(i % self.data_type.base_size, i))
    for i in range(self.data_type.base_size):
      self.code_block.write("{0}[{1}] = max_tmp[{2}];".format(dst_ptr, self.data_type.index(offset, inc, i, True), i))

  def propagate_into(self, dst_vars, src_ptr, offset, inc):
    if self.data_type.is_complex:
      src_ptr = "(({0}*){1})".format(self.data_type.base_type.name, src_ptr)

    if self.data_type.is_complex:
      broadcast = {"float": "(__m128)_mm_load1_pd((double *)({0}));", "double": "_mm_loadu_pd({0});"}[self.data_type.base_type.name]
    else:
      broadcast = "_mm_load1_p{0}({{0}});".format(self.data_type.name_char)
    self.code_block.write(" = ".join(dst_vars) + " = " + broadcast.format(mix("+", src_ptr, self.data_type.index(offset, inc, 0), paren=False)))

  def add_BLP_into(self, dst, src, blp, width):
    assert len(dst) >= width
    assert len(src) >= width
    assert len(blp) >= width

    self.include_BLP_vars()
    for i in range(width):
      self.code_block.write("{0} = _mm_add_p{1}({2}, _mm_or_p{1}({3}, blp_mask));".format(dst[i], self.data_type.base_type.name_char, src[i], blp[i]))

  def load(self, src_ptr, offset, inc, n, align=False):
    assert n > 0, "n must be nonzero"
    assert n % self.type_size == 0, "n must be a multiple of the number of types that fit in a vector"

    if self.data_type.is_complex:
      src_ptr = "(({0}*){1})".format(self.data_type.base_type.name, src_ptr)

    result = []
    for i in range(n//self.type_size):
      if inc == 1 or self.type_size == 1:
        if align:
          result += ["_mm_load_p{0}({1})".format(self.data_type.base_type.name_char, mix("+", src_ptr, self.data_type.index(offset, inc, i * self.base_size), paren=False))]
        else:
          result += ["_mm_loadu_p{0}({1})".format(self.data_type.base_type.name_char, mix("+", src_ptr, self.data_type.index(offset, inc, i * self.base_size), paren=False))]
      else:
        if self.data_type.base_type.name == "float":
          result += ["_mm_set_ps({0}[{4}], {0}[{3}], {0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, i * self.base_size), self.data_type.index(offset, inc, i * self.base_size + 1), self.data_type.index(offset, inc, i * self.base_size + 2), self.data_type.index(offset, inc, i * self.base_size + 3))]
        if self.data_type.base_type.name == "double":
          result += ["_mm_set_pd({0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, i * self.base_size), self.data_type.index(offset, inc, i * self.base_size + 1))]
    return result

  def load_partial(self, src_ptr, offset, inc, n):
    if self.data_type.is_complex:
      src_ptr = "(({0}*){1})".format(self.data_type.base_type.name, src_ptr)

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
    return ["_mm_and_p{0}({1}, abs_mask)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]

  def mul(self, src_vars, amt_vars):
    return ["_mm_mul_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def max(self, src1_vars, src2_vars):
    return ["_mm_max_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src1_var, src2_var) for (src1_var, src2_var) in zip(src1_vars, src2_vars)]

  def conj(self, src_vars):
    if self.data_type.is_complex:
      self.include_CONJ_vars()
      return ["_mm_xor_p{0}({1}, conj_mask)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]
    else:
      return src_vars

  def nconj(self, src_vars):
    if self.data_type.is_complex:
      self.include_NCONJ_vars()
      return ["_mm_xor_p{0}({1}, nconj_mask)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]
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
    super(AVX, self).__init__(code_block, data_type_class)
    self.bit_size = 256
    self.byte_size = 32 
    self.type_name = {"float": "__m256", "double": "__m256d"}[self.data_type.base_type.name]
    self.base_size = self.bit_size//self.data_type.base_type.bit_size
    self.type_size = self.bit_size//self.data_type.bit_size
    self.zero = "_mm256_setzero_p{0}()".format(self.data_type.base_type.name_char)


  def consolidate_into(self, dst_ptr, offset, inc, src_vars, common_summand_ptr, common_summand_offset, common_summand_inc):
    self.include_consolidation_vars()

    if self.data_type.is_complex:
      dst_ptr = "(({0}*){1})".format(self.data_type.base_type.name, dst_ptr)
      common_summand_ptr = "(({0}*){1})".format(self.data_type.base_type.name, common_summand_ptr)

    if self.data_type.name == "float":
      self.code_block.write("{0} = _mm256_sub_ps({0}, _mm256_set_ps({1}[{2}], {1}[{2}], {1}[{2}], {1}[{2}], {1}[{2}], {1}[{2}], {1}[{2}], 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
    elif self.data_type.name == "double":
      self.code_block.write("{0} = _mm256_sub_pd({0}, _mm256_set_pd({1}[{2}], {1}[{2}], {1}[{2}], 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0)))
    elif self.data_type.name == "float complex":
      self.code_block.write("{0} = _mm256_sub_ps({0}, _mm256_set_ps({1}[{3}], {1}[{2}], {1}[{3}], {1}[{2}], {1}[{3}], {1}[{2}], 0, 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0), self.data_type.index(common_summand_offset, common_summand_inc, 1)))
    elif self.data_type.name == "double complex":
      self.code_block.write("{0} = _mm256_sub_pd({0}, _mm256_set_pd({1}[{3}], {1}[{2}], 0, 0));".format(src_vars[0], common_summand_ptr, self.data_type.index(common_summand_offset, common_summand_inc, 0), self.data_type.index(common_summand_offset, common_summand_inc, 1)))
    if len(src_vars) > 1:
      self.propagate_into(["cons_tmp"], common_summand_ptr, common_summand_offset, common_summand_inc)
      for src_var in src_vars[1:]:
        self.code_block.write("{0} = _mm256_add_p{1}({0}, _mm256_sub_p{1}({2}, cons_tmp));".format(src_vars[0], self.data_type.base_type.name_char, src_var))
    self.code_block.write("_mm256_store_p{0}(cons_buffer, {1});".format(self.data_type.base_type.name_char, src_vars[0]))
    if self.data_type.is_complex:
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0), " + ".join(["cons_buffer[{0}]".format(2 * i) for i in range(self.type_size)])))
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 1), " + ".join(["cons_buffer[{0}]".format(2 * i + 1) for i in range(self.type_size)])))
    else:
      self.code_block.write("{0}[{1}] = {2};".format(dst_ptr, self.data_type.index(offset, inc, 0), " + ".join(["cons_buffer[{0}]".format(i) for i in range(self.base_size)])))

  def max_into(self, dst_ptr, offset, inc, src_vars): 
    self.include_max_vars()

    if self.data_type.is_complex:
      dst_ptr = "(({0}*){1})".format(self.data_type.base_type.name, dst_ptr)

    for src_var in src_vars[1:]:
      self.set_equal(self.max(src_var[0], src_var));
    self.code_block.write("_mm256_store_p{0}(max_tmp, {1});".format(self.data_type.base_type.name_char, src_vars[0]))
    for i in range(self.data_type.base_size, self.base_size):
      self.code_block.write("max_tmp[{0}] = (max_tmp[{0}] > max_tmp[{1}] ? max_tmp[{0}]: max_tmp[{1}]);".format(i % self.data_type.base_size, i))
    for i in range(self.data_type.base_size):
      self.code_block.write("{0}[{1}] = max_tmp[{2}];".format(dst_ptr, self.data_type.index(offset, inc, i, True), i))

  def propagate_into(self, dst_vars, src_ptr, offset, inc):
    if self.data_type.is_complex:
      src_ptr = "(({0}*){1})".format(self.data_type.base_type.name, src_ptr)

    if self.data_type.is_complex:
      broadcast = {"float": "(__m256)_mm256_broadcast_sd((double *)({0}));", "double": "_mm256_broadcast_pd((__m128d *)({0}));"}[self.data_type.base_type.name]
    else:
      broadcast = "_mm256_broadcast_s{0}({{0}});".format(self.data_type.name_char)
    self.code_block.write(" = ".join(dst_vars) + " = " + broadcast.format(mix("+", src_ptr, self.data_type.index(offset, inc, 0), paren=False)))

  def add_BLP_into(self, dst, src, blp, width):
    assert len(dst) >= width
    assert len(src) >= width
    assert len(blp) >= width
    self.include_BLP_vars()
    for i in range(width):
      self.code_block.write("{0} = _mm256_add_p{1}({2}, _mm256_or_p{1}({3}, blp_mask));".format(dst[i], self.data_type.base_type.name_char, src[i], blp[i]))

  def load(self, src_ptr, offset, inc, n, align=False):
    assert n > 0, "n must be nonzero"
    assert n % self.type_size == 0, "n must be a multiple of the number of types that fit in a vector"

    if self.data_type.is_complex:
      src_ptr = "(({0}*){1})".format(self.data_type.base_type.name, src_ptr)

    result = []
    for i in range(n//self.type_size):
      if inc == 1 or self.type_size == 1:
        if align:
          result += ["_mm256_load_p{0}({1})".format(self.data_type.base_type.name_char, mix("+", src_ptr, self.data_type.index(offset, inc, i * self.base_size), paren=False))]
        else:
          result += ["_mm256_loadu_p{0}({1})".format(self.data_type.base_type.name_char, mix("+", src_ptr, self.data_type.index(offset, inc, i * self.base_size), paren=False))]
      else:
        if self.data_type.base_type.name == "float":
          result += ["_mm256_set_ps({0}[{8}], {0}[{7}], {0}[{6}], {0}[{5}], {0}[{4}], {0}[{3}], {0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, i * self.base_size), self.data_type.index(offset, inc, i * self.base_size + 1), self.data_type.index(offset, inc, i * self.base_size + 2), self.data_type.index(offset, inc, i * self.base_size + 3), self.data_type.index(offset, inc, i * self.base_size + 4), self.data_type.index(offset, inc, i * self.base_size + 5), self.data_type.index(offset, inc, i * self.base_size + 6), self.data_type.index(offset, inc, i * self.base_size + 7))]
        elif self.data_type.base_type.name == "double":
          result += ["_mm256_set_pd({0}[{4}], {0}[{3}], {0}[{2}], {0}[{1}])".format(src_ptr, self.data_type.index(offset, inc, i * self.base_size), self.data_type.index(offset, inc, i * self.base_size + 1), self.data_type.index(offset, inc, i * self.base_size + 2), self.data_type.index(offset, inc, i * self.base_size + 3))]
    return result

  def load_partial(self, src_ptr, offset, inc, n):
    if self.data_type.is_complex:
      src_ptr = "(({0}*){1})".format(self.data_type.base_type.name, src_ptr)

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
    return ["_mm256_and_p{0}({1}, abs_mask)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]

  def mul(self, src_vars, amt_vars):
    return ["_mm256_mul_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src_var, amt_var) for (src_var, amt_var) in zip(src_vars, amt_vars)]

  def max(self, src1_vars, src2_vars):
    return ["_mm256_max_p{0}({1}, {2})".format(self.data_type.base_type.name_char, src1_var, src2_var) for (src1_var, src2_var) in zip(src1_vars, src2_vars)]

  def conj(self, src_vars):
    if self.data_type.is_complex:
      self.include_CONJ_vars()
      return ["_mm256_xor_p{0}({1}, conj_mask)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]
    else:
      return src_vars

  def nconj(self, src_vars):
    if self.data_type.is_complex:
      self.include_NCONJ_vars()
      return ["_mm256_xor_p{0}({1}, nconj_mask)".format(self.data_type.base_type.name_char, src_var) for src_var in src_vars]
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

vectorization_lookup = {"SISD":SISD, "SSE":SSE, "AVX":AVX}

#all_vectorizations = [AVX, SSE, SISD]
all_vectorizations = [AVX, SSE, SISD]

def iterate_all_vectorizations(f, code_block):
  for (i, vectorization) in enumerate(all_vectorizations):
    if i == 0 and len(all_vectorizations) > 1:
      code_block.write("#ifdef {}".format(vectorization.defined_macro))
    elif i < len(all_vectorizations) - 1:
      code_block.write("#elif defined({})".format(vectorization.defined_macro))
    elif len(all_vectorizations) > 1:
      code_block.write("#else")
    code_block.indent()
    f(vectorization, code_block.sub_block())
    code_block.dedent()
  if len(all_vectorizations) > 1:
    code_block.write("#endif")
