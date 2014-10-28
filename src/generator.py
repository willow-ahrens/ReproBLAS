import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
from utils import *
from dataTypes import *
from vectorizations import *


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
    code_block.write("}")

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
    code_block.srcFile.include("#include <stdio.h>")
    code_block.srcFile.include("#include <stdlib.h>")
    code_block.srcFile.include("#include <float.h>")
    code_block.srcFile.include("#include <math.h>")
    code_block.srcFile.include("#include \"config.h\"")
    code_block.srcFile.include("#include \"Common/Common.h\"")
    #code_block.srcFile.include("#include \"IndexedFP/" + self.data_type.base_type.name_char + "Indexed.h\"")
    #code_block.srcFile.include("#include \"rblas1.h\"")

  def write_body(self, code_block):
    code_block.write("SET_DAZ_FLAG;")
    if len(self.settings) == 1:
      self.write_cores(-1, self.settings[0][1])
      code_block.write("break;")
    else:
      code_block.write("switch(fold){")
      code_block.indent()
      for (fold, max_unroll_width) in self.settings:
        if fold == -1:
          code_block.write("default:{")
        else:
          code_block.write("case " + str(fold) + ":{")
        code_block.indent()
        self.write_cores(code_block, fold, max_unroll_width)
        #code_block.write("break;")
        code_block.write("RESET_DAZ_FLAG")
        code_block.write("return;")
        code_block.dedent()
        code_block.write("}")
      code_block.dedent()
      code_block.write("}")

  def write_cores(self, code_block, fold, max_unroll):
#    if self.vec.name == "AVX":
#      self.code_block.write('printf("Hi im avx {0}\\n");'.format(self.name))
#    elif self.vec.name == "SSE":
#      self.code_block.write('printf("Hi im sse {0}\\n");'.format(self.name))
#    else:
#      self.code_block.write('printf("Hi im sisd {0}\\n");'.format(self.name))

# max_unroll_width is the maximum number of elements of the input vector that are processed in each iteration.
# unroll_width is the number of elements of the input vector that are processed in each iteration.
# process_width is the number of variables that are processed in each iteration.
    max_process_width = self.compute_process_width(max_unroll)
    if fold == -1:
      code_block.write("int i, j;")
    else:
      code_block.write("int i;")
    code_block.new_line()
    sum_ptr = "sum"
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* sum_base = (" + self.data_type.base_type.name + "*) sum;")
      sum_ptr = "sum_base"
    self.define_load_ptrs(code_block, max_process_width)
    self.define_load_vars(code_block, max_process_width)
    if fold == -1:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(max_process_width)]
      code_block.define_vars(self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_" + str(i) for i in range(max_process_width)]]
      code_block.define_vars(self.vec.type_name, self.s_vars[0])
    else:
      #define q variables
      self.q_vars = ["q_" + str(i) for i in range(self.vec.suf_width)]
      code_block.define_vars(self.vec.type_name, self.q_vars)
      #define s variables
      self.s_vars = [["s_{0}_{1}".format(j, i) for i in range(self.vec.suf_width)] for j in range(fold)]
      for j in range(fold):
        code_block.define_vars(self.vec.type_name, self.s_vars[j])

    if fold == -1:
      code_block.write("{0} s_buffer[{1}];".format(self.vec.type_name, mix("*", max_process_width, "MAX_FOLD")))
      self.buffer_vars = ["s_buffer[{0}]".format(mix("+", mix("*", "j", max_process_width), i)) for i in range(max_process_width)]

    code_block.new_line()

    #propagate sum to buffer
    if fold == -1:
      code_block.write("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.propagate_into(self.buffer_vars, sum_ptr, "j", 1)
      code_block.dedent()
      code_block.write("}")
    else:
      for j in range(fold):
        self.vec.propagate_into(self.s_vars[j], sum_ptr, j, 1)

    code_block.write("if(" + " && ".join([inc + " == 1" for inc in self.standard_incs]) + "){")
    code_block.indent()
    self.write_core(code_block, fold, max_process_width, max_unroll, [1 for inc in self.standard_incs])
    code_block.dedent()
    code_block.write("}else{")
    code_block.indent()
    self.write_core(code_block, fold, max_process_width, max_unroll, self.standard_incs)
    code_block.dedent()
    code_block.write("}")
    #consolidate
    if fold == -1:
      code_block.write("for(j = 0; j < fold; j += 1){")
      code_block.indent()
      self.vec.consolidate_into("sum", "j", 1, self.buffer_vars, sum_ptr, "j", 1, self.q_vars[0])
      code_block.dedent()
      code_block.write("}")
    else:
      for j in range(fold):
        self.vec.consolidate_into("sum", j, 1, self.s_vars[j], sum_ptr, j, 1, self.q_vars[0])

  def write_core(self, code_block, fold, max_process_width, max_unroll, incs):
    code_block.new_line()
    def body(unroll):
      if type(unroll) == str:
        process_width = self.compute_process_width(self.vec.type_size)
        self.preprocess(code_block, self.vec.type_size, incs, unroll)
        self.process(code_block, fold, process_width)
      else:
        process_width = self.compute_process_width(unroll)
        self.preprocess(code_block, unroll, incs)
        self.process(code_block, fold, process_width)
    self.vec.iterate_unrolled("i", "n", self.load_ptrs, incs, max_unroll, 1, body)

  def define_load_vars(self, code_block, process_width):
    raise(NotImplementedError())

  def define_load_ptrs(self, code_block, process_width):
    raise(NotImplementedError())

  def preprocess(self, code_block, unroll, incs, partial=""):
    raise(NotImplementedError())

  def process(self, code_block, fold, process_width):
    if(fold == -1):
      code_block.write("for(j = 0; j < fold - 1; j++){")
      code_block.indent()
      code_block.set_equal(self.s_vars[0], self.buffer_vars[:process_width])
      self.vec.add_BLP_into(self.q_vars, self.s_vars[0], self.load_vars[0], process_width)
      code_block.set_equal(self.buffer_vars, self.q_vars[:process_width])
      code_block.set_equal(self.q_vars, self.vec.sub(self.s_vars[0], self.q_vars[:process_width]))
      code_block.set_equal(self.load_vars[0], self.vec.add(self.load_vars[0], self.q_vars[:process_width]))
      code_block.dedent()
      code_block.write("}")
      self.vec.add_BLP_into(self.buffer_vars, self.buffer_vars, self.load_vars[0], process_width)
    else:
      for i in range(process_width // self.vec.suf_width):
        for j in range(fold - 1):
          code_block.set_equal(self.q_vars, self.s_vars[j])
          self.vec.add_BLP_into(self.s_vars[j], self.s_vars[j], self.load_vars[0][i * self.vec.suf_width:], self.vec.suf_width)
          code_block.set_equal(self.q_vars, self.vec.sub(self.q_vars, self.s_vars[j]))
          code_block.set_equal(self.load_vars[0][i * self.vec.suf_width:], self.vec.add(self.load_vars[0][i * self.vec.suf_width:], self.q_vars))      
        self.vec.add_BLP_into(self.s_vars[fold - 1], self.s_vars[fold - 1], self.load_vars[0][i * self.vec.suf_width:], self.vec.suf_width)
    
  def compute_process_width(self, unroll):
    raise(NotImplementedError())

class NonDotOneDimensionalAccumulation(OneDimensionalAccumulation):
  standard_incs = ["incv"]

  def __init__(self, data_type_class, vectorization_class):
    super().__init__(data_type_class, vectorization_class)

  def define_load_vars(self, code_block, process_width):
    self.load_vars = [["v_" + str(i) for i in range(process_width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])

  def define_load_ptrs(self, code_block, process_width):
    if self.data_type.is_complex:
      code_block.write(self.data_type.base_type.name + "* v_base = (" + self.data_type.base_type.name + "*) v;")
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
      code_block.write(self.data_type.base_type.name + "* v_base = (" + self.data_type.base_type.name + "*) v;")
      code_block.write(self.data_type.base_type.name + "* y_base = (" + self.data_type.base_type.name + "*) y;")
      self.load_ptrs = ["v_base", "y_base"]
    else:
      self.load_ptrs = ["v", "y"]

  def define_load_vars(self, code_block, process_width):
    if self.data_type.is_complex:
      self.load_vars = [["v_" + str(i) for i in range(process_width)], ["y_" + str(i) for i in range(process_width//2)]]
    else:
      self.load_vars = [["v_" + str(i) for i in range(process_width)], ["y_" + str(i) for i in range(process_width)]]
    code_block.define_vars(self.vec.type_name, self.load_vars[0])
    code_block.define_vars(self.vec.type_name, self.load_vars[1])

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
    code_block.write("void {0}sumI2(int n, {1}* v, int incv, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))

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
    code_block.write("void {0}{1}asumI2(int n, {2}* v, int incv, int fold, {2}* sum){{".format(redundant_char, self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load(self.load_ptrs[0], 0, incs[0], unroll)))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.abs(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial)))

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
    code_block.write("void {0}{1}nrm2I2(int n, {2}* v, int incv, {3} scale, int fold, {2}* sum){{".format(redundant_char, self.data_type.name_char, self.data_type.name, self.data_type.base_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    code_block.include("{0} scale_mask = {1};".format(self.vec.type_name, self.vec.set("scale")[0]))
    process_width = self.compute_process_width(unroll)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.mul(self.vec.load(self.load_ptrs[0], 0, incs[0], unroll), ["scale_mask"] * process_width))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.mul(self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial), ["scale_mask"] * process_width))
    code_block.set_equal(self.load_vars[0], self.vec.mul(self.load_vars[0], self.load_vars[0][:process_width]))

class DotI(DotOneDimensionalAccumulation):
  def __init__(self, data_type, vectorization):
    assert not data_type.is_complex, "dot is only for real types"
    super().__init__(data_type, vectorization)

  @classmethod
  def file_name(cls, data_type_class):
    return "{0}dotI2.c".format(data_type_class.name_char)

  def write_declaration(self, code_block):
    super().write_declaration(code_block)
    code_block.write("void {0}dotI2(int n, {1}* v, int incv, {1}* y, int incy, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def define_preprocess_vars(self):
    return

  def preprocess(self, code_block, unroll, incs, partial=""):
    process_width = self.compute_process_width(unroll)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], unroll))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    
    code_block.set_equal(self.load_vars[0], self.vec.mul(self.load_vars[0], self.load_vars[1][:process_width]))

class DotUI(DotOneDimensionalAccumulation):
  def __init__(self, data_type_class, vectorization_class):
    assert data_type_class.is_complex, "dotu is only for complex types"
    super().__init__(data_type_class, vectorization_class)

  @classmethod
  def file_name(cls, data_type_class):
    return "{0}dotuI2.c".format(data_type_class.name_char)

  def write_declaration(self, code_block):
    super().write_declaration(code_block)
    code_block.write("void {0}dotuI2(int n, {1}* v, int incv, {1}* y, int incy, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    process_width = self.compute_process_width(unroll)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], unroll))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    code_block.set_equal(self.load_vars[0][process_width//2:process_width], self.vec.nconj(self.vec.mul(self.vec.swap_pairwise(self.load_vars[0]), self.vec.rep_odds(self.load_vars[1]))))
    code_block.set_equal(self.load_vars[0][:process_width//2], self.vec.mul(self.load_vars[0], self.vec.rep_evens(self.load_vars[1])))

class DotCI(DotOneDimensionalAccumulation):
  def __init__(self, data_type_class, vectorization_class):
    assert data_type_class.is_complex, "dotc is only for complex types"
    super().__init__(data_type_class, vectorization_class)

  @classmethod
  def file_name(self, data_type_class):
    return "{0}dotcI2.c".format(data_type_class.name_char)

  def write_declaration(self, code_block):
    super().write_declaration(code_block)
    code_block.write("void {0}dotcI2(int n, {1}* v, int incv, {1}* y, int incy, int fold, {1}* sum){{".format(self.data_type.name_char, self.data_type.name))

  def preprocess(self, code_block, unroll, incs, partial=""):
    process_width = self.compute_process_width(unroll)
    if partial == "":
      code_block.set_equal(self.load_vars[0], self.vec.load(self.load_ptrs[0], 0, incs[0], unroll))
      code_block.set_equal(self.load_vars[1], self.vec.load(self.load_ptrs[1], 0, incs[1], unroll))
    else:
      code_block.set_equal(self.load_vars[0], self.vec.load_partial(self.load_ptrs[0], 0, incs[0], partial))
      code_block.set_equal(self.load_vars[1], self.vec.load_partial(self.load_ptrs[1], 0, incs[1], partial))
    code_block.set_equal(self.load_vars[0][process_width//2:process_width], self.vec.conj(self.vec.mul(self.vec.swap_pairwise(self.load_vars[0]), self.vec.rep_odds(self.load_vars[1]))))
    code_block.set_equal(self.load_vars[0][:process_width//2], self.vec.mul(self.load_vars[0], self.vec.rep_evens(self.load_vars[1])))
  
def IBLASFile(target_function, data_type_class, implementations):
  output_file = SrcFile(target_function.file_name(data_type_class))
  output_block = output_file.sub_block()

  output_block.new_line()

  for (i, (vectorization, unroll_settings)) in enumerate(implementations):
    func = target_function(data_type_class, vectorization)
    if len(implementations) > 1:
      if i == 0:
        output_block.write("#if defined( " + vectorization.defined_macro + " )")
      elif i < len(implementations) - 1:
        output_block.write("#elif defined( " + vectorization.defined_macro + " )")
      else:
        output_block.write("#else")
      output_block.indent()
    func.set_settings(unroll_settings)
    func.write(output_block)
    if len(implementations) > 1:
      output_block.dedent()
  if len(implementations) > 1:
    output_block.write("#endif")

  output_file.dump()
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
