
class Nrm2I(NonDotOneDimensionalAccumulation):
  name = "nrm2I"

  def __init__(self, data_type_class, vectorization_class):
    super(Nrm2I, self).__init__(data_type_class, vectorization_class)

  def write_declaration(self, code_block):
    super(Nrm2I, self).write_declaration(code_block)
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
    super(DotI, self).__init__(data_type, vectorization)

  def write_declaration(self, code_block):
    super(DotI, self).write_declaration(code_block)
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
    super(DotUI, self).__init__(data_type_class, vectorization_class)

  def write_declaration(self, code_block):
    super(DotUI, self).write_declaration(code_block)
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
    super(DotCI, self).__init__(data_type_class, vectorization_class)

  def write_declaration(self, code_block):
    super(DotCI, self).write_declaration(code_block)
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
  def tgt(code_block, vectorization, settings):
    tgt = target_function(data_type_class, vectorization)
    tgt.set_settings(settings)
    tgt.write(code_block)
  generate(, __file__, generate.get_settings())


vec = 0
        
IBLASFile(SumI, Float, [("AVX", [(3, 32), (-1, 8)]), ("SSE", [(3, 8), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(SumI, Double, [("AVX", [(3, 16), (-1, 8)]), ("SSE", [(3, 8), (-1, 4)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(SumI, FloatComplex, [("AVX", [(3, 16), (-1, 8)]), ("SSE", [(3, 4), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(SumI, DoubleComplex, [("AVX", [(3, 8), (-1, 8)]), ("SSE", [(3, 4), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(ASumI, Float, [("AVX", [(3, 32), (-1, 8)]), ("SSE", [(3, 8), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(ASumI, Double, [("AVX", [(3, 16), (-1, 8)]), ("SSE", [(3, 4), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(ASumI, FloatComplex, [("AVX", [(3, 16), (-1, 8)]), ("SSE", [(3, 4), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(ASumI, DoubleComplex, [("AVX", [(3, 8), (-1, 8)]), ("SSE", [(3, 2), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(Nrm2I, Float, [("AVX", [(3, 32), (-1, 8)]), ("SSE", [(3, 8), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(Nrm2I, Double, [("AVX", [(3, 16), (-1, 8)]), ("SSE", [(3, 4), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(Nrm2I, FloatComplex, [("AVX", [(3, 16), (-1, 8)]), ("SSE", [(3, 4), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(Nrm2I, DoubleComplex, [("AVX", [(3, 8), (-1, 8)]), ("SSE", [(3, 2), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotI, Float, [("AVX", [(3, 32), (-1, 8)]), ("SSE", [(3, 8), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotI, Double, [("AVX", [(3, 16), (-1, 8)]), ("SSE", [(3, 4), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotUI, FloatComplex, [("AVX", [(3, 16), (-1, 8)]), ("SSE", [(3, 4), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotUI, DoubleComplex, [("AVX", [(3, 8), (-1, 8)]), ("SSE", [(3, 2), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotCI, FloatComplex, [("AVX", [(3, 16), (-1, 8)]), ("SSE", [(3, 4), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
IBLASFile(DotCI, DoubleComplex, [("AVX", [(3, 8), (-1, 8)]), ("SSE", [(3, 2), (-1, 8)]), ("SISD", [(3,     2), (-1, 2)])][vec:]) 
