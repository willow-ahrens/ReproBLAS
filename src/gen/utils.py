################################################################################
# utils.py                                                                     #
#                                                                              #
#     Helps with the actual writing (indentation, formatting, logisitcs) of    #
# generated code.                                                              #
#                                                                              #
#                                                            Peter Ahrens 2014 #
################################################################################

def mix(op, *args, **kwargs):
  paren = True
  if "paren" in kwargs:
    paren = kwargs["paren"]
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

class CodeBlock(object):
  def __init__(self, base_indent_level = 0):
    self.base_indent_level = base_indent_level
    self.indent_level = base_indent_level
    self.blocks = []
    self.included = set()
    self.includes = []

  def indent(self):
    self.indent_level += 1

  def dedent(self):
    assert self.indent_level > self.base_indent_level, "attempting to dedent beyond writable area"
    self.indent_level -= 1

  def write(self, chunk):
    self.blocks += ["\n".join(["  " * self.indent_level + line for line in chunk.split("\n")])]

  def include(self, chunk):
    if chunk not in self.included:
      self.includes += ["\n".join(["  " * self.base_indent_level + line for line in chunk.split("\n")])]
      self.included.add(chunk)

  def define_vars(self, type_name, variables):
    self.write(type_name + " " + ", ".join(variables) + ";")

  def set_equal(self, a_vars, b_vars):
    for (a_var, b_var) in zip(a_vars, b_vars):
      self.write("{0} = {1};".format(a_var, b_var))

  def new_line(self):
    self.blocks += [""]

  def sub_block(self):
    block = CodeBlock(base_indent_level = self.indent_level)
    self.blocks += [block]
    return block

  def __str__(self):
    return "\n".join([str(block) for block in (self.includes + self.blocks)])
