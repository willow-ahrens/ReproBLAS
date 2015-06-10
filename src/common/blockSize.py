import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
from utils import *
from dataTypes import *
from vectorizations import *
from generate import *
from scripts import terminal
import itertools

class BlockSize(Target):
  def __init__(self, name, var_name, minimum, maximum, default, metrics):
    super(BlockSize, self).__init__()
    self.name = name
    self.var_name = var_name
    self.minimum = minimum
    self.maximum = maximum
    self.default = default
    self.metrics = metrics

  def get_arguments(self):
    arguments = []
    arguments.append("{}_block_size_{}".format(self.name, self.var_name))
    return arguments

  def get_metrics(self):
    return {"{}_block_size_{}".format(self.name, self.var_name): self.metrics}

  def get_parameters(self):
    parameters = []
    parameters.append(PowerOfTwoParameter("{}_block_size_{}".format(self.name, self.var_name), {}, self.minimum, self.maximum, self.default))
    return parameters

  def write(self, code_block):
    code_block.write("#define {} {}".format(self.var_name, self.arguments["{}_block_size_{}".format(self.name, self.var_name)]))
