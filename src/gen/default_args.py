import argparse
from generate import *

parser = argparse.ArgumentParser(description="Create a default args file from params file")
parser.add_argument("-p", "--params", type=str, required=True, help="params file")
parser.add_argument("-a", "--args", type=str, required=True, help="args file")
args = parser.parse_args()

parameter_space = deserialize_parameter_space(args.params)

serialize_arguments(parameter_space.get_default_arguments(), args.args)
