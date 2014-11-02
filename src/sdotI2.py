import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import dotI2

generate.generate(dotI2.DotI2(dataTypes.Float), __file__, generate.get_tuning(__file__))
