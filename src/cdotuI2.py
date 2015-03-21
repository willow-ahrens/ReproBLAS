import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import dotuI2

generate.generate(dotuI2.DotUI2(dataTypes.FloatComplex), __file__, generate.get_tuning(__file__))
