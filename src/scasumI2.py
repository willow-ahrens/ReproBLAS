import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import asumI2

generate.generate(asumI2.ASumI2(dataTypes.FloatComplex), __file__, generate.get_tuning(__file__))
