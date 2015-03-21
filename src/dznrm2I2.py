import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import nrm2I2

generate.generate(nrm2I2.Nrm2I2(dataTypes.DoubleComplex), __file__, generate.get_tuning(__file__))
