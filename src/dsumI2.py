import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "gen"))
import generate
import dataTypes
import sumI2

generate.generate(sumI2.SumI2(dataTypes.Double), __file__, generate.get_tuning(__file__))
