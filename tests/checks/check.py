import os

import tests.checks.checks as checks
import tests.harness.harness as harness

check_dir = os.path.dirname(os.path.abspath(__file__))

check_suite = checks.CheckSuite()

#folds = [1, 2, 3, 4]
folds = [3]
inf_folds = [1, 3, 4]

check_suite.add_checks([checks.ValidateInternalDSCALETest(),\
                        checks.ValidateInternalDSCALETest()],\
                       ["N", "incX"],\
                       [[4], [1]])

"""
check_suite.add_checks([checks.ValidateInternalUFPTest(),\
                        checks.ValidateInternalUFPFTest()],\
                       ["N", "incX"],\
                       [[10], [1, 2, 4]])

check_suite.add_checks([checks.ValidateInternalDINDEXTest(),\
                        checks.ValidateInternalSINDEXTest(),\
                        checks.ValidateInternalDMINDEXTest(),\
                        checks.ValidateInternalSMINDEXTest()],\
                       ["N", "incX"],\
                       [[4], [1]])

check_suite.add_checks([checks.ValidateInternalDAMAXTest(),\
                        checks.ValidateInternalZAMAXTest(),\
                        checks.ValidateInternalSAMAXTest(),\
                        checks.ValidateInternalCAMAXTest()],\
                       ["N", "incX"],\
                       [[4095], [1, 2, 4]])

check_suite.add_checks([checks.ValidateInternalRDSUMTest(),\
                        checks.ValidateInternalDIDIADDTest(),\
                        checks.ValidateInternalDIDADDTest(),\
                        checks.ValidateInternalDIDDEPOSITTest(),\
                        checks.ValidateInternalRSSUMTest(),\
                        checks.ValidateInternalSISIADDTest(),\
                        checks.ValidateInternalSISADDTest(),\
                        checks.ValidateInternalSISDEPOSITTest(),\
                        ],\
                       ["N", "fold", "incX", "RealScaleX", "f"],\
                       [[4095], folds, [1, 4], [1.0, -1.0],\
                        ["constant",\
                         "+big",\
                         "++big",\
                         "+-big",\
                         "sine"]])

check_suite.add_checks([checks.ValidateInternalRZSUMTest(),\
                        checks.ValidateInternalZIZIADDTest(),\
                        checks.ValidateInternalZIZADDTest(),\
                        checks.ValidateInternalZIZDEPOSITTest(),\
                        checks.ValidateInternalRCSUMTest(),\
                        checks.ValidateInternalCICIADDTest(),\
                        checks.ValidateInternalCICADDTest(),\
                        checks.ValidateInternalCICDEPOSITTest()],\
                       ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "f"],\
                       [[4095], folds, [1, 4], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0],\
                        ["constant",\
                         "+big",\
                         "++big",\
                         "+-big",\
                         "sine"]])

check_suite.add_checks([checks.ValidateInternalRDNRM2Test(),\
                        checks.ValidateInternalRDASUMTest(),\
                        checks.ValidateInternalRSNRM2Test(),\
                        checks.ValidateInternalRSASUMTest(),\
                        ],\
                       ["N", "fold", "incX", "RealScaleX", "f"],\
                       [[4095], folds, [1, 4], [1.0, -1.0],\
                        ["constant",\
                         "+big",\
                         "++big",\
                         "+-big"]])

check_suite.add_checks([checks.ValidateInternalRDZNRM2Test(),\
                        checks.ValidateInternalRDZASUMTest(),\
                        checks.ValidateInternalRSCNRM2Test(),\
                        checks.ValidateInternalRSCASUMTest(),\
                        ],\
                       ["N", "fold", "incX", ("RealScaleX", "ImagScaleX"), "f"],\
                       [[4095], folds, [1, 4], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0],\
                        ["constant",\
                         "+big",\
                         "++big",\
                         "+-big"]])

check_suite.add_checks([checks.ValidateInternalRDDOTTest(),\
                        checks.ValidateInternalRSDOTTest(),\
                        ],\
                       ["N", "fold", "incX", "RealScaleX", "RealScaleY", "f", "g"],\
                       [[4095], folds, [1, 4], [1.0, -1.0], [1.0, -1.0],\
                        ["constant",\
                         "+big",\
                         "++big",\
                         "+-big"],\
                        ["constant",\
                         "+big",\
                         "++big",\
                         "+-big"]])

check_suite.add_checks([checks.ValidateInternalRZDOTUTest(),\
                        checks.ValidateInternalRZDOTCTest(),\
                        checks.ValidateInternalRCDOTUTest(),\
                        checks.ValidateInternalRCDOTCTest(),\
                        ],\
                       ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "RealScaleY", "ImagScaleY", "f", "g"],\
                       [[4095], folds, [1, 4], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0],\
                        ["constant",\
                         "+big",\
                         "++big",\
                         "+-big"],\
                        ["constant",\
                         "+big",\
                         "++big",\
                         "+-big"]])

check_suite.add_checks([checks.ValidateInternalRZDOTUTest(),\
                        checks.ValidateInternalRZDOTCTest(),\
                        checks.ValidateInternalRCDOTUTest(),\
                        checks.ValidateInternalRCDOTCTest(),\
                        ],\
                       ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "RealScaleY", "ImagScaleY", ("f", "g")],\
                       [[4095], folds, [1, 4], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0],\
                        [("constant", "sine"),\
                         ("sine", "constant")]])

check_suite.add_checks([checks.ValidateInternalRDSUMTest(),\
                        checks.ValidateInternalRDASUMTest(),\
                        checks.ValidateInternalRDNRM2Test(),\
                        checks.ValidateInternalDIDIADDTest(),\
                        checks.ValidateInternalDIDADDTest(),\
                        checks.ValidateInternalDIDDEPOSITTest(),\
                        checks.ValidateInternalRSSUMTest(),\
                        checks.ValidateInternalRSASUMTest(),\
                        checks.ValidateInternalRSNRM2Test(),\
                        checks.ValidateInternalSISIADDTest(),\
                        checks.ValidateInternalSISADDTest(),\
                        checks.ValidateInternalSISDEPOSITTest(),\
                        ],\
                       ["N", "fold", "incX", "RealScaleX", "f"],\
                       [[4095], inf_folds, [1, 4], [1.0, -1.0],\
                        ["+inf",\
                         "++inf",\
                         "+-inf",\
                         "nan",\
                         "+inf_nan",\
                         "++inf_nan",\
                         "+-inf_nan"]])

check_suite.add_checks([checks.ValidateInternalRZSUMTest(),\
                        checks.ValidateInternalRDZASUMTest(),\
                        checks.ValidateInternalRDZNRM2Test(),\
                        checks.ValidateInternalZIZIADDTest(),\
                        checks.ValidateInternalZIZADDTest(),\
                        checks.ValidateInternalZIZDEPOSITTest(),\
                        checks.ValidateInternalRCSUMTest(),\
                        checks.ValidateInternalRSCASUMTest(),\
                        checks.ValidateInternalRSCNRM2Test(),\
                        checks.ValidateInternalCICIADDTest(),\
                        checks.ValidateInternalCICADDTest(),\
                        checks.ValidateInternalCICDEPOSITTest()],\
                       ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "f"],\
                       [[4095], inf_folds, [1, 4], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0],\
                        ["+inf",\
                         "++inf",\
                         "+-inf",\
                         "nan",\
                         "+inf_nan",\
                         "++inf_nan",\
                         "+-inf_nan"]])

check_suite.add_checks([checks.ValidateInternalRDDOTTest(),\
                        checks.ValidateInternalRSDOTTest(),\
                        ],\
                       ["N", "fold", "incX", "RealScaleX", "RealScaleY", "f", "g"],\
                       [[4095], inf_folds, [1, 4], [1.0, -1.0], [1.0, -1.0],\
                        ["constant",\
                         "+inf",\
                         "++inf",\
                         "+-inf",\
                         "nan",\
                         "+inf_nan",\
                         "++inf_nan",\
                         "+-inf_nan"],\
                        ["constant",\
                         "+inf",\
                         "++inf",\
                         "+-inf",\
                         "nan",\
                         "+inf_nan",\
                         "++inf_nan",\
                         "+-inf_nan"]])

check_suite.add_checks([checks.ValidateInternalRZDOTUTest(),\
                        checks.ValidateInternalRZDOTCTest(),\
                        checks.ValidateInternalRCDOTUTest(),\
                        checks.ValidateInternalRCDOTCTest(),\
                        ],\
                       ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "RealScaleY", "ImagScaleY", "f", "g"],\
                       [[4095], inf_folds, [1, 4], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0],
                        ["constant",\
                         "+inf",\
                         "++inf",\
                         "+-inf",\
                         "nan",\
                         "+inf_nan",\
                         "++inf_nan",\
                         "+-inf_nan"],\
                        ["constant",\
                         "+inf",\
                         "++inf",\
                         "+-inf",\
                         "nan",\
                         "+inf_nan",\
                         "++inf_nan",\
                         "+-inf_nan"]])

check_suite.add_checks([checks.VerifyRDSUMTest(),\
                        checks.VerifyRDASUMTest(),\
                        checks.VerifyDIDIADDTest(),\
                        checks.VerifyDIDADDTest(),\
                        checks.VerifyDIDDEPOSITTest(),\
                        checks.VerifyRZSUMTest(),\
                        checks.VerifyRDZASUMTest(),\
                        checks.VerifyZIZIADDTest(),\
                        checks.VerifyZIZADDTest(),\
                        checks.VerifyZIZDEPOSITTest(),\
                        checks.VerifyRSSUMTest(),\
                        checks.VerifyRSASUMTest(),\
                        checks.VerifySISIADDTest(),\
                        checks.VerifySISADDTest(),\
                        checks.VerifySISDEPOSITTest(),\
                        checks.VerifyRCSUMTest(),\
                        checks.VerifyRSCASUMTest(),\
                        checks.VerifyCICIADDTest(),\
                        checks.VerifyCICADDTest(),\
                        checks.VerifyCICDEPOSITTest()],\
                       ["N", "fold", "B", "incX", "RealScaleX", "f"],\
                       [[4095], folds, [256], [1, 4], [0],\
                        ["rand"]])

check_suite.add_checks([checks.VerifyRDSUMTest(),\
                        checks.VerifyRDASUMTest(),\
                        checks.VerifyRDNRM2Test(),\
                        checks.VerifyDIDIADDTest(),\
                        checks.VerifyDIDADDTest(),\
                        checks.VerifyDIDDEPOSITTest(),\
                        checks.VerifyRZSUMTest(),\
                        checks.VerifyRDZASUMTest(),\
                        checks.VerifyRDZNRM2Test(),\
                        checks.VerifyZIZIADDTest(),\
                        checks.VerifyZIZADDTest(),\
                        checks.VerifyZIZDEPOSITTest(),\
                        checks.VerifyRSSUMTest(),\
                        checks.VerifyRSASUMTest(),\
                        checks.VerifyRSNRM2Test(),\
                        checks.VerifySISIADDTest(),\
                        checks.VerifySISADDTest(),\
                        checks.VerifySISDEPOSITTest(),\
                        checks.VerifyRCSUMTest(),\
                        checks.VerifyRSCASUMTest(),\
                        checks.VerifyRSCNRM2Test(),\
                        checks.VerifyCICIADDTest(),\
                        checks.VerifyCICADDTest(),\
                        checks.VerifyCICDEPOSITTest()],\
                       ["N", "fold", "B", "incX", "f"],\
                       [[4095], folds, [256], [1, 4],\
                        ["rand",\
                         "rand+(rand-1)",\
                         "sine",\
                         "small+grow*big"]])

check_suite.add_checks([checks.VerifyRDDOTTest(),\
                        checks.VerifyRZDOTUTest(),\
                        checks.VerifyRZDOTCTest(),\
                        checks.VerifyRSDOTTest(),\
                        checks.VerifyRCDOTUTest(),\
                        checks.VerifyRCDOTCTest()],\
                       ["N", "fold", "incX", "incY", "f", "g"],\
                       [[4095], folds, [1, 4], [1, 4],\
                        ["rand",\
                         "rand+(rand-1)",\
                         "sine",\
                         "small+grow*big"],\
                        ["rand",\
                         "rand+(rand-1)",\
                         "sine",\
                         "small+grow*big"]])

DBL_BIN_WIDTH=40
DBL_MAX_EXP=1024
DBL_BIG_EXP=27
DBL_SMALL_EXP=-27
DBL_MIN_EXP=-1021

for i in range(DBL_BIN_WIDTH + 2):
  check_suite.add_checks([checks.ValidateInternalRDSUMTest(),\
                          checks.ValidateInternalDIDIADDTest(),\
                          checks.ValidateInternalDIDADDTest(),\
                          checks.ValidateInternalDIDDEPOSITTest(),\
                          checks.ValidateInternalRDASUMTest(),\
                          #checks.ValidateInternalRDNRM2Test(),\
                          checks.ValidateInternalRDDOTTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "f", "g"],\
                         [[32], folds, [1, 4],\
                         [1.5 * 2**(DBL_MAX_EXP - DBL_BIG_EXP - 6 - i), 0.75 * 2**(DBL_MIN_EXP - DBL_SMALL_EXP + i)],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"],\
                          ["constant"]])

  check_suite.add_checks([checks.ValidateInternalRDDOTTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "f", "g"],\
                         [[32], folds, [1, 4],\
                         [1.5 * 2**(DBL_MAX_EXP - DBL_BIG_EXP - 6 - i), 0.75 * 2**(DBL_MIN_EXP - DBL_SMALL_EXP + i)],\
                          ["constant"],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"]])

  check_suite.add_checks([checks.ValidateInternalRDDOTTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "f", "g"],\
                         [[32], folds, [1, 4],\
                         [1.5 * 2**(DBL_MAX_EXP - 2 * DBL_BIG_EXP - 6 - i), 0.75 * 2**(DBL_MIN_EXP - 2 * DBL_SMALL_EXP + i)],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"]])

  check_suite.add_checks([checks.ValidateInternalRZSUMTest(),\
                          checks.ValidateInternalZIZIADDTest(),\
                          checks.ValidateInternalZIZADDTest(),\
                          checks.ValidateInternalZIZDEPOSITTest(),\
                          checks.ValidateInternalRDZASUMTest(),\
                          #checks.ValidateInternalRDZNRM2Test(),\
                          checks.ValidateInternalRZDOTUTest(),\
                          checks.ValidateInternalRZDOTCTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "f", "g"],\
                         [[16], folds, [1, 4],\
                          [1.5 * 2**(DBL_MAX_EXP - DBL_BIG_EXP - 6 - i), 0.75 * 2**(DBL_MIN_EXP - DBL_SMALL_EXP + i)],\
                          [1.5 * 2**(DBL_MAX_EXP - DBL_BIG_EXP - 6 - i), 0.75 * 2**(DBL_MIN_EXP - DBL_SMALL_EXP + i)],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"],\
                          ["constant"]])

  check_suite.add_checks([checks.ValidateInternalRZDOTUTest(),\
                          checks.ValidateInternalRZDOTCTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "f", "g"],\
                         [[16], folds, [1, 4],\
                          [1.5 * 2**(DBL_MAX_EXP - DBL_BIG_EXP - 6 - i), 0.75 * 2**(DBL_MIN_EXP - DBL_SMALL_EXP + i)],\
                          [1.5 * 2**(DBL_MAX_EXP - DBL_BIG_EXP - 6 - i), 0.75 * 2**(DBL_MIN_EXP - DBL_SMALL_EXP + i)],\
                          ["constant"],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"]])

  check_suite.add_checks([checks.ValidateInternalRZDOTUTest(),\
                          checks.ValidateInternalRZDOTCTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "f", "g"],\
                         [[16], folds, [1, 4],\
                          [1.5 * 2**(DBL_MAX_EXP - 2 * DBL_BIG_EXP - 6 - i), 0.75 * 2**(DBL_MIN_EXP - 2 * DBL_SMALL_EXP + i)],\
                          [1.5 * 2**(DBL_MAX_EXP - 2 * DBL_BIG_EXP - 6 - i), 0.75 * 2**(DBL_MIN_EXP - 2 * DBL_SMALL_EXP + i)],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"]])



FLT_BIN_WIDTH=13
FLT_MAX_EXP=128
FLT_BIG_EXP=13
FLT_SMALL_EXP=-12
FLT_MIN_EXP=-125

for i in range(FLT_BIN_WIDTH + 2):
  check_suite.add_checks([checks.ValidateInternalRSSUMTest(),\
                          checks.ValidateInternalSISIADDTest(),\
                          checks.ValidateInternalSISADDTest(),\
                          checks.ValidateInternalSISDEPOSITTest(),\
                          checks.ValidateInternalRSASUMTest(),\
                          #checks.ValidateInternalRSNRM2Test(),\
                          checks.ValidateInternalRSDOTTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "f", "g"],\
                         [[32], folds, [1, 4],\
                         [1.5 * 2**(FLT_MAX_EXP - FLT_BIG_EXP - 6 - i), 0.75 * 2**(FLT_MIN_EXP - FLT_SMALL_EXP + i)],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"],\
                          ["constant"]])

  check_suite.add_checks([checks.ValidateInternalRSDOTTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "f", "g"],\
                         [[32], folds, [1, 4],\
                         [1.5 * 2**(FLT_MAX_EXP - FLT_BIG_EXP - 6 - i), 0.75 * 2**(FLT_MIN_EXP - FLT_SMALL_EXP + i)],\
                          ["constant"],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"]])

  check_suite.add_checks([checks.ValidateInternalRSDOTTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "f", "g"],\
                         [[32], folds, [1, 4],\
                         [1.5 * 2**(FLT_MAX_EXP - 2 * FLT_BIG_EXP - 6 - i), 0.75 * 2**(FLT_MIN_EXP - 2 * FLT_SMALL_EXP + i)],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"]])

  check_suite.add_checks([checks.ValidateInternalRCSUMTest(),\
                          checks.ValidateInternalCICIADDTest(),\
                          checks.ValidateInternalCICADDTest(),\
                          checks.ValidateInternalCICDEPOSITTest(),\
                          checks.ValidateInternalRSCASUMTest(),\
                          #checks.ValidateInternalRSCNRM2Test(),\
                          checks.ValidateInternalRCDOTUTest(),\
                          checks.ValidateInternalRCDOTCTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "f", "g"],\
                         [[16], folds, [1, 4],\
                          [1.5 * 2**(FLT_MAX_EXP - FLT_BIG_EXP - 6 - i), 0.75 * 2**(FLT_MIN_EXP - FLT_SMALL_EXP + i)],\
                          [1.5 * 2**(FLT_MAX_EXP - FLT_BIG_EXP - 6 - i), 0.75 * 2**(FLT_MIN_EXP - FLT_SMALL_EXP + i)],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"],\
                          ["constant"]])

  check_suite.add_checks([checks.ValidateInternalRCDOTUTest(),\
                          checks.ValidateInternalRCDOTCTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "f", "g"],\
                         [[16], folds, [1, 4],\
                          [1.5 * 2**(FLT_MAX_EXP - FLT_BIG_EXP - 6 - i), 0.75 * 2**(FLT_MIN_EXP - FLT_SMALL_EXP + i)],\
                          [1.5 * 2**(FLT_MAX_EXP - FLT_BIG_EXP - 6 - i), 0.75 * 2**(FLT_MIN_EXP - FLT_SMALL_EXP + i)],\
                          ["constant"],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"]])

  check_suite.add_checks([checks.ValidateInternalRCDOTUTest(),\
                          checks.ValidateInternalRCDOTCTest(),\
                          ],\
                         ["N", "fold", "incX", "RealScaleX", "ImagScaleX", "f", "g"],\
                         [[16], folds, [1, 4],\
                          [1.5 * 2**(FLT_MAX_EXP - 2 * FLT_BIG_EXP - 6 - i), 0.75 * 2**(FLT_MIN_EXP - 2 * FLT_SMALL_EXP + i)],\
                          [1.5 * 2**(FLT_MAX_EXP - 2 * FLT_BIG_EXP - 6 - i), 0.75 * 2**(FLT_MIN_EXP - 2 * FLT_SMALL_EXP + i)],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"],\
                          ["constant",\
                           "+big",\
                           "++big",\
                           "+-big"]])
"""

#check_suite.add_checks([checks.VerifyRDGEMVTest()],\
#                       ["O", "T", "N", "M", "lda", "incX", "incY", "f", "g", "j"],\
#                       [["RowMajor", "ColMajor"], ["NoTrans", "Trans"], [1023], [1023], [1023, 1025], [1, 4], [1, 4],\
#                        ["rand",\
#                         "small+grow*big"],\
#                        ["rand",\
#                         "small+grow*big"],\
#                        ["rand",\
#                         "small+grow*big"]])

check_harness = harness.Harness("check")
check_harness.add_suite(check_suite)
check_harness.run()
