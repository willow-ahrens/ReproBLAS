import os

import tests.checks.checks as checks
import tests.harness.harness as harness

check_dir = os.path.dirname(os.path.abspath(__file__))

check_suite = checks.CheckSuite()

check_suite.add_checks([checks.ValidateInternalRCSUMTest(),\
                        checks.ValidateInternalCICIADDTest(),\
                        checks.ValidateInternalCICADDTest(),\
                        checks.ValidateInternalCICDEPOSITTest(),\
                        ],\
                       ["N", "incX", ("RealScaleX", "ImagScaleX"), "f"],\
                       [[4095], [1, 4], [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],\
                        ["constant"]])

check_suite.add_checks([checks.ValidateInternalRSCNRM2Test(),\
                        ],\
                       ["N", "incX", ("RealScaleX", "ImagScaleX"), "f"],\
                       [[4095], [1, 4], [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],\
                        ["constant"]])

check_suite.add_checks([checks.ValidateInternalRCDOTUTest(),\
                        checks.ValidateInternalRCDOTCTest(),\
                        ],\
                       ["N", "incX", ("RealScaleX", "ImagScaleX"), ("RealScaleY", "ImagScaleY"), "f", "g"],\
                       [[4095], [1, 4], [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],\
                        ["constant"],\
                        ["constant"]])

check_harness = harness.Harness("check")
check_harness.add_suite(check_suite)
check_harness.run()
