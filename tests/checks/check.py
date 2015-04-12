import os

import tests.checks.checks as checks
import tests.harness.harness as harness

check_dir = os.path.dirname(os.path.abspath(__file__))

check_suite = checks.CheckSuite()

check_suite.add_checks([checks.ValidateInternalUFPTest(),\
                        checks.ValidateInternalUFPFTest()],\
                       ["N", "incX"],\
                       [[10], [1, 2, 4]])

check_suite.add_checks([checks.ValidateInternalDAMAXTest(),\
                        checks.ValidateInternalZAMAXTest(),\
                        checks.ValidateInternalSAMAXTest(),\
                        checks.ValidateInternalCAMAXTest()],\
                       ["N", "incX"],\
                       [[4095], [1, 2, 4]])

check_suite.add_checks([checks.ValidateInternalRDBLAS1Test(),\
                        checks.ValidateInternalRZBLAS1Test(),\
                        checks.ValidateInternalRSBLAS1Test(),\
                        checks.ValidateInternalRCBLAS1Test()],\
                       ["N", "incX", "incY"],\
                       [[4095], [1, 4], [1, 4]])

check_suite.add_checks([checks.VerifyRDSUMTest(),\
                        checks.VerifyRDASUMTest(),\
                        checks.VerifyRDNRM2Test(),\
                        checks.VerifyRZSUMTest(),\
                        checks.VerifyRDZASUMTest(),\
                        checks.VerifyRDZNRM2Test(),\
                        checks.VerifyRSSUMTest(),\
                        checks.VerifyRSASUMTest(),\
                        checks.VerifyRSNRM2Test(),\
                        checks.VerifyRCSUMTest(),\
                        checks.VerifyRSCASUMTest(),\
                        checks.VerifyRSCNRM2Test()],\
                       ["N", "incX", "f"],\
                       [[4095], [1, 4],\
                        ["rand",\
                         "2*rand-1",\
                         "rand+(rand-1)",\
                         "normal",\
                         "sine",\
                         "small+grow*big",\
                         "small+rand*big",\
                         "rand_cond"]])

check_suite.add_checks([checks.VerifyRDDOTTest(),\
                        checks.VerifyRZDOTUTest(),\
                        checks.VerifyRZDOTCTest(),\
                        checks.VerifyRSDOTTest(),\
                        checks.VerifyRCDOTUTest(),\
                        checks.VerifyRCDOTCTest()],\
                       ["N", "incX", "incY", "f"],\
                       [[4095], [1, 4], [1, 4],\
                        ["rand",\
                         "2*rand-1",\
                         "rand+(rand-1)",\
                         "normal",\
                         "sine",\
                         "small+grow*big",\
                         "small+rand*big",\
                         "rand_cond"]])

check_suite.add_checks([checks.ValidateExternalRDSUMTest(),\
                        checks.ValidateExternalRDASUMTest(),\
                        checks.ValidateExternalRDNRM2Test(),\
                        checks.ValidateExternalRDDOTTest()],\
                       ["i"],\
                       [[os.path.join(check_dir, "data/d_2_rand-1_N4095.dat"),\
                         os.path.join(check_dir, "data/d_normal_N4095.dat"),\
                         os.path.join(check_dir, "data/d_rand+(rand-1)_N4095.dat"),\
                         os.path.join(check_dir, "data/d_rand_N4095.dat"),\
                         os.path.join(check_dir, "data/d_rand_cond_N4095.dat"),\
                         os.path.join(check_dir, "data/d_sine_N4095.dat"),\
                         os.path.join(check_dir, "data/d_small+grow_big_N4095.dat"),\
                         os.path.join(check_dir, "data/d_small+rand_big_N4095.dat")]])

check_suite.add_checks([checks.ValidateExternalRZSUMTest(),\
                        checks.ValidateExternalRDZASUMTest(),\
                        checks.ValidateExternalRDZNRM2Test(),\
                        checks.ValidateExternalRZDOTUTest(),\
                        checks.ValidateExternalRZDOTCTest()],\
                       ["i"],\
                       [[os.path.join(check_dir, "data/z_2_rand-1_N4095.dat"),\
                         os.path.join(check_dir, "data/z_normal_N4095.dat"),\
                         os.path.join(check_dir, "data/z_rand+(rand-1)_N4095.dat"),\
                         os.path.join(check_dir, "data/z_rand_N4095.dat"),\
                         os.path.join(check_dir, "data/z_rand_cond_N4095.dat"),\
                         os.path.join(check_dir, "data/z_sine_N4095.dat"),\
                         os.path.join(check_dir, "data/z_small+grow_big_N4095.dat"),\
                         os.path.join(check_dir, "data/z_small+rand_big_N4095.dat")]])

check_suite.add_checks([checks.ValidateExternalRSSUMTest(),\
                        checks.ValidateExternalRSASUMTest(),\
                        checks.ValidateExternalRSNRM2Test(),\
                        checks.ValidateExternalRSDOTTest()],\
                       ["i"],\
                       [[os.path.join(check_dir, "data/s_2_rand-1_N4095.dat"),\
                         os.path.join(check_dir, "data/s_normal_N4095.dat"),\
                         os.path.join(check_dir, "data/s_rand+(rand-1)_N4095.dat"),\
                         os.path.join(check_dir, "data/s_rand_N4095.dat"),\
                         os.path.join(check_dir, "data/s_rand_cond_N4095.dat"),\
                         os.path.join(check_dir, "data/s_sine_N4095.dat"),\
                         os.path.join(check_dir, "data/s_small+grow_big_N4095.dat"),\
                         os.path.join(check_dir, "data/s_small+rand_big_N4095.dat")]])

check_suite.add_checks([checks.ValidateExternalRCSUMTest(),\
                        checks.ValidateExternalRSCASUMTest(),\
                        checks.ValidateExternalRSCNRM2Test(),\
                        checks.ValidateExternalRCDOTUTest(),\
                        checks.ValidateExternalRCDOTCTest()],\
                       ["i"],\
                       [[os.path.join(check_dir, "data/c_2_rand-1_N4095.dat"),\
                         os.path.join(check_dir, "data/c_normal_N4095.dat"),\
                         os.path.join(check_dir, "data/c_rand+(rand-1)_N4095.dat"),\
                         os.path.join(check_dir, "data/c_rand_N4095.dat"),\
                         os.path.join(check_dir, "data/c_rand_cond_N4095.dat"),\
                         os.path.join(check_dir, "data/c_sine_N4095.dat"),\
                         os.path.join(check_dir, "data/c_small+grow_big_N4095.dat"),\
                         os.path.join(check_dir, "data/c_small+rand_big_N4095.dat")]])

check_harness = harness.Harness("check")
check_harness.add_suite(check_suite)
check_harness.run()
