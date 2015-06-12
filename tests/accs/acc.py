import tests.accs.accs as accs
import tests.harness.harness as harness

acc_harness = harness.Harness("acc")
params = ["N", "a", "u", "fold", "f", "RealScaleX"]
ranges = [[4096], [1000], [0], [4], ("normal", "rand_cond"), [2.0**5, 2.0**10, 2.0**20, 2.0**40, 2.0**80]]
attribute = "med_ratio(e)"
#acc_harness.add_suite(accs.AccSuite([accs.AccRDSUMTest(), accs.AccRDASUMTest(), accs.AccRDNRM2Test(), accs.AccRDDOTTest()], params, ranges, attribute))
#acc_harness.add_suite(accs.AccSuite([accs.AccRSSUMTest(), accs.AccRSASUMTest(), accs.AccRSNRM2Test(), accs.AccRSDOTTest()], params, ranges, attribute))
#acc_harness.add_suite(accs.AccSuite([accs.AccRZSUMTest(), accs.AccRDZASUMTest(), accs.AccRDZNRM2Test(), accs.AccRZDOTUTest(), accs.AccRZDOTCTest()], params, ranges, attribute))
#acc_harness.add_suite(accs.AccSuite([accs.AccRCSUMTest(), accs.AccRSCASUMTest(), accs.AccRSCNRM2Test(), accs.AccRCDOTUTest(), accs.AccRCDOTCTest()], params, ranges, attribute))
acc_harness.add_suite(accs.AccSuite([accs.AccRDSUMTest(), accs.AccRZSUMTest()], params, ranges, attribute))
acc_harness.run()
