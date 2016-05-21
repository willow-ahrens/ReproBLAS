import tests.accs.accs as accs
import tests.harness.harness as harness

acc_harness = harness.Harness("acc")
params = ["N", "a", "u", "fold", "f", "RealScaleX"]
ranges = [[32768], [100], [0], [4], ("rand_cond", "N_cond+rand"), [2.0**6, 2.0**12, 2.0**24, 2.0**48]]

attribute = "max_ratio(e)"
#attribute = "max_ratio"
#attribute = "max_ratio(e)"
#attribute = "med_ratio"
#attribute = "med_ratio(e)"
#attribute = "min_ratio"
#attribute = "min_ratio(e)"

print("Measuring: {}".format(attribute))

acc_harness.add_suite(accs.AccSuite([accs.AccRDSUMTest(), accs.AccRDASUMTest(), accs.AccRDNRM2Test(), accs.AccRDDOTTest()], params, ranges, attribute))
acc_harness.add_suite(accs.AccSuite([accs.AccRSSUMTest(), accs.AccRSASUMTest(), accs.AccRSNRM2Test(), accs.AccRSDOTTest()], params, ranges, attribute))
acc_harness.add_suite(accs.AccSuite([accs.AccRZSUMTest(), accs.AccRDZASUMTest(), accs.AccRDZNRM2Test(), accs.AccRZDOTUTest(), accs.AccRZDOTCTest()], params, ranges, attribute))
acc_harness.add_suite(accs.AccSuite([accs.AccRCSUMTest(), accs.AccRSCASUMTest(), accs.AccRSCNRM2Test(), accs.AccRCDOTUTest(), accs.AccRCDOTCTest()], params, ranges, attribute))
acc_harness.run()
