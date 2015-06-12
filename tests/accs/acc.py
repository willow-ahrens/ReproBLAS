import tests.accs.accs as accs
import tests.harness.harness as harness

acc_harness = harness.Harness("acc")
params = ["N", "a", "u", "fold", "f"]
ranges = [[4096], [100], [0], [3], ["+-big", "normal", "rand"]]
attribute = "ratio"
#acc_harness.add_suite(accs.AccSuite([accs.AccRDSUMTest(), accs.AccRDASUMTest(), accs.AccRDNRM2Test(), accs.AccRDDOTTest()], params, ranges, attribute))
#acc_harness.add_suite(accs.AccSuite([accs.AccRSSUMTest(), accs.AccRSASUMTest(), accs.AccRSNRM2Test(), accs.AccRSDOTTest()], params, ranges, attribute))
#acc_harness.add_suite(accs.AccSuite([accs.AccRZSUMTest(), accs.AccRDZASUMTest(), accs.AccRDZNRM2Test(), accs.AccRZDOTUTest(), accs.AccRZDOTCTest()], params, ranges, attribute))
#acc_harness.add_suite(accs.AccSuite([accs.AccRCSUMTest(), accs.AccRSCASUMTest(), accs.AccRSCNRM2Test(), accs.AccRCDOTUTest(), accs.AccRCDOTCTest()], params, ranges, attribute))
acc_harness.add_suite(accs.AccSuite([accs.AccRDSUMTest(), accs.AccRZSUMTest(), accs.AccRSSUMTest(), accs.AccRCSUMTest()], params, ranges, attribute))
acc_harness.run()
