import tests.benchs.benchs as benchs
import tests.harness.harness as harness

bench_harness = harness.Harness("bench")
params = [("N", "M", "K"), "TransA", "TransB", "a", "FillA", "FillB", "FillC","fold", "preN"]
ranges = [[(2048, 2048, 2048), (256, 256, 256)], ["Trans", "NoTrans"], ["Trans", "NoTrans"], [1], ["rand"], ["rand"], ["rand"], [3], [1024]]
attribute = "%peak"
#attribute = "time"
#attribute = "freq"
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDSUMTest(), benchs.BenchRDASUMTest(), benchs.BenchRDNRM2Test(), benchs.BenchRDDOTTest()], params, ranges, attribute))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRSSUMTest(), benchs.BenchRSASUMTest(), benchs.BenchRSNRM2Test(), benchs.BenchRSDOTTest()], params, ranges, attribute))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRZSUMTest(), benchs.BenchRDZASUMTest(), benchs.BenchRDZNRM2Test(), benchs.BenchRZDOTUTest(), benchs.BenchRZDOTCTest()], params, ranges, attribute))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRCSUMTest(), benchs.BenchRSCASUMTest(), benchs.BenchRSCNRM2Test(), benchs.BenchRCDOTUTest(), benchs.BenchRCDOTCTest()], params, ranges, attribute))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDDICONVTest(), benchs.BenchZZICONVTest(), benchs.BenchSSICONVTest(), benchs.BenchCCICONVTest()], params, ranges, attribute))
#bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDGEMVTest(), benchs.BenchDGEMVTest()], params + ["Order"], ranges + [["ColMajor", "RowMajor"]], attribute))
#bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDGEMMTest(), benchs.BenchDGEMMTest()], params + ["Order"], ranges + [["RowMajor"]], attribute))
#bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDGEMMTest(), benchs.BenchDGEMMTest()], params + ["Order"], ranges + [["ColMajor"]], attribute))
#bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRZGEMVTest(), benchs.BenchZGEMVTest()], params + ["Order"], ranges + [["RowMajor"]], attribute))
#bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDIDIADDTest(), benchs.BenchRDSUMTest()], params, ranges, attribute))
#bench_harness.add_suite(benchs.BenchSuite([benchs.BenchZIZIADDTest(), benchs.BenchRZSUMTest()], params, ranges, attribute))
#bench_harness.add_suite(benchs.BenchSuite([benchs.BenchSISIADDTest(), benchs.BenchRSSUMTest()], params, ranges, attribute))
#bench_harness.add_suite(benchs.BenchSuite([benchs.BenchCICIADDTest(), benchs.BenchRCSUMTest()], params, ranges, attribute))
bench_harness.run()
