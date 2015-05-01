import tests.benchs.benchs as benchs
import tests.harness.harness as harness

bench_harness = harness.Harness("bench")
params = ["N", "a", "u"]
ranges = [[4], [100000], [0]]
attribute = "peak"
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDGEMVTest(), benchs.BenchDGEMVTest()], ["N", "M", "O", "T"], [[512], [512], ["RowMajor"], ["NoTrans"]], attribute))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDSUMTest(), benchs.BenchRDASUMTest(), benchs.BenchRDNRM2Test(), benchs.BenchRDDOTTest()], params, ranges, attribute))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRSSUMTest(), benchs.BenchRSASUMTest(), benchs.BenchRSNRM2Test(), benchs.BenchRSDOTTest()], params, ranges, attribute))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRZSUMTest(), benchs.BenchRDZASUMTest(), benchs.BenchRDZNRM2Test(), benchs.BenchRZDOTUTest(), benchs.BenchRZDOTCTest()], params, ranges, attribute))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRCSUMTest(), benchs.BenchRSCASUMTest(), benchs.BenchRSCNRM2Test(), benchs.BenchRCDOTUTest(), benchs.BenchRCDOTCTest()], params, ranges, attribute))

bench_harness.run()
