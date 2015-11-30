import tests.benchs.benchs as benchs
import tests.accs.accs as accs
import tests.harness.harness as harness

bench_harness = harness.Harness("bench")

fold = 3

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDSUMTest(), benchs.BenchRSSUMTest()], ["N", "fold", "FillX", "a"], [[8 * 1024], [fold], ["rand", "half_range"], [1000]], "time"))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDSUMTest(), benchs.BenchRSSUMTest()], ["N", "fold", "FillX", "a"], [[8 * 1024], [fold], ["rand", "half_range"], [1000]], "peak"))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDSUMTest(), benchs.BenchRDSUMTest()], ["N", "fold", "FillX"], [[2**i for i in range(6, 13)], [fold], ["normal"]], "peak_time"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDSUMTest(), benchs.BenchRDSUMTest()], ["N", "fold", "FillX"], [[2**i for i in range(6, 13)], [fold], ["normal"]], "time"))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDDOTTest(), benchs.BenchRDDOTTest()], ["N", "fold", "FillX", "FillY"], [[2**i for i in range(6, 13)], [fold], ["normal"], ["normal"]], "peak_time"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDDOTTest(), benchs.BenchRDDOTTest()], ["N", "fold", "FillX", "FillY"], [[2**i for i in range(6, 13)], [fold], ["normal"], ["normal"]], "time"))

trials = 100

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMVTest(), benchs.BenchRDGEMVTest()], [("N", "M"), "TransA", "fold", "FillA", "FillX", "FillY", 'a'], [[(2**i, 2**i) for i in range(6, 13)], ["NoTrans"], [3], ["normal"], ["constant"], ["constant"], [0]], "peak_time", silent_flags="--Order ColMajor"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMVTest(), benchs.BenchRDGEMVTest()], [("N", "M"), "TransA", "fold", "FillA", "FillX", "FillY", 'a'], [[(2**i, 2**i) for i in range(6, 13)], ["NoTrans"], [3], ["normal"], ["constant"], ["constant"], [trials]], "time", silent_flags="--Order ColMajor"))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMVTest(), benchs.BenchRDGEMVTest()], [("N", "M"), "TransA", "fold", "FillA", "FillB", "FillC", 'a'], [[(2**i, 2**i) for i in range(6, 13)], ["Trans"], [3], ["normal"], ["constant"], ["constant"], [0]], "peak_time", silent_flags="--Order ColMajor"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMVTest(), benchs.BenchRDGEMVTest()], [("N", "M"), "TransA", "fold", "FillA", "FillB", "FillC", 'a'], [[(2**i, 2**i) for i in range(6, 13)], ["Trans"], [3], ["normal"], ["constant"], ["constant"], [trials]], "time", silent_flags="--Order ColMajor"))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMMTest(), benchs.BenchRDGEMMTest()], [("N", "M", "K"), ("TransA", "TransB"), "fold", "FillA", 'a'], [[(2**i, 2**i, 2**i) for i in range(6, 13)], [("NoTrans", "NoTrans")], [3], ["normal"], [0]], "peak_time", silent_flags="--Order ColMajor"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMMTest(), benchs.BenchRDGEMMTest()], [("N", "M", "K"), ("TransA", "TransB"), "fold", "FillA", 'a'], [[(2**i, 2**i, 2**i) for i in range(6, 13)], [("NoTrans", "NoTrans")], [3], ["normal"], [trials]], "time", silent_flags="--Order ColMajor"))

bench_harness.run()
