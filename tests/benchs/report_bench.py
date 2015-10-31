import tests.benchs.benchs as benchs
import tests.accs.accs as accs
import tests.harness.harness as harness

bench_harness = harness.Harness("bench")

fold = 3

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDSUMTest(), benchs.BenchRSSUMTest()], ["N", "fold", "FillX", "a"], [[16 * 1024], [fold], ["rand", "full_range"], [1000]], "time"))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDSUMTest(), benchs.BenchRDSUMTest()], ["N", "fold", "FillX"], [[2**i for i in range(6, 13)], [fold], ["normal"]], "peak_time"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDSUMTest(), benchs.BenchRDSUMTest()], ["N", "fold", "FillX"], [[2**i for i in range(6, 13)], [fold], ["normal"]], "time"))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDDOTTest(), benchs.BenchRDDOTTest()], ["N", "fold", "FillX"], [[2**i for i in range(6, 13)], [fold], ["normal"]], "peak_time"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDDOTTest(), benchs.BenchRDDOTTest()], ["N", "fold", "FillX"], [[2**i for i in range(6, 13)], [fold], ["normal"]], "time"))

trials = 10

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMVTest(), benchs.BenchRDGEMVTest()], [("N", "M"), "TransA", "fold", "FillA", 'a'], [[(2**i, 2**i) for i in range(6, 13)], ["NoTrans"], [3], ["normal"], [0]], "peak_time", silent_flags="--Order ColMajor"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMVTest(), benchs.BenchRDGEMVTest()], [("N", "M"), "TransA", "fold", "FillA", 'a'], [[(2**i, 2**i) for i in range(6, 13)], ["NoTrans"], [3], ["normal"], [trials]], "time", silent_flags="--Order ColMajor"))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMVTest(), benchs.BenchRDGEMVTest()], [("N", "M"), "TransA", "fold", "FillA", 'a'], [[(2**i, 2**i) for i in range(6, 13)], ["Trans"], [3], ["normal"], [0]], "peak_time", silent_flags="--Order ColMajor"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMVTest(), benchs.BenchRDGEMVTest()], [("N", "M"), "TransA", "fold", "FillA", 'a'], [[(2**i, 2**i) for i in range(6, 13)], ["Trans"], [3], ["normal"], [trials]], "time", silent_flags="--Order ColMajor"))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMMTest(), benchs.BenchRDGEMMTest()], [("N", "M", "K"), ("TransA", "TransB"), "fold", "FillA", 'a'], [[(2**i, 2**i, 2**i) for i in range(6, 13)], [("NoTrans", "NoTrans")], [3], ["normal"], [0]], "peak_time", silent_flags="--Order ColMajor"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchDGEMMTest(), benchs.BenchRDGEMMTest()], [("N", "M", "K"), ("TransA", "TransB"), "fold", "FillA", 'a'], [[(2**i, 2**i, 2**i) for i in range(6, 13)], [("NoTrans", "NoTrans")], [3], ["normal"], [trials]], "time", silent_flags="--Order ColMajor"))

bench_harness.run()
