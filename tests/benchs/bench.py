import tests.benchs.benchs as benchs
import tests.harness.harness as harness

bench_harness = harness.Harness("bench")
#attribute = "%peak"
attribute = "time"
#attribute = "perf"
#attribute = "freq"

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDSUMTest(), benchs.BenchDSUMTest()], ["N", "fold"], [[2**i for i in range(3,13)], [3]], attribute, silent_flags="--FillA rand --ScaleX 0"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDDOTTest(), benchs.BenchDDOTTest()], ["N", "fold"], [[2**i for i in range(3,13)], [3]], attribute, silent_flags="--FillA rand --ScaleX 0"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDGEMVTest(), benchs.BenchDGEMVTest()], [("N", "M"), "fold", "Order"], [[(2048, 2048)], [3], ["ColMajor", "RowMajor"]], attribute, silent_flags="--FillA rand --FillX rand --ScaleX 0"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDGEMMTest(), benchs.BenchDGEMMTest()], [("N", "M", "K"), "fold", "TransA", "TransB"], [[(512, 512, 4096)], [3], ["Trans", "NoTrans"], ["Trans", "NoTrans"]], attribute, silent_flags="--FillA rand --FillB rand --Order ColMajor --ScaleA 0"))

bench_harness.run()
