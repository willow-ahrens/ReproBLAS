import tests.benchs.benchs as benchs
import tests.harness.harness as harness

bench_harness = harness.Harness("bench")

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDSUMTest(), benchs.BenchRSSUMTest()], ["N", "fold", "FillX", "a"], [[1024 * 1024], [3], ["rand", "full_range"], [10]], "freq"))

bench_harness.run()
