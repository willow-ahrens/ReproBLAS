import tests.benchs.benchs as benchs
import tests.harness.harness as harness

bench_harness = harness.Harness("bench")

attribute = "time(s)"
#attribute = "peak_time(s)"
#attribute = "perf(FLOP/s)"
#attribute = "peak_perf(FLOP/s)"
#attribute = "norm(Hz)"
#attribute = "peak_norm(Hz)"
#attribute = "freq(Hz)"
#attribute = "peak_freq(Hz)"
#attribute = "peak(%)"

print("Measuring: {}".format(attribute))

bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDSUMTest(), benchs.BenchRDASUMTest(), benchs.BenchRDNRM2Test(), benchs.BenchRDDOTTest()], ["N", "fold"], [[4096], [3]], attribute, silent_flags="--FillA rand"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRZSUMTest(), benchs.BenchRDZASUMTest(), benchs.BenchRDZNRM2Test(), benchs.BenchRZDOTUTest(), benchs.BenchRZDOTCTest()], ["N", "fold"], [[4096], [3]], attribute, silent_flags="--FillA rand"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRSSUMTest(), benchs.BenchRSASUMTest(), benchs.BenchRSNRM2Test(), benchs.BenchRSDOTTest()], ["N", "fold"], [[4096], [3]], attribute, silent_flags="--FillA rand"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRCSUMTest(), benchs.BenchRSCASUMTest(), benchs.BenchRSCNRM2Test(), benchs.BenchRCDOTUTest(), benchs.BenchRCDOTCTest()], ["N", "fold"], [[4096], [3]], attribute, silent_flags="--FillA rand"))
"""
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDGEMVTest(), benchs.BenchDGEMVTest()], [("N", "M"), "fold", "Order"], [[(2048, 2048)], [3], ["ColMajor", "RowMajor"]], attribute, silent_flags="--FillA rand --FillX rand"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRZGEMVTest(), benchs.BenchZGEMVTest()], [("N", "M"), "fold", "Order"], [[(2048, 2048)], [3], ["ColMajor", "RowMajor"]], attribute, silent_flags="--FillA rand --FillX rand"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRSGEMVTest(), benchs.BenchSGEMVTest()], [("N", "M"), "fold", "Order"], [[(2048, 2048)], [3], ["ColMajor", "RowMajor"]], attribute, silent_flags="--FillA rand --FillX rand"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRCGEMVTest(), benchs.BenchCGEMVTest()], [("N", "M"), "fold", "Order"], [[(2048, 2048)], [3], ["ColMajor", "RowMajor"]], attribute, silent_flags="--FillA rand --FillX rand"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRDGEMMTest(), benchs.BenchDGEMMTest()], [("N", "M", "K"), "fold", "TransA", "TransB"], [[(512, 512, 4096)], [3], ["Trans", "NoTrans"], ["Trans", "NoTrans"]], attribute, silent_flags="--FillA rand --FillB rand --Order ColMajor"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRZGEMMTest(), benchs.BenchZGEMMTest()], [("N", "M", "K"), "fold", "TransA", "TransB"], [[(512, 512, 4096)], [3], ["Trans", "NoTrans"], ["Trans", "NoTrans"]], attribute, silent_flags="--FillA rand --FillB rand --Order ColMajor"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRSGEMMTest(), benchs.BenchSGEMMTest()], [("N", "M", "K"), "fold", "TransA", "TransB"], [[(512, 512, 4096)], [3], ["Trans", "NoTrans"], ["Trans", "NoTrans"]], attribute, silent_flags="--FillA rand --FillB rand --Order ColMajor"))
bench_harness.add_suite(benchs.BenchSuite([benchs.BenchRCGEMMTest(), benchs.BenchCGEMMTest()], [("N", "M", "K"), "fold", "TransA", "TransB"], [[(512, 512, 4096)], [3], ["Trans", "NoTrans"], ["Trans", "NoTrans"]], attribute, silent_flags="--FillA rand --FillB rand --Order ColMajor"))
"""

bench_harness.run()
