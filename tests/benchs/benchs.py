import copy
import itertools

import scripts.terminal as terminal
import tests.harness.harness as harness

class BenchSuite(harness.MetricSuite):
  pass

class BenchTest(harness.MetricTest):
  def parse_output(self):
    if self.attribute == "freq":
      self.result = (self.output["trials"] * self.output["input"]) / self.output["time"]
    elif self.attribute == "peak":
      self.result = (terminal.get_peak_time(self.output) * self.output["trials"])/self.output["time"]
    elif self.attribute == "%peak":
      self.result = (100.0 * terminal.get_peak_time(self.output) * self.output["trials"])/self.output["time"]
    else: #self.attribute == "time":
      self.result = self.output["time"] / self.output["trials"]

class BenchCAMAXTest(BenchTest):
  name = "CAMAX"
  executable = "tests/benchs/bench_camax"

class BenchCAMAXMTest(BenchTest):
  name = "CAMAXM"
  executable = "tests/benchs/bench_camaxm"

class BenchCDOTCTest(BenchTest):
  name = "CDOTC"
  executable = "tests/benchs/bench_cdotc"

class BenchCDOTUTest(BenchTest):
  name = "CDOTU"
  executable = "tests/benchs/bench_cdotu"

class BenchDAMAXTest(BenchTest):
  name = "DAMAX"
  executable = "tests/benchs/bench_damax"

class BenchDAMAXMTest(BenchTest):
  name = "DAMAXM"
  executable = "tests/benchs/bench_damaxm"

class BenchDASUMTest(BenchTest):
  name = "DASUM"
  executable = "tests/benchs/bench_dasum"

class BenchDDOTTest(BenchTest):
  name = "DDOT"
  executable = "tests/benchs/bench_ddot"

class BenchDNRM2Test(BenchTest):
  name = "DNRM2"
  executable = "tests/benchs/bench_dnrm2"

class BenchDZASUMTest(BenchTest):
  name = "DZASUM"
  executable = "tests/benchs/bench_dzasum"

class BenchDZNRM2Test(BenchTest):
  name = "DZNRM2"
  executable = "tests/benchs/bench_dznrm2"

class BenchICAMAXTest(BenchTest):
  name = "ICAMAX"
  executable = "tests/benchs/bench_icamax"

class BenchIDAMAXTest(BenchTest):
  name = "IDAMAX"
  executable = "tests/benchs/bench_idamax"

class BenchISAMAXTest(BenchTest):
  name = "ISAMAX"
  executable = "tests/benchs/bench_isamax"

class BenchIZAMAXTest(BenchTest):
  name = "IZAMAX"
  executable = "tests/benchs/bench_izamax"

class BenchRCDOTCTest(BenchTest):
  name = "RCDOTC"
  executable = "tests/benchs/bench_rcdotc"

class BenchRCDOTUTest(BenchTest):
  name = "RCDOTU"
  executable = "tests/benchs/bench_rcdotu"

class BenchRCSUMTest(BenchTest):
  name = "RCSUM"
  executable = "tests/benchs/bench_rcsum"

class BenchRDASUMTest(BenchTest):
  name = "RDASUM"
  executable = "tests/benchs/bench_rdasum"

class BenchRDDOTTest(BenchTest):
  name = "RDDOT"
  executable = "tests/benchs/bench_rddot"

class BenchRDNRM2Test(BenchTest):
  name = "RDNRM2"
  executable = "tests/benchs/bench_rdnrm2"

class BenchRDSUMTest(BenchTest):
  name = "RDSUM"
  executable = "tests/benchs/bench_rdsum"

class BenchRDZASUMTest(BenchTest):
  name = "RDZASUM"
  executable = "tests/benchs/bench_rdzasum"

class BenchRDZNRM2Test(BenchTest):
  name = "RDZNRM2"
  executable = "tests/benchs/bench_rdznrm2"

class BenchRSASUMTest(BenchTest):
  name = "RSASUM"
  executable = "tests/benchs/bench_rsasum"

class BenchRSCASUMTest(BenchTest):
  name = "RSCASUM"
  executable = "tests/benchs/bench_rscasum"

class BenchRSCNRM2Test(BenchTest):
  name = "RSCNRM2"
  executable = "tests/benchs/bench_rscnrm2"

class BenchRSDOTTest(BenchTest):
  name = "RSDOT"
  executable = "tests/benchs/bench_rsdot"

class BenchRSNRM2Test(BenchTest):
  name = "RSNRM2"
  executable = "tests/benchs/bench_rsnrm2"

class BenchRSSUMTest(BenchTest):
  name = "RSSUM"
  executable = "tests/benchs/bench_rssum"

class BenchRZDOTCTest(BenchTest):
  name = "RZDOTC"
  executable = "tests/benchs/bench_rzdotc"

class BenchRZDOTUTest(BenchTest):
  name = "RZDOTU"
  executable = "tests/benchs/bench_rzdotu"

class BenchRZSUMTest(BenchTest):
  name = "RZSUM"
  executable = "tests/benchs/bench_rzsum"

class BenchSAMAXTest(BenchTest):
  name = "SAMAX"
  executable = "tests/benchs/bench_samax"

class BenchSAMAXMTest(BenchTest):
  name = "SAMAXM"
  executable = "tests/benchs/bench_samaxm"

class BenchSASUMTest(BenchTest):
  name = "SASUM"
  executable = "tests/benchs/bench_sasum"

class BenchSCASUMTest(BenchTest):
  name = "SCASUM"
  executable = "tests/benchs/bench_scasum"

class BenchSCNRM2Test(BenchTest):
  name = "SCNRM2"
  executable = "tests/benchs/bench_scnrm2"

class BenchSDOTTest(BenchTest):
  name = "SDOT"
  executable = "tests/benchs/bench_sdot"

class BenchSNRM2Test(BenchTest):
  name = "SNRM2"
  executable = "tests/benchs/bench_snrm2"

class BenchZAMAXTest(BenchTest):
  name = "ZAMAX"
  executable = "tests/benchs/bench_zamax"

class BenchZAMAXMTest(BenchTest):
  name = "ZAMAXM"
  executable = "tests/benchs/bench_zamaxm"

class BenchZDOTCTest(BenchTest):
  name = "ZDOTC"
  executable = "tests/benchs/bench_zdotc"

class BenchZDOTUTest(BenchTest):
  name = "ZDOTU"
  executable = "tests/benchs/bench_zdotu"

class BenchZGEMVTest(BenchTest):
  name = "ZGEMV"
  executable = "tests/benchs/bench_zgemv"

class BenchDDICONVTest(BenchTest):
  name = "DDICONV"
  executable = "tests/benchs/bench_ddiconv"

class BenchZZICONVTest(BenchTest):
  name = "ZZICONV"
  executable = "tests/benchs/bench_zziconv"

class BenchSSICONVTest(BenchTest):
  name = "SSICONV"
  executable = "tests/benchs/bench_ssiconv"

class BenchCCICONVTest(BenchTest):
  name = "CCICONV"
  executable = "tests/benchs/bench_cciconv"

class BenchDIDIADDTest(BenchTest):
  name = "DIDIADD"
  executable = "tests/benchs/bench_didiadd"

class BenchZIZIADDTest(BenchTest):
  name = "ZIZIADD"
  executable = "tests/benchs/bench_ziziadd"

class BenchSISIADDTest(BenchTest):
  name = "SISIADD"
  executable = "tests/benchs/bench_sisiadd"

class BenchCICIADDTest(BenchTest):
  name = "CICIADD"
  executable = "tests/benchs/bench_ciciadd"

class BenchRDGEMVTest(BenchTest):
  name = "RDGEMV"
  executable = "tests/benchs/bench_rdgemv"

class BenchRDGEMMTest(BenchTest):
  name = "RDGEMM"
  executable = "tests/benchs/bench_rdgemm"

class BenchDGEMMTest(BenchTest):
  name = "DGEMM"
  executable = "tests/benchs/bench_dgemm"

class BenchDGEMVTest(BenchTest):
  name = "DGEMV"
  executable = "tests/benchs/bench_dgemv"

class BenchRZGEMVTest(BenchTest):
  name = "RZGEMV"
  executable = "tests/benchs/bench_rzgemv"

all_benchs = {"bench_camax": (BenchCAMAXTest, ""),\
              "bench_camaxm": (BenchCAMAXMTest, ""),\
              "bench_cciconv": (BenchCCICONVTest, ""),\
              "bench_cdotc": (BenchCDOTCTest, ""),\
              "bench_cdotu": (BenchCDOTUTest, ""),\
              "bench_ciciadd": (BenchCICIADDTest, ""),\
              "bench_damax": (BenchDAMAXTest, ""),\
              "bench_damaxm": (BenchDAMAXMTest, ""),\
              "bench_dasum": (BenchDASUMTest, ""),\
              "bench_ddiconv": (BenchDDICONVTest, ""),\
              "bench_ddot": (BenchDDOTTest, ""),\
              "bench_didiadd": (BenchDIDIADDTest, ""),\
              "bench_dnrm2": (BenchDNRM2Test, ""),\
              "bench_dzasum": (BenchDZASUMTest, ""),\
              "bench_dznrm2": (BenchDZNRM2Test, ""),\
              "bench_icamax": (BenchICAMAXTest, ""),\
              "bench_idamax": (BenchIDAMAXTest, ""),\
              "bench_isamax": (BenchISAMAXTest, ""),\
              "bench_izamax": (BenchIZAMAXTest, ""),\
              "bench_samax": (BenchSAMAXTest, ""),\
              "bench_samaxm": (BenchSAMAXMTest, ""),\
              "bench_sasum": (BenchSASUMTest, ""),\
              "bench_scasum": (BenchSCASUMTest, ""),\
              "bench_scnrm2": (BenchSCNRM2Test, ""),\
              "bench_sdot": (BenchSDOTTest, ""),\
              "bench_sisiadd": (BenchSISIADDTest, ""),\
              "bench_snrm2": (BenchSNRM2Test, ""),\
              "bench_ssiconv": (BenchSSICONVTest, ""),\
              "bench_zamax": (BenchZAMAXTest, ""),\
              "bench_zamaxm": (BenchZAMAXMTest, ""),\
              "bench_zdotc": (BenchZDOTCTest, ""),\
              "bench_zdotu": (BenchZDOTUTest, ""),\
              "bench_ziziadd": (BenchZIZIADDTest, ""),\
              "bench_zziconv": (BenchZZICONVTest, ""),\
              "bench_dgemv": (BenchDGEMVTest, ""),\
              "bench_zgemv": (BenchZGEMVTest, ""),\
             }

for i in range(terminal.get_simaxfold() + 1):
  all_benchs.update({"bench_rcdotc_fold_{}".format(i): (BenchRCDOTCTest, "--fold {}".format(i)),\
                     "bench_rcdotu_fold_{}".format(i): (BenchRCDOTUTest, "--fold {}".format(i)),\
                     "bench_rcsum_fold_{}".format(i): (BenchRCSUMTest, "--fold {}".format(i)),\
                     "bench_rsasum_fold_{}".format(i): (BenchRSASUMTest, "--fold {}".format(i)),\
                     "bench_rscasum_fold_{}".format(i): (BenchRSCASUMTest, "--fold {}".format(i)),\
                     "bench_rscnrm2_fold_{}".format(i): (BenchRSCNRM2Test, "--fold {}".format(i)),\
                     "bench_rsdot_fold_{}".format(i): (BenchRSDOTTest, "--fold {}".format(i)),\
                     "bench_rsnrm2_fold_{}".format(i): (BenchRSNRM2Test, "--fold {}".format(i)),\
                     "bench_rssum_fold_{}".format(i): (BenchRSSUMTest, "--fold {}".format(i)),\
                    })

for i in range(terminal.get_dimaxfold() + 1):
  all_benchs.update({"bench_rdasum_fold_{}".format(i): (BenchRDASUMTest, "--fold {}".format(i)),\
                     "bench_rddot_fold_{}".format(i): (BenchRDDOTTest, "--fold {}".format(i)),\
                     "bench_rdnrm2_fold_{}".format(i): (BenchRDNRM2Test, "--fold {}".format(i)),\
                     "bench_rdsum_fold_{}".format(i): (BenchRDSUMTest, "--fold {}".format(i)),\
                     "bench_rdgemv_fold_{}".format(i): (BenchRDGEMVTest, "--fold {}".format(i)),\
                     "bench_rdzasum_fold_{}".format(i): (BenchRDZASUMTest, "--fold {}".format(i)),\
                     "bench_rdznrm2_fold_{}".format(i): (BenchRDZNRM2Test, "--fold {}".format(i)),\
                     "bench_rzdotc_fold_{}".format(i): (BenchRZDOTCTest, "--fold {}".format(i)),\
                     "bench_rzdotu_fold_{}".format(i): (BenchRZDOTUTest, "--fold {}".format(i)),\
                     "bench_rzsum_fold_{}".format(i): (BenchRZSUMTest, "--fold {}".format(i)),\
                     "bench_rzgemv_fold_{}".format(i): (BenchRZGEMVTest, "--fold {}".format(i)),\
                    })
