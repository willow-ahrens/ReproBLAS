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

class BenchRDGEMVTest(BenchTest):
  name = "RDGEMV"
  executable = "tests/benchs/bench_rdgemv"

class BenchDGEMVTest(BenchTest):
  name = "DGEMV"
  executable = "tests/benchs/bench_dgemv"

class BenchPRDGEMVTest(BenchTest):
  name = "PRDGEMV"
  executable = "tests/benchs/bench_prdgemv"

class BenchPRBDGEMVTest(BenchTest):
  name = "PRBDGEMV"
  executable = "tests/benchs/bench_prbdgemv"

class BenchPDGEMVTest(BenchTest):
  name = "PDGEMV"
  executable = "tests/benchs/bench_pdgemv"

all_benchs = {"bench_camax": BenchCAMAXTest,\
              "bench_camaxm": BenchCAMAXMTest,\
              "bench_cdotc": BenchCDOTCTest,\
              "bench_cdotu": BenchCDOTUTest,\
              "bench_damax": BenchDAMAXTest,\
              "bench_damaxm": BenchDAMAXMTest,\
              "bench_dasum": BenchDASUMTest,\
              "bench_ddot": BenchDDOTTest,\
              "bench_dnrm2": BenchDNRM2Test,\
              "bench_dzasum": BenchDZASUMTest,\
              "bench_dznrm2": BenchDZNRM2Test,\
              "bench_icamax": BenchICAMAXTest,\
              "bench_idamax": BenchIDAMAXTest,\
              "bench_isamax": BenchISAMAXTest,\
              "bench_izamax": BenchIZAMAXTest,\
              "bench_rcdotc": BenchRCDOTCTest,\
              "bench_rcdotu": BenchRCDOTUTest,\
              "bench_rcsum": BenchRCSUMTest,\
              "bench_rdasum": BenchRDASUMTest,\
              "bench_rddot": BenchRDDOTTest,\
              "bench_rdnrm2": BenchRDNRM2Test,\
              "bench_rdsum": BenchRDSUMTest,\
              "bench_rdzasum": BenchRDZASUMTest,\
              "bench_rdznrm2": BenchRDZNRM2Test,\
              "bench_rsasum": BenchRSASUMTest,\
              "bench_rscasum": BenchRSCASUMTest,\
              "bench_rscnrm2": BenchRSCNRM2Test,\
              "bench_rsdot": BenchRSDOTTest,\
              "bench_rsnrm2": BenchRSNRM2Test,\
              "bench_rssum": BenchRSSUMTest,\
              "bench_rzdotc": BenchRZDOTCTest,\
              "bench_rzdotu": BenchRZDOTUTest,\
              "bench_rzsum": BenchRZSUMTest,\
              "bench_samax": BenchSAMAXTest,\
              "bench_samaxm": BenchSAMAXMTest,\
              "bench_sasum": BenchSASUMTest,\
              "bench_scasum": BenchSCASUMTest,\
              "bench_scnrm2": BenchSCNRM2Test,\
              "bench_sdot": BenchSDOTTest,\
              "bench_snrm2": BenchSNRM2Test,\
              "bench_zamax": BenchZAMAXTest,\
              "bench_zamaxm": BenchZAMAXMTest,\
              "bench_zdotc": BenchZDOTCTest,\
              "bench_zdotu": BenchZDOTUTest,\
              "bench_dgemv": BenchDGEMVTest,\
              "bench_rdgemv": BenchRDGEMVTest}



#TODO THIS IS A HACK!!! we need to support many-fold testing.
new_benchs = {}

maxs = set(["bench_damax", "bench_damaxm", "bench_zamax", "bench_zamaxm", "bench_samax", "bench_samaxm", "bench_camax", "bench_camaxm"])
for key, val in all_benchs.items():
  if key in maxs:
    new_benchs[key] = val
  else:
    new_benchs[key + "_fold_3"] = val

all_benchs = new_benchs
