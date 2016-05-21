import copy
import itertools
import sys

import scripts.terminal as terminal
import tests.harness.harness as harness

class BenchSuite(harness.MetricSuite):
  pass

class BenchTest(harness.MetricTest):
  def parse_output(self, output):
    if self.attribute == "time(s)":
      return output["time"] / output["trials"]
    elif self.attribute == "peak_time(s)":
      return terminal.get_peak_time(output)
    elif self.attribute == "perf(FLOP/s)":
      return (output["trials"] * terminal.get_flop_count(output)) / max(output["time"], sys.float_info.min)
    elif self.attribute == "peak_perf(FLOP/s)":
      return (output["trials"] * terminal.get_flop_count(output)) / max(terminal.get_peak_time(output), sys.float_info.min)
    elif self.attribute == "norm(Hz)":
      return (output["trials"] * output["normalizer"]) / max(output["time"], sys.float_info.min)
    elif self.attribute == "peak_norm(Hz)":
      return (output["trials"] * output["normalizer"]) / max(terminal.get_peak_time(output), sys.float_info.min)
    elif self.attribute == "freq(Hz)":
      return (output["trials"] * output["input"]) / max(output["time"], sys.float_info.min)
    elif self.attribute == "peak_freq(Hz)":
      return (output["trials"] * output["input"]) / max(terminal.get_peak_time(output), sys.float_info.min)
    elif self.attribute == "peak(%)":
      return (100.0 * terminal.get_peak_time(output) * output["trials"])/max(output["time"], sys.float_info.min)

class BenchCAMAXTest(BenchTest):
  name = "CAMAX"
  executable = "tests/benchs/bench_camax"

class BenchCAMAXMTest(BenchTest):
  name = "CAMAXM"
  executable = "tests/benchs/bench_camaxm"

class BenchCSUMTest(BenchTest):
  name = "CSUM"
  executable = "tests/benchs/bench_csum"

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

class BenchDSUMTest(BenchTest):
  name = "DSUM"
  executable = "tests/benchs/bench_dsum"

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

class BenchSSUMTest(BenchTest):
  name = "SSUM"
  executable = "tests/benchs/bench_ssum"

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

class BenchZSUMTest(BenchTest):
  name = "ZSUM"
  executable = "tests/benchs/bench_zsum"

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

class BenchRZGEMVTest(BenchTest):
  name = "RZGEMV"
  executable = "tests/benchs/bench_rzgemv"

class BenchRSGEMVTest(BenchTest):
  name = "RSGEMV"
  executable = "tests/benchs/bench_rsgemv"

class BenchRCGEMVTest(BenchTest):
  name = "RCGEMV"
  executable = "tests/benchs/bench_rcgemv"

class BenchRDGEMMTest(BenchTest):
  name = "RDGEMM"
  executable = "tests/benchs/bench_rdgemm"

class BenchRZGEMMTest(BenchTest):
  name = "RZGEMM"
  executable = "tests/benchs/bench_rzgemm"

class BenchRSGEMMTest(BenchTest):
  name = "RSGEMM"
  executable = "tests/benchs/bench_rsgemm"

class BenchRCGEMMTest(BenchTest):
  name = "RCGEMM"
  executable = "tests/benchs/bench_rcgemm"

class BenchDGEMVTest(BenchTest):
  name = "DGEMV"
  executable = "tests/benchs/bench_dgemv"

class BenchZGEMVTest(BenchTest):
  name = "ZGEMV"
  executable = "tests/benchs/bench_zgemv"

class BenchSGEMVTest(BenchTest):
  name = "SGEMV"
  executable = "tests/benchs/bench_sgemv"

class BenchCGEMVTest(BenchTest):
  name = "CGEMV"
  executable = "tests/benchs/bench_cgemv"

class BenchDGEMMTest(BenchTest):
  name = "DGEMM"
  executable = "tests/benchs/bench_dgemm"

class BenchZGEMMTest(BenchTest):
  name = "ZGEMM"
  executable = "tests/benchs/bench_zgemm"

class BenchSGEMMTest(BenchTest):
  name = "SGEMM"
  executable = "tests/benchs/bench_sgemm"

class BenchCGEMMTest(BenchTest):
  name = "CGEMM"
  executable = "tests/benchs/bench_cgemm"

all_benchs = {"bench_camax": (BenchCAMAXTest, [""]),\
              "bench_camaxm": (BenchCAMAXMTest, [""]),\
              "bench_cciconv": (BenchCCICONVTest, [""]),\
              "bench_cdotc": (BenchCDOTCTest, [""]),\
              "bench_cdotu": (BenchCDOTUTest, [""]),\
              "bench_cgemm": (BenchCGEMVTest, [""]),\
              "bench_cgemv": (BenchCGEMVTest, [""]),\
              "bench_ciciadd": (BenchCICIADDTest, [""]),\
              "bench_damax": (BenchDAMAXTest, [""]),\
              "bench_damaxm": (BenchDAMAXMTest, [""]),\
              "bench_dasum": (BenchDASUMTest, [""]),\
              "bench_ddiconv": (BenchDDICONVTest, [""]),\
              "bench_ddot": (BenchDDOTTest, [""]),\
              "bench_dgemm": (BenchDGEMVTest, [""]),\
              "bench_dgemv": (BenchDGEMVTest, [""]),\
              "bench_didiadd": (BenchDIDIADDTest, [""]),\
              "bench_dnrm2": (BenchDNRM2Test, [""]),\
              "bench_dzasum": (BenchDZASUMTest, [""]),\
              "bench_dznrm2": (BenchDZNRM2Test, [""]),\
              "bench_icamax": (BenchICAMAXTest, [""]),\
              "bench_idamax": (BenchIDAMAXTest, [""]),\
              "bench_isamax": (BenchISAMAXTest, [""]),\
              "bench_izamax": (BenchIZAMAXTest, [""]),\
              "bench_samax": (BenchSAMAXTest, [""]),\
              "bench_samaxm": (BenchSAMAXMTest, [""]),\
              "bench_sasum": (BenchSASUMTest, [""]),\
              "bench_scasum": (BenchSCASUMTest, [""]),\
              "bench_scnrm2": (BenchSCNRM2Test, [""]),\
              "bench_sdot": (BenchSDOTTest, [""]),\
              "bench_sgemm": (BenchSGEMVTest, [""]),\
              "bench_sgemv": (BenchSGEMVTest, [""]),\
              "bench_sisiadd": (BenchSISIADDTest, [""]),\
              "bench_snrm2": (BenchSNRM2Test, [""]),\
              "bench_ssiconv": (BenchSSICONVTest, [""]),\
              "bench_zamax": (BenchZAMAXTest, [""]),\
              "bench_zamaxm": (BenchZAMAXMTest, [""]),\
              "bench_zdotc": (BenchZDOTCTest, [""]),\
              "bench_zdotu": (BenchZDOTUTest, [""]),\
              "bench_zgemm": (BenchZGEMVTest, [""]),\
              "bench_zgemv": (BenchZGEMVTest, [""]),\
              "bench_ziziadd": (BenchZIZIADDTest, [""]),\
              "bench_zziconv": (BenchZZICONVTest, [""]),\
             }

for i in range(1, terminal.get_simaxfold() + 1):
  if i == 1:
    i = 0
    flagss = ["--fold {}".format(j) for j in range(2, terminal.get_simaxfold() + 1)]
  else:
    flagss = ["--fold {}".format(i)]
  all_benchs.update({"bench_rssum_fold_{}".format(i): (BenchRSSUMTest, flagss),\
                     "bench_rsasum_fold_{}".format(i): (BenchRSASUMTest, flagss),\
                     "bench_rsnrm2_fold_{}".format(i): (BenchRSNRM2Test, flagss),\
                     "bench_rsdot_fold_{}".format(i): (BenchRSDOTTest, flagss),\
                     "bench_rsgemv_fold_{}".format(i): (BenchRSGEMVTest, flagss),\
                     "bench_rsgemv_TransA_fold_{}".format(i): (BenchRSGEMVTest, ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rsgemv_AvgTransA_fold_{}".format(i): (BenchRSGEMVTest, flagss + \
                                                                                   ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rsgemm_fold_{}".format(i): (BenchRSGEMMTest, flagss),\
                     "bench_rsgemm_TransB_fold_{}".format(i): (BenchRSGEMMTest, ["--TransB Trans " + flags for flags in flagss]),\
                     "bench_rsgemm_TransA_fold_{}".format(i): (BenchRSGEMMTest, ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rsgemm_TransA_TransB_fold_{}".format(i): (BenchRSGEMMTest, ["--TransA Trans --TransB Trans " + flags for flags in flagss]),\
                     "bench_rsgemm_AvgTransA_AvgTransB_fold_{}".format(i): (BenchRSGEMMTest, flagss + \
                                                                                             ["--TransB Trans " + flags for flags in flagss] + \
                                                                                             ["--TransA Trans " + flags for flags in flagss] + \
                                                                                             ["--TransA Trans --TransB Trans " + flags for flags in flagss]),\
                     "bench_rcsum_fold_{}".format(i): (BenchRCSUMTest, flagss),\
                     "bench_rscasum_fold_{}".format(i): (BenchRSCASUMTest, flagss),\
                     "bench_rscnrm2_fold_{}".format(i): (BenchRSCNRM2Test, flagss),\
                     "bench_rcdotu_fold_{}".format(i): (BenchRCDOTUTest, flagss),\
                     "bench_rcdotc_fold_{}".format(i): (BenchRCDOTCTest, flagss),\
                     "bench_rcgemv_fold_{}".format(i): (BenchRCGEMVTest, flagss),\
                     "bench_rcgemv_TransA_fold_{}".format(i): (BenchRCGEMVTest, ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rcgemv_AvgTransA_fold_{}".format(i): (BenchRCGEMVTest, flagss + \
                                                                                   ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rcgemm_fold_{}".format(i): (BenchRCGEMMTest, flagss),\
                     "bench_rcgemm_TransB_fold_{}".format(i): (BenchRCGEMMTest, ["--TransB Trans " + flags for flags in flagss]),\
                     "bench_rcgemm_TransA_fold_{}".format(i): (BenchRCGEMMTest, ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rcgemm_TransA_TransB_fold_{}".format(i): (BenchRCGEMMTest, ["--TransA Trans --TransB Trans " + flags for flags in flagss]),\
                     "bench_rcgemm_AvgTransA_AvgTransB_fold_{}".format(i): (BenchRCGEMMTest, flagss + \
                                                                                             ["--TransB Trans " + flags for flags in flagss] + \
                                                                                             ["--TransA Trans " + flags for flags in flagss] + \
                                                                                             ["--TransA Trans --TransB Trans " + flags for flags in flagss]),\
                    })

for i in range(1, terminal.get_dimaxfold() + 1):
  if i == 1:
    i = 0
    flagss = ["--fold {}".format(j) for j in range(2, terminal.get_dimaxfold() + 1)]
  else:
    flagss = ["--fold {}".format(i)]
  all_benchs.update({"bench_rdsum_fold_{}".format(i): (BenchRDSUMTest, flagss),\
                     "bench_rdasum_fold_{}".format(i): (BenchRDASUMTest, flagss),\
                     "bench_rdnrm2_fold_{}".format(i): (BenchRDNRM2Test, flagss),\
                     "bench_rddot_fold_{}".format(i): (BenchRDDOTTest, flagss),\
                     "bench_rdgemv_fold_{}".format(i): (BenchRDGEMVTest, flagss),\
                     "bench_rdgemv_TransA_fold_{}".format(i): (BenchRDGEMVTest, ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rdgemv_AvgTransA_fold_{}".format(i): (BenchRDGEMVTest, flagss + \
                                                                                   ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rdgemm_fold_{}".format(i): (BenchRDGEMMTest, flagss),\
                     "bench_rdgemm_TransB_fold_{}".format(i): (BenchRDGEMMTest, ["--TransB Trans " + flags for flags in flagss]),\
                     "bench_rdgemm_TransA_fold_{}".format(i): (BenchRDGEMMTest, ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rdgemm_TransA_TransB_fold_{}".format(i): (BenchRDGEMMTest, ["--TransA Trans --TransB Trans " + flags for flags in flagss]),\
                     "bench_rdgemm_AvgTransA_AvgTransB_fold_{}".format(i): (BenchRDGEMMTest, flagss + \
                                                                                             ["--TransB Trans " + flags for flags in flagss] + \
                                                                                             ["--TransA Trans " + flags for flags in flagss] + \
                                                                                             ["--TransA Trans --TransB Trans " + flags for flags in flagss]),\
                     "bench_rdgemm_TransB_fold_{}".format(i): (BenchRDGEMMTest, ["--TransB Trans " + flags for flags in flagss]),\
                     "bench_rdgemm_TransA_fold_{}".format(i): (BenchRDGEMMTest, ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rdgemm_TransA_TransB_fold_{}".format(i): (BenchRDGEMMTest, ["--TransA Trans --TransB Trans " + flags for flags in flagss]),\
                     "bench_rzsum_fold_{}".format(i): (BenchRZSUMTest, flagss),\
                     "bench_rdzasum_fold_{}".format(i): (BenchRDZASUMTest, flagss),\
                     "bench_rdznrm2_fold_{}".format(i): (BenchRDZNRM2Test, flagss),\
                     "bench_rzdotc_fold_{}".format(i): (BenchRZDOTCTest, flagss),\
                     "bench_rzdotu_fold_{}".format(i): (BenchRZDOTUTest, flagss),\
                     "bench_rzgemv_fold_{}".format(i): (BenchRZGEMVTest, flagss),\
                     "bench_rzgemv_TransA_fold_{}".format(i): (BenchRZGEMVTest, ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rzgemv_AvgTransA_fold_{}".format(i): (BenchRZGEMVTest, flagss + \
                                                                                   ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rzgemm_fold_{}".format(i): (BenchRZGEMMTest, flagss),\
                     "bench_rzgemm_TransB_fold_{}".format(i): (BenchRZGEMMTest, ["--TransB Trans " + flags for flags in flagss]),\
                     "bench_rzgemm_TransA_fold_{}".format(i): (BenchRZGEMMTest, ["--TransA Trans " + flags for flags in flagss]),\
                     "bench_rzgemm_TransA_TransB_fold_{}".format(i): (BenchRZGEMMTest, ["--TransA Trans --TransB Trans " + flags for flags in flagss]),\
                     "bench_rzgemm_AvgTransA_AvgTransB_fold_{}".format(i): (BenchRZGEMMTest, flagss + \
                                                                                             ["--TransB Trans " + flags for flags in flagss] + \
                                                                                             ["--TransA Trans " + flags for flags in flagss] + \
                                                                                             ["--TransA Trans --TransB Trans " + flags for flags in flagss]),\
                    })
