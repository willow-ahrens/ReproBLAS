import copy
import itertools

import scripts.terminal as terminal
import tests.harness.harness as harness

class AccSuite(harness.MetricSuite):
  pass

class AccTest(harness.MetricTest):
  def parse_output(self):
    if self.attribute == "ratio":
      self.result = self.output["ratio"]/self.output["trials"]

class AccCDOTCTest(AccTest):
  name = "CDOTC"
  executable = "tests/accs/acc_cdotc"

class AccCDOTUTest(AccTest):
  name = "CDOTU"
  executable = "tests/accs/acc_cdotu"

class AccDASUMTest(AccTest):
  name = "DASUM"
  executable = "tests/accs/acc_dasum"

class AccDDOTTest(AccTest):
  name = "DDOT"
  executable = "tests/accs/acc_ddot"

class AccDNRM2Test(AccTest):
  name = "DNRM2"
  executable = "tests/accs/acc_dnrm2"

class AccDZASUMTest(AccTest):
  name = "DZASUM"
  executable = "tests/accs/acc_dzasum"

class AccDZNRM2Test(AccTest):
  name = "DZNRM2"
  executable = "tests/accs/acc_dznrm2"

class AccRCDOTCTest(AccTest):
  name = "RCDOTC"
  executable = "tests/accs/acc_rcdotc"

class AccRCDOTUTest(AccTest):
  name = "RCDOTU"
  executable = "tests/accs/acc_rcdotu"

class AccRCSUMTest(AccTest):
  name = "RCSUM"
  executable = "tests/accs/acc_rcsum"

class AccRDASUMTest(AccTest):
  name = "RDASUM"
  executable = "tests/accs/acc_rdasum"

class AccRDDOTTest(AccTest):
  name = "RDDOT"
  executable = "tests/accs/acc_rddot"

class AccRDNRM2Test(AccTest):
  name = "RDNRM2"
  executable = "tests/accs/acc_rdnrm2"

class AccRDSUMTest(AccTest):
  name = "RDSUM"
  executable = "tests/accs/acc_rdsum"

class AccRDZASUMTest(AccTest):
  name = "RDZASUM"
  executable = "tests/accs/acc_rdzasum"

class AccRDZNRM2Test(AccTest):
  name = "RDZNRM2"
  executable = "tests/accs/acc_rdznrm2"

class AccRSASUMTest(AccTest):
  name = "RSASUM"
  executable = "tests/accs/acc_rsasum"

class AccRSCASUMTest(AccTest):
  name = "RSCASUM"
  executable = "tests/accs/acc_rscasum"

class AccRSCNRM2Test(AccTest):
  name = "RSCNRM2"
  executable = "tests/accs/acc_rscnrm2"

class AccRSDOTTest(AccTest):
  name = "RSDOT"
  executable = "tests/accs/acc_rsdot"

class AccRSNRM2Test(AccTest):
  name = "RSNRM2"
  executable = "tests/accs/acc_rsnrm2"

class AccRSSUMTest(AccTest):
  name = "RSSUM"
  executable = "tests/accs/acc_rssum"

class AccRZDOTCTest(AccTest):
  name = "RZDOTC"
  executable = "tests/accs/acc_rzdotc"

class AccRZDOTUTest(AccTest):
  name = "RZDOTU"
  executable = "tests/accs/acc_rzdotu"

class AccRZSUMTest(AccTest):
  name = "RZSUM"
  executable = "tests/accs/acc_rzsum"

class AccSASUMTest(AccTest):
  name = "SASUM"
  executable = "tests/accs/acc_sasum"

class AccSCASUMTest(AccTest):
  name = "SCASUM"
  executable = "tests/accs/acc_scasum"

class AccSCNRM2Test(AccTest):
  name = "SCNRM2"
  executable = "tests/accs/acc_scnrm2"

class AccSDOTTest(AccTest):
  name = "SDOT"
  executable = "tests/accs/acc_sdot"

class AccSNRM2Test(AccTest):
  name = "SNRM2"
  executable = "tests/accs/acc_snrm2"

class AccZDOTCTest(AccTest):
  name = "ZDOTC"
  executable = "tests/accs/acc_zdotc"

class AccZDOTUTest(AccTest):
  name = "ZDOTU"
  executable = "tests/accs/acc_zdotu"

class AccDDICONVTest(AccTest):
  name = "DDICONV"
  executable = "tests/accs/acc_ddiconv"

class AccZZICONVTest(AccTest):
  name = "ZZICONV"
  executable = "tests/accs/acc_zziconv"

class AccSSICONVTest(AccTest):
  name = "SSICONV"
  executable = "tests/accs/acc_ssiconv"

class AccCCICONVTest(AccTest):
  name = "CCICONV"
  executable = "tests/accs/acc_cciconv"

class AccRDGEMVTest(AccTest):
  name = "RDGEMV"
  executable = "tests/accs/acc_rdgemv"

class AccDGEMVTest(AccTest):
  name = "DGEMV"
  executable = "tests/accs/acc_dgemv"

class AccPRDGEMVTest(AccTest):
  name = "PRDGEMV"
  executable = "tests/accs/acc_prdgemv"

class AccPRBDGEMVTest(AccTest):
  name = "PRBDGEMV"
  executable = "tests/accs/acc_prbdgemv"

class AccPDGEMVTest(AccTest):
  name = "PDGEMV"
  executable = "tests/accs/acc_pdgemv"

all_accs = {"acc_cdotc": (AccCDOTCTest, ""),\
            "acc_cdotu": (AccCDOTUTest, ""),\
            "acc_dasum": (AccDASUMTest, ""),\
            "acc_ddot": (AccDDOTTest, ""),\
            "acc_dnrm2": (AccDNRM2Test, ""),\
            "acc_dzasum": (AccDZASUMTest, ""),\
            "acc_dznrm2": (AccDZNRM2Test, ""),\
            "acc_sasum": (AccSASUMTest, ""),\
            "acc_scasum": (AccSCASUMTest, ""),\
            "acc_scnrm2": (AccSCNRM2Test, ""),\
            "acc_sdot": (AccSDOTTest, ""),\
            "acc_snrm2": (AccSNRM2Test, ""),\
            "acc_zdotc": (AccZDOTCTest, ""),\
            "acc_zdotu": (AccZDOTUTest, "")}

for i in range(terminal.get_max_fold() + 1):
  all_accs.update({"acc_rcdotc_fold_{}".format(i): (AccRCDOTCTest, "--fold {}".format(i)),\
                   "acc_rcdotu_fold_{}".format(i): (AccRCDOTUTest, "--fold {}".format(i)),\
                   "acc_rcsum_fold_{}".format(i): (AccRCSUMTest, "--fold {}".format(i)),\
                   "acc_rdasum_fold_{}".format(i): (AccRDASUMTest, "--fold {}".format(i)),\
                   "acc_rddot_fold_{}".format(i): (AccRDDOTTest, "--fold {}".format(i)),\
                   "acc_rdnrm2_fold_{}".format(i): (AccRDNRM2Test, "--fold {}".format(i)),\
                   "acc_rdsum_fold_{}".format(i): (AccRDSUMTest, "--fold {}".format(i)),\
                   "acc_rdzasum_fold_{}".format(i): (AccRDZASUMTest, "--fold {}".format(i)),\
                   "acc_rdznrm2_fold_{}".format(i): (AccRDZNRM2Test, "--fold {}".format(i)),\
                   "acc_rsasum_fold_{}".format(i): (AccRSASUMTest, "--fold {}".format(i)),\
                   "acc_rscasum_fold_{}".format(i): (AccRSCASUMTest, "--fold {}".format(i)),\
                   "acc_rscnrm2_fold_{}".format(i): (AccRSCNRM2Test, "--fold {}".format(i)),\
                   "acc_rsdot_fold_{}".format(i): (AccRSDOTTest, "--fold {}".format(i)),\
                   "acc_rsnrm2_fold_{}".format(i): (AccRSNRM2Test, "--fold {}".format(i)),\
                   "acc_rssum_fold_{}".format(i): (AccRSSUMTest, "--fold {}".format(i)),\
                   "acc_rzdotc_fold_{}".format(i): (AccRZDOTCTest, "--fold {}".format(i)),\
                   "acc_rzdotu_fold_{}".format(i): (AccRZDOTUTest, "--fold {}".format(i)),\
                   "acc_rzsum_fold_{}".format(i): (AccRZSUMTest, "--fold {}".format(i))})
