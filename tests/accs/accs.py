import copy
import itertools

import scripts.terminal as terminal
import tests.harness.harness as harness

class AccSuite(harness.MetricSuite):
  pass

class AccTest(harness.MetricTest):
  def parse_output(self):
    if self.attribute == "max_ratio":
      return output["max_ratio"]
    elif self.attribute == "max_ratio(e)":
      return output["max_ratio"]/output["e"]
    elif self.attribute == "med_ratio":
      return output["med_ratio"]
    elif self.attribute == "med_ratio(e)":
      return output["med_ratio"]/output["e"]
    elif self.attribute == "min_ratio":
      return output["min_ratio"]
    elif self.attribute == "min_ratio(e)":
      return output["min_ratio"]/output["e"]

class AccRCSUMTest(AccTest):
  name = "RCSUM"
  executable = "tests/accs/acc_rcsum"

class AccRDSUMTest(AccTest):
  name = "RDSUM"
  executable = "tests/accs/acc_rdsum"

class AccRSSUMTest(AccTest):
  name = "RSSUM"
  executable = "tests/accs/acc_rssum"

class AccRZSUMTest(AccTest):
  name = "RZSUM"
  executable = "tests/accs/acc_rzsum"

class AccDDICONVTest(AccTest):
  name = "DDICONV"
  executable = "tests/accs/acc_ddbconv"

class AccZZICONVTest(AccTest):
  name = "ZZICONV"
  executable = "tests/accs/acc_zziconv"

class AccSSICONVTest(AccTest):
  name = "SSICONV"
  executable = "tests/accs/acc_ssbconv"

class AccCCICONVTest(AccTest):
  name = "CCICONV"
  executable = "tests/accs/acc_cciconv"
