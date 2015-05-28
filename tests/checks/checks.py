import copy
import itertools

import scripts.terminal as terminal
import tests.harness.harness as harness
import config

class CheckSuite(harness.Suite):

  def __init__(self):
    self.checks = []
    self.check_rows = []
    self.args = []
    self.params = []

  def add_checks(self, checks, params, ranges):
    for args in itertools.product(*ranges):
      check_row = copy.deepcopy(checks)
      self.check_rows.append(check_row)
      self.checks += check_row
      self.params.append(params)
      self.args.append(args)

  def setup(self, verbose = "false", **kwargs):
    self.verbose = verbose;
    for check_row, params, args in zip(self.check_rows, self.params, self.args):
      for check in check_row:
        check.setup(flags = terminal.flags(params, args), verbose=verbose, **kwargs)

  def get_command_list(self):
    command_list = []
    for check in self.checks:
      command_list += check.get_command_list()
    return command_list

  def parse_output_list(self, output_list):
    for check in self.checks:
      check.parse_output_list(output_list[:len(check.get_command_list())])
      output_list = output_list[len(check.get_command_list()):]

  def get_header(self):
    return ["Check", "Res"]

  def get_align(self):
    return ["l", "c"]

  def get_dtype(self):
    return ["t", "t"]

  def get_cols_width(self, max_width):
    return [max_width - 1 - 1 - 4 - 1, 4]

  def get_rows(self):
    passed = 0
    failed = 0
    na = 0
    rows = []
    for check in self.checks:
      if check.get_result() == 0:
        if self.verbose == "true":
          rows.append([check.get_name(), "Pass"])
        passed += 1
      elif check.get_result() == 125:
        rows.append([check.get_name() + "\n" + check.get_output(), "N/A"])
        na += 1
      else:
        rows.append([check.get_name() + "\n" + check.get_output(), "Fail"])
        failed += 1
    emoticon = ":("
    if passed == len(self.checks):
      emoticon = ":D"
    rows.append(["Passed: {0}/{3} Failed: {1}/{3} N/A: {2}/{3}".format(passed, failed, na, len(self.checks)), emoticon])
    return rows

  def get_output(self):
    return "\n".join(self.get_rows)

  def get_result(self):
    return "\n".join(self.get_rows)

class CheckTest(harness.ExecutableTest):

  def get_name(self):
    """
    return the name of the test
    """
    return self.name

  def get_command_list(self):
    """
    return a list of commands that constitute the test to be run on the
    target architecture
    """
    return ["{} {} {}".format(self.executable_output, self.base_flags, self.flags), "{} {} {}".format(self.executable_output, self.base_flags, self.flags + " -p")]

  def parse_output_list(self, output_list):
    """
    parse the output of the command set. The output will be given as a list of
    (return code, output)
    """
    assert len(output_list) == 2, "ReproBLAS error: unexpected test output"
    self.output = output_list[0][1]
    self.result = output_list[0][0]
    self.name = output_list[1][1]

  def get_output(self):
    """
    return all relevant output (mostly for debugging)
    """
    return self.output

  def get_result(self):
    """
    return test result
    """
    return self.result

class ValidateInternalUFPTest(CheckTest):
  executable = "tests/checks/validate_internal_ufp"
  name = "validate_internal_ufp"

class ValidateInternalUFPFTest(CheckTest):
  executable = "tests/checks/validate_internal_ufpf"
  name = "validate_internal_ufpf"

class ValidateInternalDAMAXTest(CheckTest):
  executable = "tests/checks/validate_internal_damax"
  name = "validate_internal_damax"

class ValidateInternalZAMAXTest(CheckTest):
  executable = "tests/checks/validate_internal_zamax"
  name = "validate_internal_zamax"

class ValidateInternalSAMAXTest(CheckTest):
  executable = "tests/checks/validate_internal_samax"
  name = "validate_internal_samax"

class ValidateInternalCAMAXTest(CheckTest):
  executable = "tests/checks/validate_internal_camax"
  name = "validate_internal_camax"

class ValidateInternalRDSUMTest(CheckTest):
  base_flags = "-w rdsum"
  executable = "tests/checks/validate_internal_daugsum"
  name = "validate_internal_rdsum"

class ValidateInternalRDASUMTest(CheckTest):
  base_flags = "-w rdasum"
  executable = "tests/checks/validate_internal_daugsum"
  name = "validate_internal_rdasum"

class ValidateInternalRDNRM2Test(CheckTest):
  base_flags = "-w rdnrm2"
  executable = "tests/checks/validate_internal_daugsum"
  name = "validate_internal_rdrnm2"

class ValidateInternalRDDOTTest(CheckTest):
  base_flags = "-w rddot"
  executable = "tests/checks/validate_internal_daugsum"
  name = "validate_internal_rddot"

class ValidateInternalDIDIADDTest(CheckTest):
  base_flags = "-w didiadd"
  executable = "tests/checks/validate_internal_daugsum"
  name = "validate_internal_didiadd"

class ValidateInternalDIDADDTest(CheckTest):
  base_flags = "-w didadd"
  executable = "tests/checks/validate_internal_daugsum"
  name = "validate_internal_didadd"

class ValidateInternalDIDDEPOSITTest(CheckTest):
  base_flags = "-w diddeposit"
  executable = "tests/checks/validate_internal_daugsum"
  name = "validate_internal_diddeposit"

class ValidateInternalRZSUMTest(CheckTest):
  base_flags = "-w rzsum"
  executable = "tests/checks/validate_internal_zaugsum"
  name = "validate_internal_rzsum"

class ValidateInternalRDZASUMTest(CheckTest):
  base_flags = "-w rdzasum"
  executable = "tests/checks/validate_internal_zaugsum"
  name = "validate_internal_rdzasum"

class ValidateInternalRDZNRM2Test(CheckTest):
  base_flags = "-w rdznrm2"
  executable = "tests/checks/validate_internal_zaugsum"
  name = "validate_internal_rdzrnm2"

class ValidateInternalRZDOTUTest(CheckTest):
  base_flags = "-w rzdotu"
  executable = "tests/checks/validate_internal_zaugsum"
  name = "validate_internal_rzdotu"

class ValidateInternalRZDOTCTest(CheckTest):
  base_flags = "-w rzdotc"
  executable = "tests/checks/validate_internal_zaugsum"
  name = "validate_internal_rzdotc"

class ValidateInternalZIZIADDTest(CheckTest):
  base_flags = "-w ziziadd"
  executable = "tests/checks/validate_internal_zaugsum"
  name = "validate_internal_ziziadd"

class ValidateInternalZIZADDTest(CheckTest):
  base_flags = "-w zizadd"
  executable = "tests/checks/validate_internal_zaugsum"
  name = "validate_internal_zizadd"

class ValidateInternalZIZDEPOSITTest(CheckTest):
  base_flags = "-w zizdeposit"
  executable = "tests/checks/validate_internal_zaugsum"
  name = "validate_internal_zizdeposit"

class ValidateInternalRSSUMTest(CheckTest):
  base_flags = "-w rssum"
  executable = "tests/checks/validate_internal_saugsum"
  name = "validate_internal_rssum"

class ValidateInternalRSASUMTest(CheckTest):
  base_flags = "-w rsasum"
  executable = "tests/checks/validate_internal_saugsum"
  name = "validate_internal_rsasum"

class ValidateInternalRSNRM2Test(CheckTest):
  base_flags = "-w rsnrm2"
  executable = "tests/checks/validate_internal_saugsum"
  name = "validate_internal_rsrnm2"

class ValidateInternalRSDOTTest(CheckTest):
  base_flags = "-w rsdot"
  executable = "tests/checks/validate_internal_saugsum"
  name = "validate_internal_rsdot"

class ValidateInternalSISIADDTest(CheckTest):
  base_flags = "-w sisiadd"
  executable = "tests/checks/validate_internal_saugsum"
  name = "validate_internal_sisiadd"

class ValidateInternalSISADDTest(CheckTest):
  base_flags = "-w sisadd"
  executable = "tests/checks/validate_internal_saugsum"
  name = "validate_internal_sisadd"

class ValidateInternalSISDEPOSITTest(CheckTest):
  base_flags = "-w sisdeposit"
  executable = "tests/checks/validate_internal_saugsum"
  name = "validate_internal_sisdeposit"

class ValidateInternalRCSUMTest(CheckTest):
  base_flags = "-w rcsum"
  executable = "tests/checks/validate_internal_caugsum"
  name = "validate_internal_rcsum"

class ValidateInternalRSCASUMTest(CheckTest):
  base_flags = "-w rscasum"
  executable = "tests/checks/validate_internal_caugsum"
  name = "validate_internal_rscasum"

class ValidateInternalRSCNRM2Test(CheckTest):
  base_flags = "-w rscnrm2"
  executable = "tests/checks/validate_internal_caugsum"
  name = "validate_internal_rscrnm2"

class ValidateInternalRCDOTUTest(CheckTest):
  base_flags = "-w rcdotu"
  executable = "tests/checks/validate_internal_caugsum"
  name = "validate_internal_rcdotu"

class ValidateInternalRCDOTCTest(CheckTest):
  base_flags = "-w rcdotc"
  executable = "tests/checks/validate_internal_caugsum"
  name = "validate_internal_rcdotc"

class ValidateInternalCICIADDTest(CheckTest):
  base_flags = "-w ciciadd"
  executable = "tests/checks/validate_internal_caugsum"
  name = "validate_internal_ciciadd"

class ValidateInternalCICADDTest(CheckTest):
  base_flags = "-w cicadd"
  executable = "tests/checks/validate_internal_caugsum"
  name = "validate_internal_cicadd"

class ValidateInternalCICDEPOSITTest(CheckTest):
  base_flags = "-w cicdeposit"
  executable = "tests/checks/validate_internal_caugsum"
  name = "validate_internal_cicdeposit"

class ValidateInternalRZBLAS1Test(CheckTest):
  executable = "tests/checks/validate_internal_rzblas1"
  name = "validate_internal_rzblas1"

class ValidateInternalRSBLAS1Test(CheckTest):
  executable = "tests/checks/validate_internal_rsblas1"
  name = "validate_internal_rsblas1"

class ValidateInternalRCBLAS1Test(CheckTest):
  executable = "tests/checks/validate_internal_rcblas1"
  name = "validate_internal_rcblas1"

class ValidateExternalRDSUMTest(CheckTest):
  base_flags = "-w rdsum"
  executable = "tests/checks/validate_external_rdblas1"
  name = "validate_external_rdsum"

class ValidateExternalRDASUMTest(CheckTest):
  base_flags = "-w rdasum"
  executable = "tests/checks/validate_external_rdblas1"
  name = "validate_external_rdasum"

class ValidateExternalRDNRM2Test(CheckTest):
  base_flags = "-w rdnrm2"
  executable = "tests/checks/validate_external_rdblas1"
  name = "validate_external_rdnrm2"

class ValidateExternalRDDOTTest(CheckTest):
  base_flags = "-w rddot"
  executable = "tests/checks/validate_external_rdblas1"
  name = "validate_external_rddot"

class ValidateExternalRZSUMTest(CheckTest):
  base_flags = "-w rzsum"
  executable = "tests/checks/validate_external_rzblas1"
  name = "validate_external_rzsum"

class ValidateExternalRDZASUMTest(CheckTest):
  base_flags = "-w rdzasum"
  executable = "tests/checks/validate_external_rzblas1"
  name = "validate_external_rdzasum"

class ValidateExternalRDZNRM2Test(CheckTest):
  base_flags = "-w rdznrm2"
  executable = "tests/checks/validate_external_rzblas1"
  name = "validate_external_rdznrm2"

class ValidateExternalRZDOTUTest(CheckTest):
  base_flags = "-w rzdotu"
  executable = "tests/checks/validate_external_rzblas1"
  name = "validate_external_rzdotu"

class ValidateExternalRZDOTCTest(CheckTest):
  base_flags = "-w rzdotc"
  executable = "tests/checks/validate_external_rzblas1"
  name = "validate_external_rzdotc"

class ValidateExternalRSSUMTest(CheckTest):
  base_flags = "-w rssum"
  executable = "tests/checks/validate_external_rsblas1"
  name = "validate_external_rssum"

class ValidateExternalRSASUMTest(CheckTest):
  base_flags = "-w rsasum"
  executable = "tests/checks/validate_external_rsblas1"
  name = "validate_external_rsasum"

class ValidateExternalRSNRM2Test(CheckTest):
  base_flags = "-w rsnrm2"
  executable = "tests/checks/validate_external_rsblas1"
  name = "validate_external_rsnrm2"

class ValidateExternalRSDOTTest(CheckTest):
  base_flags = "-w rsdot"
  executable = "tests/checks/validate_external_rsblas1"
  name = "validate_external_rsdot"

class ValidateExternalRCSUMTest(CheckTest):
  base_flags = "-w rcsum"
  executable = "tests/checks/validate_external_rcblas1"
  name = "validate_external_rcsum"

class ValidateExternalRSCASUMTest(CheckTest):
  base_flags = "-w rscasum"
  executable = "tests/checks/validate_external_rcblas1"
  name = "validate_external_rscasum"

class ValidateExternalRSCNRM2Test(CheckTest):
  base_flags = "-w rscnrm2"
  executable = "tests/checks/validate_external_rcblas1"
  name = "validate_external_rscnrm2"

class ValidateExternalRCDOTUTest(CheckTest):
  base_flags = "-w rcdotu"
  executable = "tests/checks/validate_external_rcblas1"
  name = "validate_external_rcdotu"

class ValidateExternalRCDOTCTest(CheckTest):
  base_flags = "-w rcdotc"
  executable = "tests/checks/validate_external_rcblas1"
  name = "validate_external_rcdotc"

class VerifyRDSUMTest(CheckTest):
  base_flags = "-w rdsum"
  executable = "tests/checks/verify_daugsum"
  name = "verify_rdsum"

class VerifyRDASUMTest(CheckTest):
  base_flags = "-w rdasum"
  executable = "tests/checks/verify_daugsum"
  name = "verify_rdasum"

class VerifyRDNRM2Test(CheckTest):
  base_flags = "-w rdnrm2"
  executable = "tests/checks/verify_daugsum"
  name = "verify_rdnrm2"

class VerifyRDDOTTest(CheckTest):
  base_flags = "-w rddot"
  executable = "tests/checks/verify_daugsum"
  name = "verify_rddot"

class VerifyDIDIADDTest(CheckTest):
  base_flags = "-w didiadd"
  executable = "tests/checks/verify_daugsum"
  name = "verify_didiadd"

class VerifyDIDADDTest(CheckTest):
  base_flags = "-w didadd"
  executable = "tests/checks/verify_daugsum"
  name = "verify_didadd"

class VerifyDIDDEPOSITTest(CheckTest):
  base_flags = "-w diddeposit"
  executable = "tests/checks/verify_daugsum"
  name = "verify_diddeposit"

class VerifyRZSUMTest(CheckTest):
  base_flags = "-w rzsum"
  executable = "tests/checks/verify_zaugsum"
  name = "verify_rzsum"

class VerifyRDZASUMTest(CheckTest):
  base_flags = "-w rdzasum"
  executable = "tests/checks/verify_zaugsum"
  name = "verify_rdzasum"

class VerifyRDZNRM2Test(CheckTest):
  base_flags = "-w rdznrm2"
  executable = "tests/checks/verify_zaugsum"
  name = "verify_rdznrm2"

class VerifyRZDOTUTest(CheckTest):
  base_flags = "-w rzdotu"
  executable = "tests/checks/verify_zaugsum"
  name = "verify_rzdotu"

class VerifyRZDOTCTest(CheckTest):
  base_flags = "-w rzdotc"
  executable = "tests/checks/verify_zaugsum"
  name = "verify_rzdotc"

class VerifyZIZIADDTest(CheckTest):
  base_flags = "-w ziziadd"
  executable = "tests/checks/verify_zaugsum"
  name = "verify_ziziadd"

class VerifyZIZADDTest(CheckTest):
  base_flags = "-w zizadd"
  executable = "tests/checks/verify_zaugsum"
  name = "verify_zizadd"

class VerifyZIZDEPOSITTest(CheckTest):
  base_flags = "-w zizdeposit"
  executable = "tests/checks/verify_zaugsum"
  name = "verify_zizdeposit"

class VerifyRSSUMTest(CheckTest):
  base_flags = "-w rssum"
  executable = "tests/checks/verify_saugsum"
  name = "verify_rssum"

class VerifyRSASUMTest(CheckTest):
  base_flags = "-w rsasum"
  executable = "tests/checks/verify_saugsum"
  name = "verify_rsasum"

class VerifyRSNRM2Test(CheckTest):
  base_flags = "-w rsnrm2"
  executable = "tests/checks/verify_saugsum"
  name = "verify_rsnrm2"

class VerifyRSDOTTest(CheckTest):
  base_flags = "-w rsdot"
  executable = "tests/checks/verify_saugsum"
  name = "verify_rsdot"

class VerifySISIADDTest(CheckTest):
  base_flags = "-w sisiadd"
  executable = "tests/checks/verify_saugsum"
  name = "verify_sisiadd"

class VerifySISADDTest(CheckTest):
  base_flags = "-w sisadd"
  executable = "tests/checks/verify_saugsum"
  name = "verify_sisadd"

class VerifySISDEPOSITTest(CheckTest):
  base_flags = "-w sisdeposit"
  executable = "tests/checks/verify_saugsum"
  name = "verify_sisdeposit"

class VerifyRCSUMTest(CheckTest):
  base_flags = "-w rcsum"
  executable = "tests/checks/verify_caugsum"
  name = "verify_rcsum"

class VerifyRSCASUMTest(CheckTest):
  base_flags = "-w rscasum"
  executable = "tests/checks/verify_caugsum"
  name = "verify_rscasum"

class VerifyRSCNRM2Test(CheckTest):
  base_flags = "-w rscnrm2"
  executable = "tests/checks/verify_caugsum"
  name = "verify_rscnrm2"

class VerifyRCDOTUTest(CheckTest):
  base_flags = "-w rcdotu"
  executable = "tests/checks/verify_caugsum"
  name = "verify_rcdotu"

class VerifyRCDOTCTest(CheckTest):
  base_flags = "-w rcdotc"
  executable = "tests/checks/verify_caugsum"
  name = "verify_rcdotc"

class VerifyCICIADDTest(CheckTest):
  base_flags = "-w ciciadd"
  executable = "tests/checks/verify_caugsum"
  name = "verify_ciciadd"

class VerifyCICADDTest(CheckTest):
  base_flags = "-w cicadd"
  executable = "tests/checks/verify_caugsum"
  name = "verify_cicadd"

class VerifyCICDEPOSITTest(CheckTest):
  base_flags = "-w cicdeposit"
  executable = "tests/checks/verify_caugsum"
  name = "verify_cicdeposit"

class VerifyDINDEXTest(CheckTest):
  base_flags = ""
  executable = "tests/checks/verify_dindex"
  name = "verify_dindex"

class VerifySINDEXTest(CheckTest):
  base_flags = ""
  executable = "tests/checks/verify_sindex"
  name = "verify_sindex"

class VerifyDMINDEXTest(CheckTest):
  base_flags = ""
  executable = "tests/checks/verify_dmindex"
  name = "verify_dmindex"

class VerifySMINDEXTest(CheckTest):
  base_flags = ""
  executable = "tests/checks/verify_smindex"
  name = "verify_smindex"

#class VerifyRDGEMVTest(CheckTest):
#  base_flags = ""
#  executable = "tests/checks/verify_rdgemv"
#  name = "verify_rdgemv"

all_checks = {"validate_internal_ufp": ValidateInternalUFPTest,\
              "validate_internal_ufpf": ValidateInternalUFPFTest,\
              "validate_internal_damax": ValidateInternalDAMAXTest,\
              "validate_internal_zamax": ValidateInternalZAMAXTest,\
              "validate_internal_samax": ValidateInternalSAMAXTest,\
              "validate_internal_camax": ValidateInternalCAMAXTest,\
              "validate_internal_rdsum": ValidateInternalRDSUMTest,\
              "validate_internal_rdasum": ValidateInternalRDASUMTest,\
              "validate_internal_rdnrm2": ValidateInternalRDNRM2Test,\
              "validate_internal_rddot": ValidateInternalRDDOTTest,\
              "validate_internal_didiadd": ValidateInternalDIDIADDTest,\
              "validate_internal_didadd": ValidateInternalDIDADDTest,\
              "validate_internal_diddeposit": ValidateInternalDIDDEPOSITTest,\
              "validate_internal_rzsum": ValidateInternalRZSUMTest,\
              "validate_internal_rdzasum": ValidateInternalRDZASUMTest,\
              "validate_internal_rdznrm2": ValidateInternalRDZNRM2Test,\
              "validate_internal_rzdotu": ValidateInternalRZDOTUTest,\
              "validate_internal_rzdotc": ValidateInternalRZDOTCTest,\
              "validate_internal_ziziadd": ValidateInternalZIZIADDTest,\
              "validate_internal_zizadd": ValidateInternalZIZADDTest,\
              "validate_internal_zizdeposit": ValidateInternalZIZDEPOSITTest,\
              "validate_internal_rssum": ValidateInternalRSSUMTest,\
              "validate_internal_rsasum": ValidateInternalRSASUMTest,\
              "validate_internal_rsnrm2": ValidateInternalRSNRM2Test,\
              "validate_internal_rsdot": ValidateInternalRSDOTTest,\
              "validate_internal_sisiadd": ValidateInternalSISIADDTest,\
              "validate_internal_sisadd": ValidateInternalSISADDTest,\
              "validate_internal_sisdeposit": ValidateInternalSISDEPOSITTest,\
              "validate_internal_rcsum": ValidateInternalRCSUMTest,\
              "validate_internal_rscasum": ValidateInternalRSCASUMTest,\
              "validate_internal_rscnrm2": ValidateInternalRSCNRM2Test,\
              "validate_internal_rcdotu": ValidateInternalRCDOTUTest,\
              "validate_internal_rcdotc": ValidateInternalRCDOTCTest,\
              "validate_internal_ciciadd": ValidateInternalCICIADDTest,\
              "validate_internal_cicadd": ValidateInternalCICADDTest,\
              "validate_internal_cicdeposit": ValidateInternalCICDEPOSITTest,\
              "validate_external_rdsum": ValidateExternalRDSUMTest,\
              "validate_external_rdasum": ValidateExternalRDASUMTest,\
              "validate_external_rdnrm2": ValidateExternalRDNRM2Test,\
              "validate_external_rddot": ValidateExternalRDDOTTest,\
              "validate_external_rzsum": ValidateExternalRZSUMTest,\
              "validate_external_rdzasum": ValidateExternalRDZASUMTest,\
              "validate_external_rdznrm2": ValidateExternalRDZNRM2Test,\
              "validate_external_rzdotu": ValidateExternalRZDOTUTest,\
              "validate_external_rzdotc": ValidateExternalRZDOTCTest,\
              "validate_external_rssum": ValidateExternalRSSUMTest,\
              "validate_external_rsasum": ValidateExternalRSASUMTest,\
              "validate_external_rsnrm2": ValidateExternalRSNRM2Test,\
              "validate_external_rsdot": ValidateExternalRSDOTTest,\
              "validate_external_rcsum": ValidateExternalRCSUMTest,\
              "validate_external_rscasum": ValidateExternalRSCASUMTest,\
              "validate_external_rscnrm2": ValidateExternalRSCNRM2Test,\
              "validate_external_rcdotu": ValidateExternalRCDOTUTest,\
              "validate_external_rcdotc": ValidateExternalRCDOTCTest,\
              "verify_rdsum": VerifyRDSUMTest,\
              "verify_rdasum": VerifyRDASUMTest,\
              "verify_rdnrm2": VerifyRDNRM2Test,\
              "verify_rddot": VerifyRDDOTTest,\
              "verify_didiadd": VerifyDIDIADDTest,\
              "verify_didadd": VerifyDIDADDTest,\
              "verify_diddeposit": VerifyDIDDEPOSITTest,\
              "verify_rzsum": VerifyRZSUMTest,\
              "verify_rdzasum": VerifyRDZASUMTest,\
              "verify_rdznrm2": VerifyRDZNRM2Test,\
              "verify_rzdotu": VerifyRZDOTUTest,\
              "verify_rzdotc": VerifyRZDOTCTest,\
              "verify_ziziadd": VerifyZIZIADDTest,\
              "verify_zizadd": VerifyZIZADDTest,\
              "verify_zizdeposit": VerifyZIZDEPOSITTest,\
              "verify_rssum": VerifyRSSUMTest,\
              "verify_rsasum": VerifyRSASUMTest,\
              "verify_rsnrm2": VerifyRSNRM2Test,\
              "verify_rsdot": VerifyRSDOTTest,\
              "verify_sisiadd": VerifySISIADDTest,\
              "verify_sisadd": VerifySISADDTest,\
              "verify_sisdeposit": VerifySISDEPOSITTest,\
              "verify_rcsum": VerifyRCSUMTest,\
              "verify_rscasum": VerifyRSCASUMTest,\
              "verify_rscnrm2": VerifyRSCNRM2Test,\
              "verify_rcdotu": VerifyRCDOTUTest,\
              "verify_rcdotc": VerifyRCDOTCTest,\
              "verify_ciciadd": VerifyCICIADDTest,\
              "verify_cicadd": VerifyCICADDTest,\
              "verify_cicdeposit": VerifyCICDEPOSITTest,\
              "verify_dindex": VerifyDINDEXTest,\
              "verify_sindex": VerifySINDEXTest,\
              "verify_dmindex": VerifyDMINDEXTest,\
              "verify_smindex": VerifySMINDEXTest}
