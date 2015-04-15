import itertools
import os

from scripts import terminal
from tests.checks import checks
import config

fill_types=["rand", "2*rand-1", "rand+(rand-1)", "normal", "sine", "small+grow*big", "small+rand*big", "rand_cond"]
data_types=["double", "float", "double_complex", "float_complex"]

d_data = ["data/d_2_rand-1_N4095.dat", "data/d_normal_N4095.dat", "data/d_rand+(rand-1)_N4095.dat", "data/d_rand_N4095.dat", "data/d_rand_cond_N4095.dat", "data/d_sine_N4095.dat", "data/d_small+grow_big_N4095.dat", "data/d_small+rand_big_N4095.dat"]
s_data = ["data/s_2_rand-1_N4095.dat", "data/s_normal_N4095.dat", "data/s_rand+(rand-1)_N4095.dat", "data/s_rand_N4095.dat", "data/s_rand_cond_N4095.dat", "data/s_sine_N4095.dat", "data/s_small+grow_big_N4095.dat", "data/s_small+rand_big_N4095.dat"]
z_data = ["data/z_2_rand-1_N4095.dat", "data/z_normal_N4095.dat", "data/z_rand+(rand-1)_N4095.dat", "data/z_rand_N4095.dat", "data/z_rand_cond_N4095.dat", "data/z_sine_N4095.dat", "data/z_small+grow_big_N4095.dat", "data/z_small+rand_big_N4095.dat"]
c_data = ["data/c_2_rand-1_N4095.dat","data/c_normal_N4095.dat", "data/c_rand+(rand-1)_N4095.dat", "data/c_rand_N4095.dat", "data/c_rand_cond_N4095.dat", "data/c_sine_N4095.dat", "data/c_small+grow_big_N4095.dat", "data/c_small+rand_big_N4095.dat"]


test_data = terminal.make("tests/common/test_data")
for (data_type, fill_type) in itertools.product(data_types, fill_types):
  terminal.call("cd {}; {} {}".format(os.path.join(terminal.top, "tests/checks/data"), test_data, terminal.flags(["N", "d", "f"], [4095, data_type, fill_type])))

d_checks = [checks.ValidateExternalRDSUMTest,\
            checks.ValidateExternalRDASUMTest,\
            checks.ValidateExternalRDNRM2Test,\
            checks.ValidateExternalRDDOTTest]
for data in d_data:
  for check in d_checks:
    check = check()
    check.setup(flags=terminal.flags(["i", "r"], [data, ""]))
    config.run(check.get_command_list())

z_checks = [checks.ValidateExternalRZSUMTest,\
            checks.ValidateExternalRDZASUMTest,\
            checks.ValidateExternalRDZNRM2Test,\
            checks.ValidateExternalRZDOTUTest,\
            checks.ValidateExternalRZDOTCTest]

for data in z_data:
  for check in z_checks:
    check = check()
    check.setup(flags=terminal.flags(["i", "r"], [data, ""]))
    config.run(check.get_command_list())

s_checks = [checks.ValidateExternalRSSUMTest,\
            checks.ValidateExternalRSASUMTest,\
            checks.ValidateExternalRSNRM2Test,\
            checks.ValidateExternalRSDOTTest]

for data in s_data:
  for check in s_checks:
    check = check()
    check.setup(flags=terminal.flags(["i", "r"], [data, ""]))
    config.run(check.get_command_list())

c_checks = [checks.ValidateExternalRCSUMTest,\
            checks.ValidateExternalRSCASUMTest,\
            checks.ValidateExternalRSCNRM2Test,\
            checks.ValidateExternalRCDOTUTest,\
            checks.ValidateExternalRCDOTCTest]

for data in c_data:
  for check in c_checks:
    check = check()
    check.setup(flags=terminal.flags(["i", "r"], [data, ""]))
    config.run(check.get_command_list())
