#!/usr/bin/python
import argparse
import copy
import os

import opentuner
from opentuner.api import TuningRunManager
from opentuner.resultsdb.models import Result, TuningRun
from opentuner.search import manipulator, objective
from opentuner import measurement
from opentuner.measurement import inputmanager
import logging

import config
from src.gen import generate
from tests.benchs import benchs
from scripts import terminal

argparser = argparse.ArgumentParser(parents=opentuner.argparsers())
argparser.add_argument('--params', help='parameter space file')
argparser.add_argument('--args', help='output arguments file')
argparser.add_argument("--trials", type=int,default=100, help="number of trials to run for each benchmark")
argparser.add_argument('-v', '--verbose', default="false", type=str, nargs="?", help='verbose if "true"')

log = logging.getLogger(__name__)

class ReproBLASTuner(measurement.MeasurementInterface):
  def __init__(self, args, *pargs, **kwargs):
    super(ReproBLASTuner, self).__init__(project_name="ReproBLAS", program_name=args.benchmark, program_version=config.version, *pargs, **kwargs)
    self.benchmark = args.benchmark
    self.parameter_space = generate.deserialize_parameter_space(args.params)
    self.arguments_file_name = args.args
    if os.path.isfile(self.arguments_file_name):
      self.arguments = generate.deserialize_arguments(self.arguments_file_name)
    else:
      self.arguments = self.parameter_space.get_default_arguments()
    self.parallel_compile = False

  def run(self):
    pass

  def manipulator(self):
    m = manipulator.ConfigurationManipulator()
    for parameter in self.parameter_space.forward_metrics[self.benchmark]:
      parameter = self.parameter_space.parameters[parameter]
      if "vectorization" not in parameter.tags or parameter.tags["vectorization"] == terminal.get_vectorization():
        if type(parameter) == generate.IntegerParameter:
          m.add_parameter(manipulator.IntegerParameter(parameter.name, parameter.minimum/parameter.step, parameter.maximum/parameter.step))
        if type(parameter) == generate.BooleanParameter:
          m.add_parameter(manipulator.BooleanParameter(parameter.name))
        if type(parameter) == generate.PowerOfTwoParameter:
          m.add_parameter(manipulator.PowerOfTwoParameter(parameter.name, parameter.minimum, parameter.maximum))
    return m

def create_benchmark_api(benchmark):
  args = argparser.parse_args()
  args.benchmark = benchmark
  args.database = "{}_{}{}".format(os.path.splitext(args.database)[0], benchmark, os.path.splitext(args.database)[1])
  interface = ReproBLASTuner(args)
  api = TuningRunManager(interface, args)
  return api

def run(apis):
  #setup
  args = argparser.parse_args()
  parameter_space = generate.deserialize_parameter_space(args.params)
  arguments_file_name = args.args
  if os.path.isfile(arguments_file_name):
    arguments = generate.deserialize_arguments(arguments_file_name)
  else:
    arguments = parameter_space.get_default_arguments()

  for i in range(args.trials):
    print("{}/{}".format(i, args.trials))

    #create a new argument file
    test_arguments = copy.deepcopy(arguments)
    apis_to_run = []
    desired_results = []
    for api in apis:
      desired_result = api.get_next_desired_result()
      if not desired_result:
        continue
      desired_results.append(desired_result)
      for (parameter, argument) in desired_result.configuration.data.items():
        parameter = parameter_space.parameters[parameter]
        if type(parameter) == generate.IntegerParameter:
          test_arguments[parameter.name] = argument * parameter.step
        if type(parameter) == generate.BooleanParameter:
          test_arguments[parameter.name] = argument
        if type(parameter) == generate.PowerOfTwoParameter:
          test_arguments[parameter.name] = argument
      apis_to_run.append(api)

    if apis_to_run:
      #build with these arguments
      generate.serialize_arguments(test_arguments, arguments_file_name)
      terminal.make_clean("./")

      bench_tests = []
      command_list = []
      for api in apis_to_run:
        bench_test = benchs.all_benchs[api.measurement_interface.benchmark]
        bench_tests.append(bench_test[0]())
        bench_tests[-1].setup(flagss = bench_test[1], attribute="%peak", args=arguments_file_name, remake=True, verbose=args.verbose)
        command_list += bench_tests[-1].get_command_list()

      #run with these arguments
      output_list = config.run(command_list, verbose=args.verbose)

      #return the results to the apis
      for api, desired_result, bench_test in zip(apis_to_run, desired_results, bench_tests):
        bench_test.parse_output_list(output_list[:len(bench_test.get_command_list())])
        output_list = output_list[len(bench_test.get_command_list()):]
        result = Result(time=(100.0/bench_test.get_result()))
        api.report_result(desired_result, result)

  #parse the best configurations
  best_arguments = copy.deepcopy(arguments)
  for api in apis:
    api.search_driver.process_new_results()
    for (parameter, argument) in api.get_best_configuration().items():
      parameter = parameter_space.parameters[parameter]
      if type(parameter) == generate.IntegerParameter:
        best_arguments[parameter.name] = argument * parameter.step
      if type(parameter) == generate.BooleanParameter:
        best_arguments[parameter.name] = argument
      if type(parameter) == generate.PowerOfTwoParameter:
        best_arguments[parameter.name] = argument
    api.finish()
  generate.serialize_arguments(best_arguments, arguments_file_name)

def main():
  print("tuning {}".format(terminal.get_vectorization()))
  """
  apis = [create_benchmark_api(benchmark) for benchmark in ["bench_damax",
                                                            "bench_damaxm",
                                                            "bench_zamax",
                                                            "bench_zamaxm",
                                                            "bench_samax",
                                                            "bench_samaxm",
                                                            "bench_camax",
                                                            "bench_camaxm"]]
  run(apis)
  fold = 0
  apis = [create_benchmark_api(benchmark) for benchmark in ["bench_rdsum_fold_{}".format(fold),
                                                            "bench_rdasum_fold_{}".format(fold),
                                                            "bench_rdnrm2_fold_{}".format(fold),
                                                            "bench_rddot_fold_{}".format(fold),
                                                            "bench_rzsum_fold_{}".format(fold),
                                                            "bench_rdzasum_fold_{}".format(fold),
                                                            "bench_rdznrm2_fold_{}".format(fold),
                                                            "bench_rzdotu_fold_{}".format(fold),
                                                            "bench_rzdotc_fold_{}".format(fold),
                                                            "bench_rssum_fold_{}".format(fold),
                                                            "bench_rsasum_fold_{}".format(fold),
                                                            "bench_rsnrm2_fold_{}".format(fold),
                                                            "bench_rsdot_fold_{}".format(fold),
                                                            "bench_rcsum_fold_{}".format(fold),
                                                            "bench_rscasum_fold_{}".format(fold),
                                                            "bench_rscnrm2_fold_{}".format(fold),
                                                            "bench_rcdotu_fold_{}".format(fold),
                                                            "bench_rcdotc_fold_{}".format(fold)]]
  run(apis)
  apis = []
  for fold in range(2, terminal.get_dimaxfold() + 1):
    apis += [create_benchmark_api(benchmark) for benchmark in ["bench_rdsum_fold_{}".format(fold),
                                                               "bench_rdasum_fold_{}".format(fold),
                                                               "bench_rdnrm2_fold_{}".format(fold),
                                                               "bench_rddot_fold_{}".format(fold),
                                                               "bench_rzsum_fold_{}".format(fold),
                                                               "bench_rdzasum_fold_{}".format(fold),
                                                               "bench_rdznrm2_fold_{}".format(fold),
                                                               "bench_rzdotu_fold_{}".format(fold),
                                                               "bench_rzdotc_fold_{}".format(fold)]]
  for fold in range(2, terminal.get_simaxfold() + 1):
    apis += [create_benchmark_api(benchmark) for benchmark in ["bench_rssum_fold_{}".format(fold),
                                                               "bench_rsasum_fold_{}".format(fold),
                                                               "bench_rsnrm2_fold_{}".format(fold),
                                                               "bench_rsdot_fold_{}".format(fold),
                                                               "bench_rcsum_fold_{}".format(fold),
                                                               "bench_rscasum_fold_{}".format(fold),
                                                               "bench_rscnrm2_fold_{}".format(fold),
                                                               "bench_rcdotu_fold_{}".format(fold),
                                                               "bench_rcdotc_fold_{}".format(fold)]]
  run(apis)
"""
  apis = [create_benchmark_api(benchmark) for benchmark in ["bench_rdgemv_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rdgemv_TransA_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rdgemm_AvgTransA_AvgTransB_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rzgemv_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rzgemv_TransA_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rzgemm_AvgTransA_AvgTransB_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rsgemv_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rsgemv_TransA_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rsgemm_AvgTransA_AvgTransB_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rcgemv_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rcgemv_TransA_fold_{}".format(terminal.get_didefaultfold()),
                                                            "bench_rcgemm_AvgTransA_AvgTransB_fold_{}".format(terminal.get_didefaultfold())]]
  run(apis)

if __name__ == '__main__':
  opentuner.init_logging()
  main()

