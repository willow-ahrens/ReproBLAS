from src.gen import generate
from tests.benchs import benchs
from opentuner.resultsdb.models import Result, TuningRun
from opentuner.search import manipulator

class ScaledIntegerParameter(manipulator.ScaledNumericParameter, manipulator.IntegerParameter):
  """
  A scaled integer
  """

  def __init__(self, name, min_value, max_value, step, **kwargs):
    kwargs['value_type'] = int
    assert min_value % step == 0  # must be multiple of step
    assert max_value % step == 0  # must be multiple of step
    self.step = step
    super(PowerOfTwoParameter, self).__init__(name, min_value/step, max_value/step,
                                              **kwargs)

  def _scale(self, v):
    return self.step * v

  def _unscale(self, v):
    return v/self.step

class ReproBLASTuner(opentuner.measurement.MeasurementInterface):
  def __init__(self, benchmark, parameters_file_name, arguments_file_name, *pargs, **kwargs):
    super(ReproBLAS, self).__init__(project_name="ReproBLAS", program_name=benchmark, program_version=config.version, objective=MinimizeTime, *pargs, **kwargs)
    self.benchmark = benchmark
    self.benchmark_class = benchs.all_benchs[benchmark]
    self.parameters = generate.deserialize_parameters(parameters_file_name)
    self.arguments_file_name = arguments_file_name
    self.arguments = generate.deserialize_arguments(arguments_file_name)
    self.parallel_compile = False

  def manipulator(self):
    m = manipulator.ConfigurationManipulator()
    for parameter in self.parameters.backward_metrics[self.benchmark]:
      parameter = self.parameters.parameters[parameter]
      if type(parameter) == generate.IntegerParameter:
        m.add_parameter(ScaledIntegerParameter(parameter.name, parameter.minimum, parameter.maximum - parameter.step, parameter.step))
      if type(parameter) == generate.BooleanParameter:
        m.add_parameter(manipulator.BooleanParameter(parameter.name))

  def run(self, desired_result, input, limit):
    pass

  def compile(self, config_data, result_id):
    terminal.make_clean("src/")
    test_arguments = copy.deepcopy(self.arguments)
    for parameter in self.parameters.backward_metrics[self.benchmark]:
      test_arguments[parameter] = config_data[parameter]
    serialize_arguments(test_arguments, self.arguments_file_name):

    build_benchmark = self.benchmark_class()

    build_benchmark.setup(args = self.arguments_file_name, id=result_id)
    return build_benchmark

  def run_precompiled(self, desired_result, input, limit, compile_result, id):
    compile_result.parse_output_lines(config.run(compile_result.get_output_lines()))
    return Result(time=compile_result.get_result())

  def save_final_config(self, config):
    """
    called at the end of autotuning with the best resultsdb.models.Configuration
    """
    test_arguments = copy.deepcopy(self.arguments)
    for parameter in self.parameters.backward_metrics[self.benchmark]:
      test_arguments[parameter] = config_data[parameter]
    serialize_arguments(test_arguments, self.arguments_file_name):

  @classmethod
  def main(cls, args, *pargs, **kwargs):
    from opentuner.tuningrunmain import TuningRunMain

    return TuningRunMain(cls(args, *pargs, **kwargs), args).main()
