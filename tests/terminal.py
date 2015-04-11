def callsafe(command):
  rc = 0
  try:
    out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True).decode(sys.stdout.encoding)
  except subprocess.CalledProcessError as e:
    rc = e.returncode
    out = e.output.decode(sys.stdout.encoding)
  return (rc, out)

def call(command):
  return subprocess.check_output(command, shell=True).decode(sys.stdout.encoding)

top = call("make top")

def make_clean(location):
  call("cd {0}; make clean".format(location).split()[-1], os.path.split(executable)[1])

def make(executable, args = None, id = None):
  executable = os.path.join(top,executable)
  result = os.path.join(call("cd {0}; make pbd".format(os.path.split(executable)[0])).split()[-1], os.path.split(executable)[1])
  env = ""
  if args:
    env = "ARGS={}".format(args)
  call("make {} {}".format(result, env))
  assert os.path.isfile(result), "Error: make unsuccessful."
  if id:
    new_result = "{}__{}{}".format(os.path.splitext(result)[1], id, os.path.splitext(result[0]))
    call("cp {} {}".format(result, new_result)
    result = new_result
  return result

def flags(params, args):
  return " ".join(['-{0} "{1}"'.format(param, arg) if len(param) == 1 else '--{0} "{1}"'.format(param, arg) for (param, arg) in zip(params, args)])


