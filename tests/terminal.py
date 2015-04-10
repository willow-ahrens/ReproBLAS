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

built = {}

def make(executable):
  global built
  executable = os.path.join(top,executable)
  if executable not in built:
    built[executable] = os.path.join(call("cd {0}; make pbd".format(os.path.split(executable)[0])).split()[-1], os.path.split(executable)[1])
    call("make {0}".format(built[executable]))
    assert os.path.isfile(built[test]), "Error: make unsuccessful."
  return built[executable]

def flags(params, args):
  return " ".join(['-{0} "{1}"'.format(param, arg) if len(param) == 1 else '--{0} "{1}"'.format(param, arg) for (param, arg) in zip(params, args)])


