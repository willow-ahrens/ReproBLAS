import scripts.terminal as terminal

def run(command_list):
  """
  A function that runs the command list (a list of string commands) on the
  target and returns a list of their results as a list of tuples of (return
  code, output)
  """
  return[terminal.callsafe(command) for command in command_list]

version = "0.0.0"
