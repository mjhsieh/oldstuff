import exceptions

class Error(exceptions.Exception):
   def __init__(self, args=None):
      self.args = args
