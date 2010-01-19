# Class to implement the AccTag enum
# Define the enum values in the dict from a SWIG generated AccTag_enum class
# Don't allow these enum attributes to be set
# (But allow other attributes to be added)

class AccTag_enum(AccTag_enum):
   def __init__(self):
      self.this = "None"
      self.__dict__.update({'interior' : AccTag_enum.interior, 'exterior' : AccTag_enum.exterior, 'undecided' : AccTag_enum.undecided, 'in_tube' : AccTag_enum.in_tube})
   def __setattr__(self, name, value):
      if (name == "interior") or (name == "exterior") or (name == "undecided") or (name == "in_tube"):
         error = "Can't set constant enum attribute " + name
         raise AttributeError, error
      else:
         self.__dict__[name] = value

# Create an instance
AccTag = AccTag_enum()
