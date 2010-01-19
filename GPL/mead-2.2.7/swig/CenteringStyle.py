# Class to implement the CenteringStyle enum
# Define the enum values in the dict from a SWIG generated CenteringStyle_enum class.
# Don't allow these enum attributes to be set
# (But allow other attributes to be added)

class CenteringStyle_enum(CenteringStyle_enum):
   def __init__(self):
      self.this = "None"
      self.__dict__.update({'ON_ORIGIN' : CenteringStyle_enum.ON_ORIGIN, 'ON_CENT_OF_INTR' : CenteringStyle_enum.ON_CENT_OF_INTR, 'ON_GEOM_CENT' : CenteringStyle_enum.ON_GEOM_CENT, 'SPECIFIC' : CenteringStyle_enum.SPECIFIC})
   def __setattr__(self, name, value):
      if (name == "ON_ORIGIN") or (name == "ON_CENT_OF_INTR") or (name == "ON_GEOM_CENT") or (name == "SPECIFIC"):
         error = "Can't set constant enum attribute " + name
         raise AttributeError, error
      else:
         self.__dict__[name] = value

# Create an instance
CenteringStyle = CenteringStyle_enum()
