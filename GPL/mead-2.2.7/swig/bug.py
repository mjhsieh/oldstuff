"""Test the ability of the MEAD wrapper to catch a C++ exception
from MEAD (MEADexecpt) and transform it into a python exception."""

from PyMead import MEAD
import Numeric


a = MEAD.AtomSet()
print "Here we go making an error!"
a.read("junk")
print "Should never get here"

