"""
The python interface to the libsmac files exists in this mini package. 
The purpose of this mini package is to offer a more robust *pythonic*
interface to the SMAC C++ distribution.

All of the original interface to the original SMAC C++ interface can
be accessed through the libsmac variable.
>>> from pysmac import libsmac
"""

import libsmac

# The pysmac interface
from pysmac.formats import BaseFormat, FormatError, HOOMDFormat
from pysmac.match import BaseCompare, ZernikeCompare, Matcher, CompareError