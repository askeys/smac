"""
Formats file wraps certain file types for ease of use.
The BaseFormat class provides the minimal interface needed
to function with the BaseMatcher class.
Supported formats
HOOMDXML.
"""
# SMAC
import libsmac

# Instantiate an interface
from abc import ABCMeta, abstractmethod, abstractproperty
from os.path import split as pathsplit
from xml.etree.ElementTree import ElementTree

class FormatError(Exception):
    """Exception is raised when there has been a formatting error"""
    def __init__(self, format, qualifier):
        """
        @args - format - A @class{BaseFormat} object
        @args - qualifier - A string qualifier - why did I fail?
        """
        self.format = format
        self.qualifier = qualifier

    def __str__(self):
        """Return the string representation of the format"""
        return "FormatError: " + repr(self.format) + self.qualifier

def getFormats():
    """
    Returns a list of supported format. If you are writing a new format, make
    sure to add it to this list
    """
    return [HOOMDFormat, PythonList]

def getFormatStrings():
    """
    Returns a python list of format.getFormatName() from getFormats()
    """
    fList = getFormats()
    fStringList = []
    for i in range(0,len(fList)):
        fStringList.append(fList[i]().getFormatName())
    return fStringList

class BaseFormat():
    """The abstract class for all other format classes"""
    __metaclass__ = ABCMeta

    # These abstract methods, properties *MUST* be implemented
    # in the classes derived from base format
    @abstractproperty
    def getFormatName(self):
        """Get the format name (HOOMD XML, XYZ, etc.)"""
        pass

    @abstractproperty
    def getText(self):
        """Get the text from a filename, or file-like variable"""

    @abstractmethod
    def getPointDictionary(self):
        """
        Function takes the extracted information from the files and places the
        information into a dictionary{dict}. The dictionary is a hash between
        a particle type, and a list of coordinate data.
        """

    @abstractmethod
    def getDimension(self):
        """Returns the dimension of the data"""

    @abstractmethod
    def isFormat(self, f):
        """Checks if a file conforms to the format"""

class HOOMDFormat(BaseFormat):
    """
    Parses the HOOMD xml format, and extracts relevant point/type
    information for shape matching
    """
    def __init__(self, f=None):
        """
        Extract the relevant information from a hoomd file
        @args - f - A path string to a HOOMD file
        """
        # Extract the information
        self.filename = f
        self.text = ""
        try:
            self.text = open(f, 'r').read()
        except Exception, e:
            if self.filename is not None:
                raise FormatError(self, "HOOMDFormat.__init__ file could not be opened %s" % (str(e)))

        # Class variables
        self.formatName = "HOOMD"
        self.pointInformation = {}
        self.dimension = 0

        # Parse the information
        if self.filename is not None:
            self.parse()

    def parse(self):
        """
        Parse the text from a hoomd file. The relevant point information will
        be stored in pointInformation.
        TODO: Box, diameter parameters in the XML
        """
        # Create the tree and parse
        tree = ElementTree()

        # Parse
        try:
            tree.parse(self.filename)
        except Exception, e:
            raise FormatError(self, "HOOMDFormat.parse file is not well-formed XML" %(str(e)))

        # Extract the point information
        elem = tree.find("hoomd_xml/position")
        pointText = elem.text
        pointList = []
        for line in pointText.splitlines():
            pointList.append([float(a) for a in line.split()])

        # Extract the dimension of the information
        if len(pointList):
            self.dimension = len(pointList[0])

        # Extract the type information
        elem = tree.find("hoomd_xml/type")
        typeText = elem.text
        typeList = []
        for line in typeText.splitlines():
            typeList.append(a)

        # After filtering copy over to the coordlist
        data = {}
        for index, type in enumerate(typeText):
            if not data.has_key(type):
                data[key] = []
            data.key.append(pointList[index])

        # Copy over the data into coords, and save into the point information
        for type, pointList in data.iteritems():
            # Copy over the coordinates to C++
            coords = libsmac.coordlist_t()
            coords[:] = pointList

            # Save the object in the pointInformation
            self.pointInformation[type] = coords

    def __str__(self):
        """String representation of the hoomd file format"""
        return self.formatName + "{" + self.filename + ":" + self.dimension + "}"
    # Implementation of the BaseFormat interface
    def getText(self):
        return self.text

    def getFormatName(self):
        return self.formatName

    def getPointDictionary(self):
        return self.pointInformation

    def getDimension(self):
        return self.dimension

    def isFormat(self, f):
        """Check if the necessary nodes to match can be extracted from the xml"""
        try:
            tree = ElementTree()
            tree.parse(f)
            elem = tree.find("hoomd_xml/position")
            elem = tree.find("hoomd_xml/type")
        except Exception, e:
            return False
        return True

class PythonList(BaseFormat):
    """
    Converts a list of points to data suitable for matching
    """
    def __init__(self,pointList=None):
        #class variables
        self.formatName = "PyList"
        self.pointInformation = {}


        # Parse the information
        #self.pointList = pointList
        self.points = pointList
        if self.points is not None:
            self.parse()


    def parse(self):
        """Parse the list into proper data format while doing error checking"""
        from numpy import array
        import string


        #check that it it is 2 or 3 dimensional
        if array(self.points).ndim is 2:
            #we have a single list of points
            tempCoords = libsmac.coord_t();
            tempCoordlist = libsmac.coordlist_t()
            for i in range(0,len(self.points)):
                tempCoords[:] = self.points[i]
                tempCoordlist.append(tempCoords)
            #copy coordlist to our object
            self.coords = tempCoordlist
            self.pointInformation = {"a":self.coords}
            self.dimension = len(self.points[0])

        elif array(self.points).ndim is 3:
            #we have a multi list!
            self.dimension = array(self.points)[0][0].size
            if len(self.points) > 26:
                raise FormatError(self, "Python multilist has a maximum of 26 types allowed")

            for i in range(0,len(self.points)):
                tempCoords = libsmac.coord_t();
                tempCoordlist = libsmac.coordlist_t()
                for j in range(0,len(self.points[0])):
                    tempCoords[:] = self.points[i][j]
                    tempCoordlist.append(tempCoords)
                self.pointInformation[string.ascii_lowercase[i]] = tempCoordlist
                del tempCoords
                del tempCoordlist
        else:
            raise FormatError(self, "PythonList.parse list must be 2 or 3 dimensional")






    def getText(self):
        ##TODO
        return ""

    def getFormatName(self):
        return self.formatName

    def getPointDictionary(self):
        return self.pointInformation

    def getDimension(self):
        return self.dimension

    def isFormat(self):
        ##TODO
        return True

class UnsupportedFormat(BaseFormat):
    """ Format for reporting format errors.  This format does not actually work for matching."""
    def getFormatName(self):
       return "UnsupportedFormat"


    def getText(self):
        """Get the text from a filename, or file-like variable"""
        return ""

    def getPointDictionary(self):
        """
        Function takes the extracted information from the files and places the
        information into a dictionary{dict}. The dictionary is a hash between
        a particle type, and a list of coordinate data.
        """
        return {'a':[]}

    def getDimension(self):
        """Returns the dimension of the data"""
        return 0

    def isFormat(self, f):
        """Checks if a file conforms to the format"""
        return False


