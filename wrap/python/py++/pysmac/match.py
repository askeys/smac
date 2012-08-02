"""
Match file wraps the comparing of two objects
into a well oiled interface.
"""
import libsmac
import os, sys
from abc import ABCMeta, abstractmethod, abstractproperty
from pysmac.formats import BaseFormat, getFormats, FormatError

class CompareError(Exception):
    """Exception raised when there has been a comparison error"""
    def __init__(self, compare, qualifier):
        """
        compare - A class{BaseCompare} object
        qualifier - A string qualifier highlighting why I failed
        """
        self.compare = compare
        self.qualifier = qualifier

    def __str__(self):
        """What I tell the user when I am raised"""
        return "CompareError:" + repr(self.compare) + self.qualifier

class BaseCompare():
    """Abstract class for all other compare classes"""
    __metaclass__ = ABCMeta

    def __init__(self,format,format2):
        """
        Compares two BaseFormat Objects."
        format - @class{BaseFormat} the point information
        format2 - @class{BaseFormat} the point information
        """
        if not isinstance(format, BaseFormat):
            raise CompareError(self, "FourierCompare.__init__ Format object is not a BaseFormat object %s" % (str(format)))
        if not isinstance(format2, BaseFormat):
            raise CompareError(self, "FourierCompare.__init__ Format object is not a BaseFormat object %s" % (str(format2)))

        # Save a reference to the base compare objects
        self.format = format
        self.format2 = format2

        # Class variables
        self.dimension = 0
        self.match = {}
        self.desc = {}
        self.desc2 = {}
        self.descInfo = None
       
        self.dirtyinfo = True #variables to force update of info if descInfo changes or match type changes
        self.dirtymatch = False
        self.MultiMatch = False

        # Check the dimension of the data
        # todo add pad 3d option to this
        if self.format.getDimension() != self.format2.getDimension():
            raise CompareError(self, " FourierCompare.__init__ Format data is not in the same dimension %s:%s" %
                              (str(self.format), str(self.format2)))
        else:
            self.dimension = self.format.getDimension()
            
        self._checkflags()
    
    @abstractproperty
    def _match(self):
        """internal function for heavy lifting of performing a match"""
        
    @abstractproperty
    def _multimatch(self):
        """internal function for heavy lifting of performing a multimatch"""
    
    @abstractproperty
    def _bestmatch(self,exclude):
        """internal function for heavy lifting of performing a bestmatch"""

    @abstractproperty
    def _checkflags(self):
        """internal function to check flags and perform internal data updates"""


    def getTypeMatch(self,type1,type2=None):
        """returns the specific value of of the match dictionary corresponding to type1 and type2
        type2 = none, returns a value based only on type1"""
        self._checkflags()
        if type2 is None or not self.MultiMatch:
            return self.match[type1]
        else:
            return self.Mmatch[type1][type2]

    def getMeanMatch(self):
        """Get the mean match between two properties over all types"""
        self._checkflags()
        if self.MulitMatch:
            return None #todo
        else:
            mData = [f for f in self.match.itervalues()]
            return sum(mData) / len(mData)
            
    def getMatch(self):
        """Return a dictionary hashing match information by type"""
        self._checkflags()
        if self.MultiMatch:
            return self.Mmatch
        else:
            return self.match
    
    def getBestMatch(self,exclude=False):
        """returns best match for the type"""
        self._checkflags()
        self._bestmatch(exclude)
        return self.BestMatch


    def setMultiMatch(self,Multi=True):
        """Allows switching of match type"""
        self.MultiMatch = bool(Multi)
        self.dirtymatch = True

    def getMatchType(self):
        """returns current match type"""
        return MatchType

    @abstractmethod
    def setDescInfo(self,info):
        """Allows setting of info paramaters"""

    @abstractmethod
    def getDescInfo(self):
        """returns dictionary of *_info data"""

class ZernikeCompare(BaseCompare):
    """compares two formats with a well defined interface using ZernikeDescriptors"""
    def _checkflags(self):
        if self.dirtyinfo:
            self._infomaker()
            self.dirtyinfo = True
            self.dirtymatch = True
        if self.dirtymatch:
            if self.MultiMatch:
                self._multimatch()
                
            self._match()
            self.dirtymatch = False

    def _getdescriptors(self):
        #internal function to get the descriptors of all the data.
        pI = self.format.getPointDictionary()
        pI2 = self.format2.getPointDictionary()

        for type1 in pI:
            shapeData = libsmac.shapedata_t(pI[type1])
            desc = libsmac.shpdesc_t()

            if self.dimension == 2:
               libsmac.zerndesc2d(shapeData,desc,self.descInfo)
            elif self.dimension == 3:
               libsmac.zerndesc3d(shapeData,desc,self.descInfo)

            self.desc[type1] = desc

        for type1 in pI2:
            shapeData = libsmac.shapedata_t(pI2[type1])
            desc = libsmac.shpdesc_t()

            if self.dimension == 2:
               libsmac.zerndesc2d(shapeData,desc,self.descInfo)
            elif self.dimension == 3:
               libsmac.zerndesc3d(shapeData,desc,self.descInfo)

            self.desc2[type1] = desc


    def _infomaker(self):
        #internal function to generate a fourier_info object
        if self.descInfo == None:
                self.descInfo = libsmac.zernike_info()
                self.infodict = {}
        self._getdescriptors()
        self.dirtymatch = True



    def _multimatch(self):
        """internal function for calculating the multimatch"""
        #so default outputs still work
        self.Mmatch = dict(dict())
        for type1, descData in self.desc.iteritems():
            tmatchdict = {}
            for type2, descData2 in self.desc2.iteritems():
                tmatchdict[type2] = libsmac.matchfundist(self.desc[type1],self.desc2[type2])
            self.Mmatch[type1] = tmatchdict
            del tmatchdict
        self.dirtymatch = False



    def _bestmatch(self,exclude):
        """Internal function for computing the best match to each type in format1 """
        if self.MultiMatch:
            self.BestMatch = {}
            for type1 in self.Mmatch:
                bestval = -1000
                besttype = None
                for type2 in self.Mmatch[type1]:
                    if (type1 is not type2) or not exclude:
                        if self.Mmatch[type1][type2] > bestval:
                            besttype = type2
                            bestval = self.Mmatch[type1][type2]
                self.BestMatch[type1] = {besttype:self.Mmatch[type1][besttype]}
                    

    def _match(self):
        """Internal function used to perform the dirty work of matching the objects type-to-type"""
        try:
            for type1 in self.desc:
                self.match[type1] = libsmac.matchfundist(self.desc[type1], self.desc2[type1])
            self.dirtymatch = False
        except KeyError, e:
            raise CompareError(self, "FourierCompare._match Format data does not contain the same type information %s: %s. \n perhaps you meant to use multimatch?\n" %
                               (str(self.format), str(self.format2)))
    def getDescInfo(self,key=None):
        """gets the current setting 'key' from the descInfo object.  If key is None return the whole dictionary"""
        if key is not None:
            return self.infodict[key]
        else:
            return self.infodict
        
        
    def setDescInfo(self,valDict):
        """set parameters of the descInfo object
        moment - list - Zernike moments to compute
        invariant - bool -Compute invariant moments
        normethod - int - how to normalize the shape
        rmax - int - maximum radius for normalized shape"""
        for type1 in valDict:
            if type1 == "moment":
                self.infodict[type1] = valDict[type1]
                self.descInfo.moment[:] = valDict[type1]
            elif type1 == "bool":
                self.infodict[type1] = valDict[type1]
                self.descInfo.invariant = valDict[type1]
            elif type1 == "invariant":
                self.infodict[type1] = valDict[type1]
                self.descInfo.invariant = valDict[type1]
            elif type1 == "normmethod":
                self.infodict[type1] = valDict[type1]
                self.descInfo.normmethod = valDict[type1]
            elif type1 == "rmax":
                self.infodict[type1] = valDict[type1]
                self.descInfo.rmax = valDict[type1]
            else:
                raise ValueWarning("dictionary contains extraneous information")
            self.dirtyinfo = True
    def getDescriptorDict(self):
        """returns a dictionary of Descriptor data for argument1"""
        return self.infodict

class FourierCompare(BaseCompare):
    """compares two formats with a well defined interface using Fourier Descriptors"""
    def _checkflags(self):
        if self.dirtyinfo:
            self._infomaker()
            self.dirtyinfo = True
            self.dirtymatch = True
        if self.dirtymatch:
            if self.MultiMatch:
                self._multimatch()
                
            self._match()
            self.dirtymatch = False

    def _getdescriptors(self):
        #internal function to get the descriptors of all the data.
        pI = self.format.getPointDictionary()
        pI2 = self.format2.getPointDictionary()

        for type1 in pI:
            shapeData = libsmac.shapedata_t(pI[type1])
            desc = libsmac.shpdesc_t()

            if self.dimension == 2:
               libsmac.fourierdesc2d(shapeData,desc,self.descInfo)
            elif self.dimension == 3:
               libsmac.fourierdesc3d(shapeData,desc,self.descInfo)

            self.desc[type1] = desc

        for type1 in pI2:
            shapeData = libsmac.shapedata_t(pI2[type1])
            desc = libsmac.shpdesc_t()

            if self.dimension == 2:
               libsmac.fourierdesc2d(shapeData,desc,self.descInfo)
            elif self.dimension == 3:
               libsmac.fourierdesc3d(shapeData,desc,self.descInfo)

            self.desc2[type1] = desc


    def _infomaker(self):
        #internal function to generate a fourier_info object
        if self.descInfo == None:
                self.descInfo = libsmac.fourier_info()
                self.infodict = {"frequency":range(0,13), "invariant":1}
                self.descInfo.frequency[:] = range(0,13)
                self.descInfo.invariant = 1
                self.descInfo.nsectors = 0
                self.descInfo.trigtablesize = 0
        self._getdescriptors()
        self.dirtymatch = True



    def _multimatch(self):
        """internal function for calculating the multimatch"""
        #so default outputs still work
        self.Mmatch = dict(dict())
        for type1, descData in self.desc.iteritems():
            tmatchdict = {}
            for type2, descData2 in self.desc2.iteritems():
                tmatchdict[type2] = libsmac.matchfundist(self.desc[type1],self.desc2[type2])
            self.Mmatch[type1] = tmatchdict
            del tmatchdict
        self.dirtymatch = False



    def _bestmatch(self,exclude):
        """Internal function for computing the best match to each type in format1 """
        if self.MultiMatch:
            self.BestMatch = {}
            for type1 in self.Mmatch:
                bestval = -1000
                besttype = None
                for type2 in self.Mmatch[type1]:
                    if (type1 is not type2) or not exclude:
                        if self.Mmatch[type1][type2] > bestval:
                            besttype = type2
                            bestval = self.Mmatch[type1][type2]
                self.BestMatch[type1] = {besttype:self.Mmatch[type1][besttype]}
                    

    def _match(self):
        """Internal function used to perform the dirty work of matching the objects type-to-type"""
        try:
            for type1 in self.desc:
                self.match[type1] = libsmac.matchfundist(self.desc[type1], self.desc2[type1])
            self.dirtymatch = False
        except KeyError, e:
            raise CompareError(self, "FourierCompare._match Format data does not contain the same type information %s: %s. \n perhaps you meant to use multimatch?\n" %
                               (str(self.format), str(self.format2)))
    def getDescInfo(self,key=None):
        """gets the current setting 'key' from the descInfo object.  If key is None return the whole dictionary"""
        if key is not None:
            return self.infodict[key]
        else:
            return self.infodict
        
        
    def setDescInfo(self,valDict):
        """set parameters of the descInfo object
        frequency[0..13] - list - frequencies to take fourier descriptor off
        shells[None] - list - define shells for psudo radial dependence
        invariant[2] - int{0,1,2} - 0:non-invariant 1:rotational invariance 2:wigner-seitz invariance
        trigtablesize[None] - int - #TODO
        nsectors[None] - int -#TODO"""
        for type1 in valDict:
            if type1 == "frequency":
                self.infodict[type1] = valDict[type1]
                self.descInfo.frequency[:] = valDict[type1]
            elif type1 == "shells":
                self.infodict[type1] = valDict[type1]
                self.descInfo.shells[:] = valDict[type1]
            elif type1 == "invariant":
                self.infodict[type1] = valDict[type1]
                self.descInfo.invariant = valDict[type1]
            elif type1 == "trigtablesize":
                self.infodict[type1] = valDict[type1]
                self.descInfo.trigtablesize = valDict[type1]
            elif type1 == "nsectors":
                self.infodict[type1] = valDict[type1]
                self.descInfo.nsectors = valDict[type1]
            else:
                raise ValueWarning("dictionary contains extraneous information")
            self.dirtyinfo = True
    def getDescriptorDict(self):
        """returns a dictionary of Descriptor data for argument1"""
        return self.infodict



class Matcher(object):
    """
    Class matches two different structures using a well-defined interface
    """
    def __init__(self, f, f2, compare):
        """
        Construct/Validate matcher object
        f - A file name containing point data
        f2 - Another file name containing point data
        compare - @class{BaseCompare} which implements compares the two files
        """
        # The format objects associated with the matcher
        self.format = None
        self.format2 = None

        # Extract the format information
        formatList = getFormats()
        for fM in formatList:
            if fM.isFormat(f):
                self.format = fM(f)
            if fM.isFormat(f2):
                self.format2 = fM(f2)
        if self.format is None:
            raise FormatError(self, "Matcher.__init__ unknown format in file %s" % (str(f)))
        if self.format2 is None:
            raise FormatError(self, "Matcher.__init__ unknown format in file %s" % (str(f2)))

        # Check the compare
        if not callable(compare):
            raise CompareError(self, "Matcher.__init__ object type is not callable" % (str(compare)))
        if not issubclass(compare, BaseCompare):
            raise CompareError(self, "Matcher.__init__ unknown compare type %s" % (str(compare)))

        # Create the compare object
        self.compare = compare(self.format, self.format2)

        # Get the match
        self.match = self.compare.getMeanMatch()

    def __cmp__(self, other):
        """Compare two matcher objects. Allows for sorting of matcher objects"""
        if not isisntance(other, Matcher):
            return False
        return cmp(self.match, other.match)



