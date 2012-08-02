#!/usr/bin/env python
#
#  matchlib.py
#  SMAC
#
#  Created by Ross Smith on 23/01/2011.
#  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
#
## File for automagic-ish interfaces

import pysmac
import pysmac.match as matchs
import pysmac.formats as formats

def ZernikeMatch(obj1,obj2,format='PyList'):
    """Abstraction function of the underlying matching algorithm and format objects"""
    " obj1, obj2 are two items formatted according to format"
    " format = all supported formats from formats.getFormatStrings()"
    " output = ZernikeMatch Object"
    ##todo implement a function to auto determine format of obj1 and obj2

    #get list of supported formats and add auto to them.
    SFormats = formats.getFormatStrings()
    SFormats.append('auto')

    #make sure we got a good format string
    try:
        formatIndex = SFormats.index(format)
    except ValueError:
        raise formats.FormatError(formats.UnsupportedFormat, "format %s is undefined" % format)

    # make format object out of the input
    if (format != 'auto'):
        fClassList = formats.getFormats()
        FObj1 = fClassList[formatIndex](obj1)
        FObj2 = fClassList[formatIndex](obj2)
    else :
        raise ValueError("not implemented")

    zernikeMatched = pysmac.match.ZernikeCompare(FObj1,FObj2);
    return zernikeMatched


def FourierMatch(obj1,obj2,format='PyList'):
    """abstraction function of the underlying matching algorithm and format objects"""
    " obj1, obj2, are two items formatted according to format"
    " format = all supported formats from formats.getFormatStrings()"
    " output = FourierMatch object"
    ##todo implement a function to auto determine format of obj1 and obj2

    #get list of supported formats and add auto to them.
    SFormats = formats.getFormatStrings()
    SFormats.append('auto')

    #make sure we got a good format string
    try:
        formatIndex = SFormats.index(format)
    except ValueError:
        raise formats.FormatError(formats.UnsupportedFormat, "format %s is undefined" % format)

    # make format object out of the input
    if (format != 'auto'):
        fClassList = formats.getFormats()
        FObj1 = fClassList[formatIndex](obj1)
        FObj2 = fClassList[formatIndex](obj2)
    else :
        raise ValueError("not implemented")

    FourierMatched = pysmac.match.FourierCompare(FObj1,FObj2);
    return FourierMatched
