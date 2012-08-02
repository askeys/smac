"""Python script builds the bindings between the SMAC code and the C++"""
import os, sys
from pyplusplus import module_builder

def createBindings(pathOfCode):
    """
    Create bindings around any code
    """
    #list_of_files = [pathOfCode +"/" + f for f in os.listdir(pathOfCode) if '.h' in f or '.hpp' in f or '.c' in f or '.cpp' in f]
    #print list_of_files
    print pathOfCode
    headerFile = "../../../src/smac.h"
    mb = module_builder.module_builder_t( [headerFile]
                                      , gccxml_path=r"/usr/local/bin/gccxml" 
                                      , working_directory=pathOfCode
                                      , include_paths=[pathOfCode, pathOfCode]
                                      , define_symbols=[] )


    #Well, don't you want to see what is going on?
    mb.print_declarations()
    
    #Creating code creator. After this step you should not modify/customize declarations.
    mb.build_code_creator( module_name='libsmac' )
    
    #Writing code to file.
    mb.write_module( './automatic_bindings.cpp' )

if __name__ == '__main__':
    pathOfCode = sys.argv[1]
    createBindings(pathOfCode)