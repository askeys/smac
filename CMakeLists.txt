# The name of our project is "smac".  CMakeLists files in this project can
# refer to the root source directory of the project as ${SMAC_SOURCE_DIR} and
# to the root binary directory of the project as ${SMAC_BINARY_DIR}.

cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

PROJECT(SMAC)

# Compiler and path settings:
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZATION_FLAGS}" CACHE INTERNAL "Compiler optimization flags.")
SET(EXECUTABLE_OUTPUT_PATH ${SMAC_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for executables.")
SET(LIBRARY_OUTPUT_PATH ${SMAC_BINARY_DIR}/lib CACHE INTERNAL "Single output directory for libraries.")
SET(OPTIMIZATION_FLAGS "-O3 -funroll-loops -Wall" CACHE STRING "Compiler optimization flags.")
SET(SWIG_DIR /usr/local CACHE STRING "directory containing SWIG.")

#force a default build type if it's not set.
IF(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RELEASE")
ENDIF()

# Wrappers:

# Define statements:
OPTION(ENABLE_DOXYGEN "use Doxgen to build developer documentation" 1)
OPTION(ENABLE_ICP "use the ICP object registration package" 0)
OPTION(ENABLE_PCA "use the PCA object registration package" 0)
OPTION(ENABLE_HUNGARIAN "use Hungarian" 0)
OPTION(ENABLE_IPGZ "use IPGZ for filtering/image analysis" 0)
OPTION(ENABLE_GEOM "use geometry package for manipulating geometrical data" 0)
OPTION(ENABLE_VMD "use the VMD molecular visualizer" 0)
OPTION(WRAP_PYTHON "use structure library from python" 0)
#OPTION(WRAP_JAVA "use structure library from java" 0)

# Advanced options in the ccmake gui:

#Set up definitions and libraries:

#################################

#write a configure file 
FILE(WRITE ${SMAC_SOURCE_DIR}/src/config.h "//configure file for smac\n")
FILE(APPEND ${SMAC_SOURCE_DIR}/src/config.h "#define _SMAC_SRC_PATH \"${SMAC_SOURCE_DIR}\"\n")

#################################

IF(ENABLE_VMD)
	SET(VMD_EXECUTABLE_PATH "/Applications/VMD\\ 1.8.6.app/Contents/vmd/vmd*" CACHE STRING "path to the program VMD")
	FILE(APPEND ${SMAC_SOURCE_DIR}/src/config.h "#define _USE_VMD 1\n")
	FILE(APPEND ${SMAC_SOURCE_DIR}/src/config.h "#define _VMD_EXECUTABLE ${VMD_EXECUTABLE_PATH}\n")	
ENDIF(ENABLE_VMD)

IF(ENABLE_IPGZ)
	FILE(APPEND ${SMAC_SOURCE_DIR}/src/config.h "#define _USE_IPGZ 1\n")	
	#this is all very half-assed:
	AUX_SOURCE_DIRECTORY(${SMAC_SOURCE_DIR}/3rdparty/ipgz IPGZ_SRC)
	SET(LIBRARIES_USED ${LIBRARIES_USED};tiff)
	INCLUDE_DIRECTORIES(/sw/include)
ENDIF(ENABLE_IPGZ)

IF(ENABLE_GEOM)
	AUX_SOURCE_DIRECTORY(${SMAC_SOURCE_DIR}/3rdparty/geometry GEOM_SRC)
	FILE(APPEND ${SMAC_SOURCE_DIR}/src/config.h "#define _USE_GEOM 1\n")	
ENDIF(ENABLE_GEOM)

IF(ENABLE_HUNGARIAN)
MESSAGE(WARNING "Enabling the Hungarian library causes the program to fall under the jurisdiction of the GNU public license.  Disable Hungarian to revert to original license.")
AUX_SOURCE_DIRECTORY(${SMAC_SOURCE_DIR}/3rdparty/hungarian HUNGARIAN_SRC)
	FILE(APPEND ${SMAC_SOURCE_DIR}/src/config.h "#define _USE_HUNGARIAN 1\n")	
ENDIF(ENABLE_HUNGARIAN)

IF(ENABLE_ICP)
	AUX_SOURCE_DIRECTORY(${SMAC_SOURCE_DIR}/3rdparty/icp ICP_SRC)
	FILE(APPEND ${SMAC_SOURCE_DIR}/src/config.h "#define _USE_ICP 1\n")	
ENDIF(ENABLE_ICP)

IF(ENABLE_PCA)
	AUX_SOURCE_DIRECTORY(${SMAC_SOURCE_DIR}/3rdparty/pca PCA_SRC)
	FILE(APPEND ${SMAC_SOURCE_DIR}/src/config.h "#define _USE_PCA 1\n")	
ENDIF(ENABLE_PCA)

INCLUDE_DIRECTORIES(${SMAC_SOURCE_DIR}/3rdparty)

SET(THIRD_PARTY_SRC ${HUNGARIAN_SRC};${ICP_SRC};${PCA_SRC};${IPGZ_SRC};${GEOM_SRC} CACHE INTERNAL "Third party libraries")

#Set up source files:
AUX_SOURCE_DIRECTORY(src SMAC_SRC)

#################################
find_package(Doxygen)
IF(DOXYGEN_FOUND)
	OPTION(ENABLE_DOXYGEN "Enables building of documentation with doxygen" 1)
	IF(ENABLE_DOXYGEN)
		ADD_SUBDIRECTORY(doc)
		#ADD_CUSTOM_TARGET(doc ALL ${DOXYGEN_EXECUTABLE} ${TPS_SOURCE_DIR}/Doxyfile)
	ENDIF(ENABLE_DOXYGEN)
ENDIF(DOXYGEN_FOUND)
#################################

#################################

#make smac.h header file:

FILE(GLOB HEADER_FILES "${SMAC_SOURCE_DIR}/src/*.h")
STRING(REPLACE "${SMAC_SOURCE_DIR}/src/smac.h;" "" HEADER_FILES "${HEADER_FILES}")
STRING(REPLACE "${SMAC_SOURCE_DIR}/src/" "" HEADER_FILES "${HEADER_FILES}")
FILE(WRITE ${SMAC_SOURCE_DIR}/src/smac.h "#ifndef _SMAC_H\n#define _SMAC_H\n\n")
FOREACH(_FILENAME ${HEADER_FILES})
	FILE(APPEND ${SMAC_SOURCE_DIR}/src/smac.h "#include \"${_FILENAME}\"\n")
ENDFOREACH(_FILENAME ${HEADER_FILES})
FILE(APPEND ${SMAC_SOURCE_DIR}/src/smac.h "\n#endif\n\n") 

#install headers:

FILE(GLOB files "${SMAC_SOURCE_DIR}/src/*.h*")
INSTALL(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/smac)

################################

# Create a library called "smac" 
# note: some of these variables will be empty if not defined

ADD_LIBRARY(smac 
	${SMAC_SRC} ${THIRD_PARTY_SRC}
)

TARGET_LINK_LIBRARIES(smac ${LIBRARIES_USED})

INSTALL(TARGETS smac
			RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
			LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
			
IF(WRAP_PYTHON)
    FIND_PACKAGE(Boost)
    FIND_PACKAGE(PythonLibs)
    IF(Boost_FOUND)
        FIND_PACKAGE(Boost COMPONENTS python REQUIRED)
        INCLUDE_DIRECTORIES("wrap/python/py++/" ${PYTHON_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
        SET(Boost_USE_STATIC_LIBS OFF)
        SET(Boost_USE_STATIC_RUNTIME OFF)
        execute_process ( COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE )
        ADD_LIBRARY(smacpy SHARED wrap/python/py++/bindings.cpp)
        TARGET_LINK_LIBRARIES(smacpy ${Boost_PYTHON_LIBRARY} ${PYTHON_LIBRARY} ${LIBRARIES_USED} smac)
        SET_TARGET_PROPERTIES(smacpy PROPERTIES PREFIX "lib" SUFFIX ".so")
        INSTALL(PROGRAMS lib/libsmac.so DESTINATION ${PYTHON_SITE_PACKAGES})
        
    ENDIF(Boost_FOUND)
ENDIF(WRAP_PYTHON)


#################################
ADD_SUBDIRECTORY(unit_tests)
#################################
