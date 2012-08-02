
OPTION(BUILD_TESTING "build testing" 1)

#################################


IF(BUILD_TESTING)

INCLUDE_DIRECTORIES(${SMAC_SOURCE_DIR}/src)
ADD_EXECUTABLE(rmsdescriptor_test rmsdescriptor_test.cpp)
TARGET_LINK_LIBRARIES(rmsdescriptor_test smac)
ADD_EXECUTABLE(fourierdescriptor2d_test fourierdescriptor2d_test.cpp)
TARGET_LINK_LIBRARIES(fourierdescriptor2d_test smac)
ADD_EXECUTABLE(fourierdescriptor3d_test fourierdescriptor3d_test.cpp)
TARGET_LINK_LIBRARIES(fourierdescriptor3d_test smac)
ADD_EXECUTABLE(localorder_test localorder_test.cpp)
TARGET_LINK_LIBRARIES(localorder_test smac)
ADD_EXECUTABLE(registericp_test registericp_test.cpp)
TARGET_LINK_LIBRARIES(registericp_test smac)
ADD_EXECUTABLE(matchfun_test matchfun_test.cpp)
TARGET_LINK_LIBRARIES(matchfun_test smac)
ADD_EXECUTABLE(loadxyz_test loadxyz_test.cpp)
TARGET_LINK_LIBRARIES(loadxyz_test smac)
ADD_EXECUTABLE(loadlammps_test loadlammps_test.cpp)
TARGET_LINK_LIBRARIES(loadlammps_test smac)

#    FILE(GLOB UNIT_TESTS "${SMAC_SOURCE_DIR}/unit_tests/*.cpp")
#    STRING(REPLACE ".cpp" "" UNIT_TESTS "${UNIT_TESTS}")
#    FOREACH(_UNIT_TEST ${UNIT_TESTS})
#        INCLUDE_DIRECTORIES(${SMAC_SOURCE_DIR}/src)
#        ADD_EXECUTABLE(${_UNIT_TEST} ${_UNIT_TEST}.cpp)
#        TARGET_LINK_LIBRARIES(${_UNIT_TEST} smac)
#    ENDFOREACH(_UNIT_TEST ${UNIT_TESTS})

ENDIF(BUILD_TESTING)

#################################
