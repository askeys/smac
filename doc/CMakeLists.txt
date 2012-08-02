if (DOXYGEN)

# disable DOT in doxygen if dot is not found
if(DOXYGEN_DOT_PATH)
    set(DOXYGEN_HAVE_DOT "YES")
else(DOXYGEN_DOT_PATH)
    set(DOXYGEN_HAVE_DOT "NO")
endif(DOXYGEN_DOT_PATH)

# configure the doxygen file
configure_file (${SMAC_SOURCE_DIR}/doc/Doxyfile.in ${SMAC_BINARY_DIR}/doc/Doxyfile)

add_custom_command (OUTPUT ${SMAC_BINARY_DIR}/doc/html/index.html
                    COMMAND ${DOXYGEN} ${SMAC_BINARY_DIR}/doc/Doxyfile
                    #DEPENDS ${SMAC_SRC}
                    ${SMAC_SOURCE_DIR}/doc/*.dox
                    ${SMAC_SOURCE_DIR}/src/*.h
                    ${SMAC_SOURCE_DIR}/src/*.cpp
                    ${SMAC_SOURCE_DIR}/doc/Doxyfile.in)
                    #${SMAC_SOURCE_DIR}/doc/DoxygenLayout.xml)

add_custom_target (doc ALL DEPENDS ${SMAC_BINARY_DIR}/doc/html/index.html)

configure_file(${SMAC_SOURCE_DIR}/resources/matching.jpg ${SMAC_BINARY_DIR}/doc/html/matching.jpg COPYONLY)

endif (DOXYGEN)
