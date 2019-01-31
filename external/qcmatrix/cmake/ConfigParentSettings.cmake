# Settings from the host program
IF(NOT "${PARENT_DEFINITIONS}" STREQUAL "")
    FOREACH (_definition ${PARENT_DEFINITIONS})
        ADD_DEFINITIONS(${_definition})
    ENDFOREACH()
ENDIF()
IF(NOT "${PARENT_INCLUDE_DIR}" STREQUAL "")
    INCLUDE_DIRECTORIES(${PARENT_INCLUDE_DIR})
ENDIF()
IF(NOT "${PARENT_MODULE_DIR}" STREQUAL "")
    SET(CMAKE_Fortran_MODULE_DIRECTORY ${PARENT_MODULE_DIR})
ENDIF()
