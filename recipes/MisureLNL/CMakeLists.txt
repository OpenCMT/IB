################################################################################
################################################################################
################################################################################

cmake_minimum_required(VERSION 2.6)

project(misureLNL)
set(SOURCES
    # put your source here
    main.cpp
)

################################################################################
################################################################################

################################################################################
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")

### MUDULE FIND PCL PACKAGE ------------------------------------------------------ ##
find_package(PCL 1.3 REQUIRED COMPONENTS common io)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

################################################################################
include(uLibDebugMacro)
include(uLibCommon)
include(uLibFindDependencies)
include(uLibConfigHeader)
include(uLibGenerateRMake)

################################################################################
## MUDULE FIND PACKAGE ------------------------------------------------------ ##
#find_package(uLib REQUIRED) # not working at the moment

## CONFIG FIND PACKAGE ------------------------------------------------------ ##
#set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}/lib/cmake")
unset(ULIB_CONFIG)
unset(ULIB_INCLUDE_DIRS)
find_package(ULIB 0.2 CONFIG
             NAMES uLib
             PATH_SUFFIXES mutom-0.2
)
if(ULIB_CONFIG)
 set(ULIB_FOUND true)
endif()

unset(IB_CONFIG)
find_package(IB 0.2 CONFIG
             PATH_SUFFIXES mutom-0.2
)
if(IB_CONFIG)
 set(IB_FOUND true)
endif()

## PROJECT EXECUTABLE ------------------------------------------------------- ##
add_executable(${PROJECT_NAME} ${SOURCES})
add_executable(IB_abtrim ../../utils/filters/IB_abtrim.cpp)

if(ULIB_FOUND AND IB_FOUND)
    include_directories(${ULIB_INCLUDE_DIRS} ${IB_INCLUDE_DIRS} ${ULIB_ROOTLINK_DIR})
 target_link_libraries(${PROJECT_NAME}
    ############################################################################
    ## ULIB MODULES ##
         mutomCore
         mutomMath
         mutomDetectors
         mutomVtk
         mutomRoot
    ############################################################################
         mutomIB
    ############################################################################
        ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES}
    ############################################################################
    #    MANUAL INLCLUDE EXTERNAL DEPENDENCIES
    #    ${Boost_SERIALIZATION_LIBRARY}
    #    ${Boost_SIGNALS_LIBRARY}
    #    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    #    ${Geant4_LIBRARIES}
    #    ${ROOT_LIBRARIES}
    #    ${VTK_LIBRARIES} # all VTK libs
    ############################################################################
    )
    target_link_libraries(IB_abtrim
        mutomCore
        mutomMath
        mutomDetectors
        mutomVtk
   ############################################################################
        mutomIB
   ############################################################################
    )
endif(ULIB_FOUND AND IB_FOUND)


