################################################################################
#list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/recipes/MisureLNL/CMake")

### MUDULE FIND PCL PACKAGE ####################################################

find_package(PCL 1.3 REQUIRED COMPONENTS common io)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

################################################################################

#unset(ULIB_CONFIG)
#unset(ULIB_INCLUDE_DIRS)
#find_package(ULIB 0.2 CONFIG REQUIRED
#             NAMES uLib
#             PATHS /usr/lib64/cmake /usr/lib/cmake
#             PATH_SUFFIXES mutom-0.2
#)
#if(NOT ULIB_CONFIG)
# message("Cannot detect uLib")
#endif()


add_executable(misureLNL main.cpp)

include_directories(${ULIB_INCLUDE_DIRS} ${IB_INCLUDE_DIRS} ${ULIB_ROOTLINK_DIR})

target_link_libraries( misureLNL
                       ${PACKAGE_LIBPREFIX}Core
                       ${PACKAGE_LIBPREFIX}Math
                       ${PACKAGE_LIBPREFIX}Vtk
                       ${PACKAGE_LIBPREFIX}IB
                       ${PCL_COMMON_LIBRARIES}
                       ${PCL_IO_LIBRARIES}
)
  ############################################################################
  #    MANUAL INLCLUDE EXTERNAL DEPENDENCIES
  #    ${Boost_SERIALIZATION_LIBRARY}
  #    ${Boost_SIGNALS_LIBRARY}
  #    ${Boost_PROGRAM_OPTIONS_LIBRARY}
  #    ${Geant4_LIBRARIES}
  #    ${ROOT_LIBRARIES}
  #    ${VTK_LIBRARIES} # all VTK libs
  ############################################################################

install(TARGETS misureLNL RUNTIME DESTINATION bin)

