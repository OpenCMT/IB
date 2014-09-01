
################################################################################
##### CMAKE LISTS ##############################################################
################################################################################

cmake_minimum_required (VERSION 2.6)

execute_process(COMMAND "clear")
message("//////////////////////////////////////////////////////////////////////")
message("/////////////////// CMAKE PJOJECT CONFIGURATION //////////////////////")
message("//////////////////////////////////////////////////////////////////////")

## -------------------------------------------------------------------------- ##

project(IB)

# The version number.
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 2)
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")
set(PROJECT_SOVERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")

set(PACKAGE_VERSION   ${PROJECT_VERSION})
set(PACKAGE_NAME      "mutom-${PROJECT_VERSION}" CACHE STRING "name of the package defined inside code (config.h)")
set(PACKAGE_LIBPREFIX "mutom" CACHE STRING "suffix for each library component name")
set(PACKAGE_URL "http:://mutom.pd.infn.it")
set(PACKAGE_AUTHOR "Andrea Rigoni Garola <andrea.rigoni@pd.infn.it>")

## IMPORT MODULES ----------------------------------------------------------- ##

if(NOT DEFINED IB_CMAKE_DIR)
execute_process(COMMAND uLib_config --cmake-modules
                OUTPUT_VARIABLE ULIB_CONFIG_COMMAND_OUTPUT)
message(STATUS "Adding of CMake directory: ")
message(STATUS "${ULIB_CONFIG_COMMAND_OUTPUT}")
set(IB_CMAKE_DIR "${PROJECT_BINARY_DIR}/CMake" CACHE PATH "CMake modules path")
endif()

set(CMAKE_MODULE_PATH ${IB_CMAKE_DIR} ${CMAKE_MODULE_PATH})
message(STATUS "Module path: ${CMAKE_MODULE_PATH}")


## MACROS ------------------------------------------------------------------- ##

include(uLibMacros)
include(uLibDebugMacro)
include(uLibTargetMacros)

## CONFIG ------------------------------------------------------------------- ##

include(uLibConfigHeader)


## GLOBALS ------------------------------------------------------------------ ##

include(uLibCommon)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ") ## Add flags here ##


## FIND PACKAGES ------------------------------------------------------------ ##

include(uLibFindDependencies)





## CONFIG FIND PACKAGE
set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}/lib/cmake")
unset(ULIB_CONFIG)
find_package(ULIB 0.2 CONFIG
             NAMES uLib
             PATH_SUFFIXES mutom-0.2
)

debug(ULIB_CONFIG)
debug(ULIB_CONSIDERED_CONFIGS)
debug(ULIB_CONSIDERED_VERSIONS)
debug(ULIB_INCLUDE_DIRS)
debug(ULIB_LIBRARIES)

if(ULIB_CONFIG)
 set(ULIB_FOUND true)
 include_directories(${ULIB_INCLUDE_DIRS})
endif()



## ADD LIBRARIES SUBDIRECTORIES --------------------------------------------- ##

include_directories(${SRC_DIR})
add_subdirectory(${SRC_DIR})


add_subdirectory(${PROJECT_SOURCE_DIR}/utils/filters)
add_subdirectory(${PROJECT_SOURCE_DIR}/utils/roc)

add_subdirectory(${PROJECT_SOURCE_DIR}/examples EXCLUDE_FROM_ALL)



## INSTALLS AND EXPORTS ----------------------------------------------------- ##

export(PACKAGE IB)
export(TARGETS ${ULIB_SHARED_LIBRARIES}
       FILE "${PROJECT_BINARY_DIR}/IBTargets.cmake"
      )

# Create the FooBarConfig.cmake and FooBarConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${PACKAGE_INSTALL_CMAKE_DIR}"
   "${PACKAGE_INSTALL_INC_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src" "${PROJECT_BINARY_DIR}")
configure_file(IBConfig.cmake.in
  "${PROJECT_BINARY_DIR}/IBConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${IB_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(IBConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/IBConfig.cmake" @ONLY)
# ... for both
configure_file(IBConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/IBConfigVersion.cmake" @ONLY)

# Install the FooBarConfig.cmake and FooBarConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/IBConfig.cmake"
  "${PROJECT_BINARY_DIR}/IBConfigVersion.cmake"
  DESTINATION "${PACKAGE_INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT "${PROJECT_NAME}Targets" DESTINATION
  "${PACKAGE_INSTALL_CMAKE_DIR}" COMPONENT dev)

#install(FILES ${PROJECT_SOURCE_DIR}/ulib.h DESTINATION ${PACKAGE_INSTALL_INC_DIR})



## -------------------------------------------------------------------------- ##



# TODO clean also subdirectories
add_custom_target(clean-cmake-files
   COMMAND ${CMAKE_COMMAND} -P ${IB_CMAKE_DIR}/clean-all.cmake
)




