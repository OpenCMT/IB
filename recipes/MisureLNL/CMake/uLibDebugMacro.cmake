include(uLibMacros)
include_guard(ULIB_DEBUG_MACRO_CMAKE)


set(CMAKE_CONFIGURE_DEBUG_ENABLE ON CACHE BOOL "cmake configure debugger")
set(CMAKE_CONFIGURE_DEBUG_LIST OFF CACHE BOOL "cmake debugger show list enable")
mark_as_advanced(
 CMAKE_CONFIGURE_DEBUG_ENABLE
 CMAKE_CONFIGURE_DEBUG_LIST
)

macro(debug str)
 if(CMAKE_CONFIGURE_DEBUG_ENABLE)
 set(var ${${str}})
 list(LENGTH var len)
 if((${len} GREATER 1) AND CMAKE_CONFIGURE_DEBUG_LIST)
  message(STATUS "[DEBUG] [${str}] list: ")
  foreach(item ${var})
  message(STATUS "        -> ${item}")
  endforeach()
 else()
  message(STATUS "[DEBUG] [${str}] ${var}")
 endif()
 endif()
endmacro()

macro(debug_list str)
 if(CMAKE_CONFIGURE_DEBUG_ENABLE)
 set(var ${${str}})
 list(LENGTH var len)
 if((${len} GREATER 1))
  message(STATUS "[DEBUG] [${str}] list: ")
  foreach(item ${var})
  message(STATUS "        -> ${item}")
  endforeach()
 else()
  message(STATUS "[DEBUG] [${str}] ${var}")
 endif()
 endif()
endmacro()

macro(debug_package str)
 if(CMAKE_CONFIGURE_DEBUG_ENABLE)
  set(have ${str}_FOUND)
  if(${have})
   debug("${str}_INCLUDE_DIRS")
   debug("${str}_LIBRARIES")
   debug("${str}_LIBRARY_DIRS")
   debug("${str}_DEFINITIONS")
  else(${have})
   message(WARNING "package ${str} not found")
  endif(${have})
 endif(CMAKE_CONFIGURE_DEBUG_ENABLE)
endmacro()
