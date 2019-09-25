# UTILS
set(UTILS IB_abtrim
          IB_abtrimu
          IB_bilateral
          IB_median
          IB_spr
          IB_scale
)

set(LIBRARIES  ${PACKAGE_LIBPREFIX}Core
               ${PACKAGE_LIBPREFIX}Math
               ${PACKAGE_LIBPREFIX}Vtk
               ${PACKAGE_LIBPREFIX}IB
)

# DEPRECATED
#uLib_add_utils(IB-filters)

foreach(tn ${UTILS})

  add_executable(${tn} ${tn}.cpp)

  target_link_libraries(${tn} ${LIBRARIES})

  install(TARGETS ${tn} RUNTIME DESTINATION bin)

endforeach(tn)

add_custom_target(all-IB-filters)
add_dependencies(all-IB-filters ${UTILS})

