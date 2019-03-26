# UTILS
set(UTILS IB_lnlROCBuilder
          IB_lnlROCBuilder2
          IB_portal3x6x5ROCBuilder
          IB_ROCScaleToMidpt
          IB_ROCScaleToMidAwo
          IB_ROCScaleToMideff
#          IB_ROCPlot             TODO fix source
)

set(LIBRARIES  ${PACKAGE_LIBPREFIX}Core
               ${PACKAGE_LIBPREFIX}Math
               ${PACKAGE_LIBPREFIX}IB
)

# DEPRECATED
#uLib_add_utils(IB-ROC-utils)

foreach(tn ${UTILS})

  add_executable(${tn} ${tn}.cpp)

  target_link_libraries(${tn} ${LIBRARIES})

  install(TARGETS ${tn} RUNTIME DESTINATION bin)

endforeach(tn)

add_custom_target(all-IB-ROC-utils)
add_dependencies(all-IB-ROC-utils ${UTILS})

