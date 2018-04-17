# UTILS
set( UTILS
        IB_lnlROCBuilder
        IB_lnlROCBuilder2
        IB_portal3x6x5ROCBuilder
        IB_ROCScaleToMidpt
        IB_ROCScaleToMidAwo
        IB_ROCScaleToMideff
        IB_ROCPlot
)

set(LIBRARIES
       ${PACKAGE_LIBPREFIX}Core
       ${PACKAGE_LIBPREFIX}Math
#       ${PACKAGE_LIBPREFIX}Vtk
       ${PACKAGE_LIBPREFIX}IB
)

uLib_add_utils(IB-ROC-utils)
