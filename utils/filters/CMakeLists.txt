# UTILS
set( UTILS
        IB_abtrim
        IB_abtrimu
        IB_bilateral
        IB_median
        IB_spr
        IB_scale
)

set(LIBRARIES
       ${PACKAGE_LIBPREFIX}Core
       ${PACKAGE_LIBPREFIX}Math
       ${PACKAGE_LIBPREFIX}Vtk
       ${PACKAGE_LIBPREFIX}Root
       ${PACKAGE_LIBPREFIX}IB
)

uLib_add_utils(IB-filters)
