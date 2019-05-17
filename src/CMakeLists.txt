set(HEADERS
                          IB.h
                          IBAbstractReader.h
                          IBAnalyzer.h
                          IBAnalyzerEM.h
                          IBAnalyzerEM3p.h
                          IBAnalyzerEMTrim.h
                          IBAnalyzerEMAlgorithm.h
                          IBAnalyzerEMAlgorithmMGA.h
                          IBAnalyzerEMAlgorithmSGA.h
                          IBAnalyzerPoca.h
                          IBAnalyzerTrackCount.h
                          IBAnalyzerTrackLengths.h
                          IBAnalyzerWPoca.h
                          IBAnalyzerWTrackLengths.h
                          IBExperiment.h
                          IBLineDistancePocaEvaluator.h
                          IBMAPUpdateDensityAlgorithms.h
                          IBMinimizationVariablesEvaluator.h
                          IBMuonCollection.h
                          IBMuonError.h
                          IBMuonEventTTreeLNLdataReader.h
                          IBMuonEventTTreeReader.h
                          IBNormalPlaneMinimizationVariablesEvaluator.h
                          IBSimpleTwoViewsMinimizationVariablesEvaluator.h
                          IBPWeightAlgorithms.h
                          IBPocaEvaluator.h
                          IBSubImageGrabber.h
                          IBTiltedAxisPocaEvaluator.h
                          IBVoxCollection.h
                          IBVoxCollectionCap.h
                          IBVoxFilters.h
                          IBVoxImageScanner.h
                          IBVoxRaytracer.h
                          IBVoxel.h
                          IBVoxImageFilterPlasmon.hpp
                          IBVoxImageFilterGradient.hpp
                          IBROC.h
                          IBROCBuilder.h
)

SET(SOURCES
                IB.cpp
                IBMinimizationVariablesEvaluator.cpp
                IBNormalPlaneMinimizationVariablesEvaluator.cpp
                IBSimpleTwoViewsMinimizationVariablesEvaluator.cpp
                IBAnalyzer.cpp
                IBPocaEvaluator.cpp
                IBTiltedAxisPocaEvaluator.cpp
                IBLineDistancePocaEvaluator.cpp
                IBAnalyzerPoca.cpp
                IBAnalyzerEM.cpp
                IBAnalyzerEM3p.cpp
                IBAnalyzerEMTrim.cpp
                IBMuonEventTTreeReader.cpp
                IBMuonEventTTreeLNLdataReader.cpp
                IBAnalyzerWPoca.cpp
                IBVoxCollection.cpp
                IBVoxFilters.cpp
                IBAnalyzerEMAlgorithm.cpp
                IBAnalyzerEMAlgorithmSGA.cpp
                IBAnalyzerEMAlgorithmMGA.cpp
                IBAnalyzerTrackCount.cpp
                IBAnalyzerTrackLengths.cpp
                IBAnalyzerWTrackLengths.cpp
                IBMuonCollection.cpp
                IBMuonError.cpp
                IBROCBuilder.cpp
)



set(LIBRARIES
       ${PACKAGE_LIBPREFIX}Core
       ${PACKAGE_LIBPREFIX}Math
       ${PACKAGE_LIBPREFIX}Root
)

# DEPRECATED
#uLib_add_shared_library(IB)
add_library(${PACKAGE_LIBPREFIX}IB SHARED ${SOURCES})

source_group("${PROJECT_NAME}\\IB" FILES ${SOURCES})

set_target_properties(${PACKAGE_LIBPREFIX}IB PROPERTIES
		                  VERSION ${PROJECT_VERSION}
		                  SOVERSION ${PROJECT_SOVERSION}
#		                  INTERFACE_LINK_LIBRARIES "${LIBRARIES}"
#		                  LINK_INTERFACE_LIBRARIES "${LIBRARIES}"          ## <- ok for cmake from  2.12 ?
#		                  IMPORTED_LINK_INTERFACE_LIBRARIES "${LIBRARIES}" ## <- breaks cmake after 2.11 ?
)

target_link_libraries(${PACKAGE_LIBPREFIX}IB ${LIBRARIES})

install(TARGETS ${PACKAGE_LIBPREFIX}IB
        EXPORT "${PROJECT_NAME}Targets"
        RUNTIME DESTINATION ${PACKAGE_INSTALL_BIN_DIR} COMPONENT bin
        LIBRARY DESTINATION ${PACKAGE_INSTALL_LIB_DIR} COMPONENT lib
)

install(FILES ${HEADERS} DESTINATION ${PACKAGE_INSTALL_INC_DIR}/IB)

#add_subdirectory(testing)
