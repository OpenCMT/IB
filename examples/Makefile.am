
include $(top_srcdir)/Common.am

#AM_DEFAULT_SOURCE_EXT = .cpp

# if HAVE_CHECK
TESTS =  LNL_SijGuess \
         LNL_SijGuess_Trim \
         IBReaderTest \
         IBVoxCollectionTest \
         IBPoCATest \
         IBMinVarTest \
         IBMinVarTest2 \
         IBAnalyzerPocaTest \
         IBRaytracerTest \
         IBAnalyzerWPocaTest \
         IBAnalyzerEMTest \
         IBAnalyzerEMTestLNL \
         IBAnalyzerTrackCountTest \
         IBRaytracerEntryTest \
         IBVoxFilterTest \
         IBVoxCollectionThresholdTest \
         IBVtkReadAndFilter \
         IBQuickBuildTest \
         IBFilterUtility \
         IBAnalyzerEMTestMGA \
         IBBulkBuilder \
         IBVoxBoxTest \
         IBScanImageTest \
         IBRocBuilder \
         IBSingleROCBuilder \
         IBlnlROCBuilder \
         IBPocaDiff \
         IBMuonEventReader3DReaderTest	\
         IBImageFindGradient \
         IBAnalyzerEM_WeightMAP \
         IBImageApplyGradient \
         IBAnalyzerEMShiftTest \
         IBAnalyzerEMWeightFactor \
         IBAnalyzerEMShiftP \
         IBAnalyzerEMShakeTest \
         IBMuonCollectionTest \
         IBAnalyzerEM_NeighbourMAP\
         IBTotalDensityTest \
         IBPlaneRocBuilder \
         IBChamberDistanceTest \
         IBChi2CutStudy

# else
# TEST =
# endif

LDADD = $(top_srcdir)/libmutomIB-${PACKAGE_VERSION}.la

check_PROGRAMS = $(TESTS)

