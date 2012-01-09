# Makefile for the ROOT RADMU programs.
# Expectation Maximization - EM
# Author: Sara Vanini 10/6/2008

ObjSuf        = o
SrcSuf        = cxx
ExeSuf        = x
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

# ROOT
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

#Valgrind client
GRIND_CLIENT = $(shell pkg-config --libs --cflags valgrind)

# Linux
# add -m32 to compile and link at 32 bit, otherwise incompatible with Hit Bank...
CXX           = g++ -m64
CC            = gcc -m64
#flags for maximum optimization (but no debugging....) 
CXXFLAGS      = -Wall -Winline -fPIC -I$(ROOTSYS)/include
CXXFLAGS     += $(ROOTCFLAGS)
LD            = g++ 
LDFLAGS       = -m64
#SOFLAGS       = -Wl,-soname,libEvent.so -shared
LIBS          = $(ROOTLIBS) -lgcc -lm -ldl
GLIBS         = $(ROOTGLIBS)


OPT =  -DNDEBUG -O2 -fopenmp #-falign-function=2 -falign-jumps=2
ifdef mode
ifeq ($(mode),debug)
	OPT = -g -D_DEBUG 
else
	OPT =  -DNDEBUG -O2 -fopenmp #-falign-function=2 -falign-jumps=2
endif
endif

#------------------------------------------------------------------------------
OBJvEM = templates.o debug.o ImgAnalyzer.o MuonCollection.o VoxCollection.o Muon.o Voxel.o IBMuon.o IBMuonCollection.o IBGeometry.o IBVoxCollection.o IBAnalyzerEM.o
OBJrunEM = templates.o debug.o IBMuon.o IBMuonCollection.o IBGeometry.o IBGrid3d.o IBVoxCollection.o IBAnalyzerEM.o IBAnalyzerEM_simple.o IBAnalyzerPoca.o IBVoxFilters.o


runEM:  ${OBJrunEM} main.C
	${CXX} main.C -o runEM ${CXXFLAGS} ${OBJrunEM} $(OPT) $(LIBS) $(GLIBS)

lib:    ${OBJrunEM}
	$(CXX)  $(CXXFLAGS) $(OPT) $(LIBS) $(GLIBS) -shared -Wl,-soname,libMutomIB.so -o libMutomIB.so $(OBJrunEM)

test_debug: ${OBJrunEM} test_debug.C
	${CXX} test_debug.C -o test_debug ${CXXFLAGS} ${OBJrunEM} $(OPT) $(LIBS) $(GLIBS)

test_timescale: ${OBJrunEM} test_timescale.C
	${CXX} test_timescale.C -o test_timescale ${CXXFLAGS} ${OBJrunEM} $(OPT) $(LIBS) $(GLIBS)

test_muons: ${OBJrunEM} test_muons.C
	${CXX} test_muons.C -o test_muons ${CXXFLAGS} ${OBJrunEM} $(OPT) $(LIBS) $(GLIBS)

test_filters: ${OBJrunEM} test_filters.C
	${CXX} test_filters.C -o test_filters ${CXXFLAGS} ${OBJrunEM} $(OPT) $(LIBS) $(GLIBS)

test_Sijcut: ${OBJrunEM} test_Sijcut.C
	${CXX} test_Sijcut.C -o test_Sijcut ${CXXFLAGS} ${OBJrunEM} $(OPT) $(LIBS) $(GLIBS)

test_gauss: ${OBJrunEM} test_gauss.C
	${CXX} test_gauss.C -o test_gauss ${CXXFLAGS} ${OBJrunEM} $(OPT) $(LIBS) $(GLIBS)

vEM: ${OBJvEM} main_v.C
	${CXX} main_v.C -o vEM ${CXXFLAGS} ${OBJvEM} $(OPT) $(LIBS) $(GLIBS)


doxy:
	doxygen doxy.conf

# ---- EM classes -----
#Voxel.o: Voxel.C Voxel.h
#	${CXX} ${CXXFLAGS} -c Voxel.C

#Muon.o: Muon.C Muon.h
#	${CXX} ${CXXFLAGS} -c Muon.C
#
#VoxCollection.o: VoxCollection.C VoxCollection.h
#	${CXX} ${CXXFLAGS} -c VoxCollection.C
#
#MuonCollection.o: MuonCollection.C MuonCollection.h
#	${CXX} ${CXXFLAGS} -c MuonCollection.C
#
#ImgAnalyzer.o: ImgAnalyzer.C ImgAnalyzer.h
#	${CXX} ${CXXFLAGS} -c ImgAnalyzer.C
#
#debug.o: debug.c
#	${CXX} ${CXXFLAGS} -c debug.c 

.C.o:
	$(CXX) $(CXXFLAGS) $(OPT) -c $<

.c.o:
	$(CC) $(CXXFLAGS) $(OPT) -c $<

# ---- THE MAIN -----
#main.o: main.C
#	${CXX} ${CXXFLAGS} $(GRIND_CLIENT) -c main.C

clean:
	rm -f runEM vEM libMutomIB.so test_muons test_timescale test_debug test_filters test_Sijcut ${OBJvEM} ${OBJrunEM} *~ *.vtk *.vti *.ply
	@echo "all cleaned up!"
