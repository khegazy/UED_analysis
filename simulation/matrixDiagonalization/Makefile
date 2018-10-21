CC=g++
#CFLAGS=$(shell root-config --cflags) -I/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel6-gcc44-opt/include -I/reg/g/psdm/sw/releases/ana-current/Eigen/include
CFLAGS=-pthread -m64 -I/reg/g/psdm/sw/external/root/5.30.06-python2.7/x86_64-rhel6-gcc44-opt/include/root -I/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel6-gcc44-opt/include -I/reg/g/psdm/sw/releases/ana-current/Eigen/include
#LINKLIBS=$(shell root-config --libs) -L/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel6-gcc44-opt/lib -lfftw3 #-L${ROOTSYS}/lib -lHist -lMatrix -lTree -lNet -lProofPlayer -lProof -lRIO -lThread -lMathCore -lPhysics
LINKLIBS=-L/reg/g/psdm/sw/external/root/5.30.06-python2.7/x86_64-rhel6-gcc44-opt/lib -lCint -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMathMore -lThread -pthread -lm -ldl -rdynamic -L/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel6-gcc44-opt/lib -lfftw3 #-L${ROOTSYS}/lib -lHist -lMatrix -lTree -lNet -lProofPlayer -lProof -lRIO -lThread -lMathCore -lPhysics

DEP_HPP_FILES := $(wildcard *Class.h) $(wildcard *Struct.h)
DEP_CPP_FILES := $(wildcard *.cpp)
DEPS := $(DEP_HPP_FILES) $(DEP_CPP_FILES)
#CPP_FILES := $(filter-out plotFunctions.cpp,$(ACPP_FILES))
CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:%.cpp=%.o)
DEP_OBJ_FILES := $(DEP_CPP_FILES:%.cpp=%.o)
EXE_FILES := $(CPP_FILES:%.cpp=%.exe)

all: $(EXE_FILES)

%.o: %.cpp $(DEPS)
	$(CC) -std=c++0x -c -o $@ $< $(CFLAGS) $(LINKLIBS)

%.exe: %.o $(DEP_OBJ_FILES)
	$(CC) -o $@ $^ $(LINKLIBS)

clean:
	rm -f *.o *.exe
