
TDTBASE=../../TDTbase

INCLUDES = -I.. \
	-I$(TDTBASE)/tracking

CXX = g++
CXXFLAGS = -c -O3 -Wall -std=c++11 $(INCLUDES)

#CXX = /opt/intel/bin/icc
#CXXFLAGS = -c -O3 $(INCLUDES)

LIBS = -L$(TDTBASE)/tracking -l tracking  -L.. -ltraveltime2d

TARGETS = timing



all : $(TARGETS)

timing : timing.o
	$(CXX) -o timing timing.o $(LIBS)

timing.o : ../traveltimefield.hpp ../coordinate.hpp ../velocityfield.hpp ../traveltimenode.hpp

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

clean :
	rm -f $(TARGETS) *.o
