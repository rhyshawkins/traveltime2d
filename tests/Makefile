
INCLUDES = -I..

CXX = g++
CXXFLAGS = -c -g -Wall -std=c++11 $(INCLUDES)

CXXFLAGS += -O3

LIBS = -lgtest -L.. -ltraveltime2d -lpthread

TARGETS = velocityfield_tests \
	coordinate_tests \
	traveltimefield_tests \
	iceland_traveltime_tests \
	subfield_tests \
	pseudopath_tests \
	gradient_tests \
	anomally_tests


all : $(TARGETS)

velocityfield_tests : velocityfield_tests.o
	$(CXX) -o velocityfield_tests velocityfield_tests.o $(LIBS)

velocityfield_tests.o : ../velocityfield.hpp

coordinate_tests : coordinate_tests.o
	$(CXX) -o coordinate_tests coordinate_tests.o $(LIBS)

coordinate_tests.o : ../coordinate.hpp

traveltimefield_tests : traveltimefield_tests.o
	$(CXX) -o traveltimefield_tests traveltimefield_tests.o $(LIBS)

gradient_tests : gradient_tests.o
	$(CXX) -o gradient_tests gradient_tests.o $(LIBS)

pseudopath_tests : pseudopath_tests.o
	$(CXX) -o pseudopath_tests pseudopath_tests.o $(LIBS)

traveltimefield_tests.o : ../traveltimefield.hpp ../coordinate.hpp ../velocityfield.hpp ../traveltimenode.hpp

iceland_traveltime_tests : iceland_traveltime_tests.o
	$(CXX) -o iceland_travelime_tests iceland_traveltime_tests.o $(LIBS)

iceland_traveltime_tests.o : ../traveltimefield.hpp ../coordinate.hpp ../velocityfield.hpp ../traveltimenode.hpp

subfield_tests : subfield_tests.o
	$(CXX) -o $@ $< $(LIBS)

anomally_tests : anomally_tests.o
	$(CXX) -o $@ $< $(LIBS)

subfield_tests.o : ../subfield.hpp 

pseudopath_tests.o : ../traveltimefield.hpp ../coordinate.hpp ../velocityfield.hpp ../traveltimenode.hpp

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

clean :
	rm -f $(TARGETS) *.o
