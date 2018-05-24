
CXX ?= g++
CFLAGS = -c -g -Wall --std=c++11 -fPIC

AR = ar
ARFLAGS = -r

INSTALL = install
INSTALLFLAGS = -D

LIBS =

TARGETS = libtraveltime2d.a

OBJS = traveltimeexception.o

SRCS = Makefile \
	circularbuffer.hpp \
	coordinate.hpp \
	subfield.hpp \
	traveltimefield.hpp \
	traveltimeexception.hpp \
	traveltimeexception.cpp \
	traveltimenode.hpp \
	velocityfield.hpp \
	velocityfieldref.hpp 

all : $(TARGETS)

libtraveltime2d.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS) 

%.o : %.cpp
	$(CXX) $(CFLAGS) -o $*.o $*.cpp

DATE = $(shell date +"%Y%m%d%H%M")
DIR = traveltime2d
TGZ = $(DIR).tar.gz

dist :
	mkdir -p $(DIR)
	echo $(DATA) > $(DIR)/Version
	for f in $(SRCS) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean : 
	rm -f $(TARGETS) *.o
