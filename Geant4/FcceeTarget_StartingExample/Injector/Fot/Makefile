OS        := $(shell uname)

ifeq ($(OS),Linux)
LIBEXT=so
endif
ifeq ($(OS),Darwin)
LIBEXT=dylib
endif
LIBNAME=libfot.$(LIBEXT)

ifeq ($(OS),Darwin)
LDFLAGS=-dynamiclib
endif

ifeq ($(OS),Linux)
LDFLAGS= -shared -g 
endif

TESTNAME = channeling
CFLAGS = -Ofast -Wall -Wno-unused-private-field -fPIC

INCLUDES =
LIBS =
OBJ_FILES = Crystal.o Evenement.o mathematics.o statistiques.o Bremsstrahlung.o ParticleCollection.o Photon.o Particle.o PhotonCollection.o Snake.o Fot.o

all : $(LIBNAME) $(TESTNAME)

#$(EXEC_NAME) : $(OBJ_FILES)
#	$(CXX) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)

$(LIBNAME) : $(OBJ_FILES)
	$(CXX) $(LDFLAGS) -o $(LIBNAME) $(OBJ_FILES) $(LIBS)

$(TESTNAME) : test_channeling.o $(OBJ_FILES)
	$(CXX) `root-config --glibs` -L. -lfot -o $(TESTNAME) $< 

test_channeling.o: test_channeling.cc
	$(CXX) $(CFLAGS) $(INCLUDES) `root-config --cflags` -o test_channeling.o -c $<

%.o: %.cc
	$(CXX) $(CFLAGS) $(INCLUDES) -o $@ -c $<

clean :
	rm -f test_channeling.o $(OBJ_FILES) $(LIBNAME) $(TESTNAME)
