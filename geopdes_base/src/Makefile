OCTFILES =$(patsubst %.cc,%.oct,$(wildcard *.cc)) 
CXXFLAGS =-fopenmp ${CXXFLAGS}
LDFLAGS  =-fopenmp ${LDFLAGS}

all: $(OCTFILES)

%.oct:  %.cc
	mkoctfile $(GEOPDES_INCLUDE_DIR) $<

clean:
	$(RM) *.o core octave-core *.oct *~
