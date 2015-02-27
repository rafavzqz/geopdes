OCTFILES =$(patsubst %.cc,%.oct,$(wildcard *.cc)) 
CXXFLAGS = ${CXXFLAGS}
LDFLAGS  = ${LDFLAGS}

all: $(OCTFILES)

%.oct:  %.cc
	mkoctfile $(GEOPDES_INCLUDE_DIR) $<

clean:
	$(RM) *.o core octave-core *.oct *~
