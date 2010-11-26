OCTFILES =$(patsubst %.cc,%.oct,$(wildcard *.cc)) 

all: $(OCTFILES)

%.oct: %.cc
	mkoctfile $(GEOPDES_INCLUDE_DIR) $<

clean:
	$(RM) *.o core octave-core *.oct *~
