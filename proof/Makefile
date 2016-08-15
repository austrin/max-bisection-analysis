CXXFLAGS=-O3 -Wall
LDFLAGS=-lmpfr

MAIN=prover
OBJECTS=alpha.o bvnl.o configurations.o gaussian.o search.o
HEADERS=intervals.h

all: $(MAIN)

%.o: %.cc %.h $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) -c -o $@ $<


$(MAIN): $(MAIN).cc $(OBJECTS) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) -o $@ $(MAIN).cc $(OBJECTS) $(LDFLAGS)

prover_linear: prover_linear.cc $(OBJECTS) $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) -o $@ prover_linear.cc $(OBJECTS) $(LDFLAGS)

clean: 
	rm -f $(OBJECTS) $(MAIN) prover_linear

