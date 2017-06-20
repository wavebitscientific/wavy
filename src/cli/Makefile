# wavy Makefile
#
.SUFFIXES: .f90 .o
.PHONY: all wavy

FC      = gfortran
FCFLAGS = -O0 -fbacktrace -C -g

INCLUDE = ../lib
LIB     = ../lib

JSONLIB = ../external/json-fortran-5.1.0/build/lib
JSONINC = ../external/json-fortran-5.1.0/build/include

FLAPLIB = ../external/flap-1.1.6/build/lib
FLAPINC = ../external/flap-1.1.6/build/modules

OBJS =

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) -I$(JSONINC) -I$(FLAPINC) $<

all: wavy-cli

wavy-cli: wavy-cli.o
	$(FC) $(FCFLAGS) $(OBJS) $^ -L$(LIB) -lwavy -L$(JSONLIB) -ljsonfortran -L$(FLAPLIB) -lFLAP -lPENF -o $@

test-cli: test-cli.o
	$(FC) $(FCFLAGS) $(OBJS) $^ -L$(LIB) -lwavy -L$(JSONLIB) -ljsonfortran -L$(FLAPLIB) -lFLAP -lPENF -o $@

wavy-cli.o: wavy-cli.f90
test-cli.o: test-cli.f90

clean:
	rm -f *.o *.mod wavy-cli
