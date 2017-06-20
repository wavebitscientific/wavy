# wavy top-level Makefile for external dependencies
#
.PHONY: all clean

DATETIME   = src/external/datetime-fortran
FUNCTIONAL = src/external/functional-fortran
JSON       = src/external/json-fortran

all: datetime functional json

datetime: 
	mkdir $(DATETIME)/build
	cd $(DATETIME)/build && FC=${FC} cmake .. && make && ctest

functional: 
	mkdir $(FUNCTIONAL)/build
	cd $(FUNCTIONAL)/build && FC=${FC} cmake .. && make && ctest

json: 
	cp src/external/pickFortranCompilerFlags.cmake $(JSON)/cmake
	mkdir $(JSON)/build_real32
	mkdir $(JSON)/build_real64
	mkdir $(JSON)/build_real128
	cd $(JSON)/build_real32 && FC=${FC} cmake .. -DSKIP_DOC_GEN:BOOL=True -DREAL=32 && make
	cd $(JSON)/build_real64 && FC=${FC} cmake .. -DSKIP_DOC_GEN:BOOL=True -DREAL=64 && make
	cd $(JSON)/build_real128 && FC=${FC} cmake .. -DSKIP_DOC_GEN:BOOL=True -DREAL=128 && make

clean:
	rm -r {$(DATETIME),$(FUNCTIONAL),$(JSON)}/build*
