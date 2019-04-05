
LIBS      =  -lm -pthread
MYLIBS    = -lhmmer$(PROGSUFFIX)
CFLAGS   += -I../src

SHIVA = alignalign_test\
	evd_test\
	masks_test\
	parsingviterbi_test\
	tophits_test\
	trace_test\
	viterbi_exercise\
	weeviterbi_test

#######
## Targets defining how to make Shiva executables.
#######

.c.o:
	$(CXX) $(CFLAGS) $(PTHREAD_CFLAGS) $(DEFS) -I../src -c $<

all: 	$(SHIVA)

$(SHIVA): %: %.o
	$(CXX) $(CFLAGS) $(PTHREAD_CFLAGS) $(DEFS) -o $@ -I../src -L../src $@.o $(MYLIBS) $(PTHREAD_LIBS) $(LIBS)

#######
## `make check` actually runs the tests.
#######

check:
	LD_LIBRARY_PATH=LD_LIBRARY_PATH:../src/ ./sqc 2 exercises.sqc . ../src


#######
## Miscellaneous
#######

clean:
	-rm -f *.o *~ Makefile.bak core $(SHIVA) TAGS gmon.out

TAGS:
	etags -t *.c *.h Makefile.in
