

MYLIBS	 =
MYLIBDIR =
MYINCDIR =

HMMLIB = libhmmer$(PROGSUFFIX).so


OBJS =	alphabet.o\
	algorithms.o\
	debug.o\
	emit.o\
	emulation.o\
	getopt.o\
	histogram.o\
	hmmio.o\
	infocontent.o\
	lsj_eweight.o\
	mathsupport.o\
	masks.o\
	matrices.o\
	misc.o\
	modelmakers.o\
	plan7.o\
	plan9.o\
	postprob.o\
	prior.o\
	tophits.o\
	trace.o\
	ugly.o\
	hsregex.o\
	vectorops.o\
	msa.o\
	squid.o\
	stockholm.o\
	selex.o\
	phylip.o\
	msf.o\
	a2m.o\
	ssi.o\
	alignio.o\
	file.o\
	clustal.o\
	sqio.o\
	weight.o\
	cluster.o\
	sre_math.o\
	stack.o\
	aligneval.o\
	stopwatch.o\
	sre_random.o\
	sre_string.o

HDRS =	config.h\
	funcs.h\
	globals.h\
	lsjfuncs.h\
	postprob.h\
	structs.h

.c.o:
	${CXX} ${CFLAGS} ${PTHREAD_CFLAGS} ${MYINCDIR} ${LIBS} -c $<

#################################################################
## Targets defining how to make HMMER executables.
##
all:	 $(PROGS)

${PROGS}: %: %.o ${OBJS}
	${CXX} ${CFLAGS} ${MYLIBDIR} -o $@ ${OBJS} ${MYLIBS} $@.o ${PTHREAD_LIBS} ${LIBS} $(LDFLAGS)


#################################################################
## Targets used in making HMMER module for testsuite compilation.
##
module: $(HMMLIB)

$(HMMLIB): $(OBJS)
	$(AR) libhmmer2.a $(OBJS)
	$(RANLIB) libhmmer2.a
	chmod 644 libhmmer2.a
	gcc -shared $(OBJS) -o $@


clean:
	rm -f *.o *~ Makefile.bak core $(PROGS) TAGS gmon.out *.so
	rm -f hmmalign hmmbuild hmmcalibrate hmmconvert hmmemit hmmfetch hmmindex hmmpfam hmmsearch
