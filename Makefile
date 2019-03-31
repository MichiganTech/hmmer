# You may want to modify the following make variables:
#   BINDIR  - where the executables will be installed by a 'make install'
#   MANDIR  - where the man pages will be installed by a 'make install'
#   CC      - which compiler to use
#   CFLAGS  - compiler flags to use
#
export prefix      = /usr
export exec_prefix = ${prefix}
export BINDIR      = ${exec_prefix}/bin
export MANDIR      = ${prefix}/share/man

export PROGSUFFIX = 2
export MANSUFFIX = 1

export CFLAGS += -pipe -Wall -Wextra -Wpedantic -Wno-sign-compare -fPIC -O0 -ggdb
export LIBS    = -lm -pthread

# The program lists below for HMMER are not necessarily
# a complete manifest. They are the list of stable programs that the
# package will install. There must be a man page for each one of them
# in the appropriate places (documentation/man for HMMER)
#
export PROGS = hmmalign\
               hmmbuild\
               hmmcalibrate\
               hmmconvert\
               hmmemit\
               hmmfetch\
               hmmindex\
               hmmpfam\
               hmmsearch


HMMER2_LIBS = libhmmer$(PROGSUFFIX).a


# all: Compile everything.
#
all:
	(cd src;       make; make module);\
	(cd testsuite; make);\


# Compiles and runs test suite in testsuite.
check:
	(cd testsuite; make check)

src/libhmmer.a:
	(cd src;   make; make module)


# install: installs the binaries in BINDIR/
#          installs man pages in MANDIR/man1/  (e.g. if MANSUFFIX is 1)
#          Creates these directories if they don't exist.
#          Prefix those paths with ${DESTDIR} (rarely used, usually null;
#          may be set on a make command line when building contrib RPMs).
install:
	mkdir -p $(BINDIR)

	for file in $(PROGS) ; do\
		cp src/$$file "$(BINDIR)/$$file""$(PROGSUFFIX)" ;\
	done
	for file in documentation/man/*.man ; do\
	   install -D $$file $(MANDIR)/man$(MANSUFFIX)/`basename $$file .man`.$(MANSUFFIX);\
	done

# Reverses the steps of "make install".  However, this should be handled
# managed by a package management system like apt, rpm, pacman, or similar.
uninstall:
	for file in $(PROGS) ; do\
	   rm "${BINDIR}/$$file""$(PROGSUFFIX)";\
	done
	for file in hmmer2 $(ls documentation/man/); do\
	   rm -f $(MANDIR)/man$(MANSUFFIX)/$$file.$(MANSUFFIX);\
	done

# after building, prep the whole directory for a binary distribution: symlink
# supported binaries in binaries subdir, remove everything but binaries and
# Makefiles.
bindist:
	mkdir binaries
	for prog in ${PROGS} do\
	   (cd binaries; ln -s ../src/$$prog .);\
	done

# "make clean" removes almost everything except configuration files.
clean:
	(cd src;       make clean)
	(cd testsuite; make clean)
	rm -f *.o *~ Makefile.bak core TAGS gmon.out
	rm -rf binaries


# doc:  build the Userguide and on-line manual
#
doc:
	(cd documentation/userguide; make)
