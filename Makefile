# MAKEFILE - Makefile for NMSAC library and main programs
#
# Author: Nils Maercklin, March 2010

# Directories for executables and libraries (change as required):
BINDIR = ./bin
LIBDIR = ./lib

# Compiler and linker flags:
CFLAGS = -O -Wall
LFLAGS = -lm -L$(LIBDIR) -lnmsac
ARFLAGS = rv

# Main programs:
PROGS =	sacapkbk sacpkaic sacbutter sacop1 sacmamp sacpolar sacrot3 \
	sacpofilt sacpolwh saclh sacsh sacpkcor

# Object files and NMSAC library file name:
LOBJ = butterworth.o nmutils.o nmsaclib.o nmpolar.o nmgeo.o nmetime.o
SACLIB = $(LIBDIR)/libnmsac.a


install:	$(BINDIR)	$(PROGS)
	-mv $(PROGS) $(BINDIR)
	-rm -f $(LOBJ)

all:	$(PROGS)

$(PROGS):	$(SACLIB)
	-$(CC) $(CFLAGS) $(@F).c $(LFLAGS) -o $@

$(SACLIB):	$(LIBDIR)	$(LOBJ)
	$(AR) $(ARFLAGS) $(SACLIB) $(LOBJ);
	ranlib $(SACLIB)

$(BINDIR):
	if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi

$(LIBDIR):
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi

clean:
	-rm -f $(PROGS) $(LOBJ)

cleanall:
	-rm -f $(PROGS) $(LOBJ) $(SACLIB)
	cd $(BINDIR); rm -f $(PROGS)

remake:
	-rm -f $(PROGS) $(SACLIB) $(LOBJ)
	$(MAKE)

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c $<

