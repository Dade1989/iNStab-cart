##############################################################################
# Rules
#
all: pasta_2d.exe
default: pasta_2d.exe

include Dependencies.inc
# Is Dependencies.inc available ?
Dependencies.inc:
	@echo "##################################################################"
	@echo BEFORE COMPILING, YOU SHOULD HAVE AN APPROPRIATE FILE
	@echo Dependencies.inc AVALAIBLE.
	@echo "##################################################################"
	@exit 1

release:
	@echo making release
	@rm -f *.mod *.o
	@cd .. && tar -czf ./backups/pasta_2d_$(NOW).tar.gz pasta_2d/

pasta_2d.exe: main.o $(OBJ) 
	$(F90LNK) $(F90FLAGS) main.o $(OBJ) $(LIBS) -o $(FOLDER)$(NAME)

clean:
	@echo cleaning
	@rm -f *.o *.mod

###############################################################################
# Define variables

NOW := $(shell date +"%c" | tr ' :' '__')

# Compilers
FC = mpif90
CC = mpicc

F90CMP = $(FC) -c
F90LNK = $(FC)
F90OPT = -O3 -fopenmp -g
F90DEB = #-g -Wall #-fcheck=all
F90FLAGS = $(F90OPT) $(F90DEB)

CCMP   = $(CC) -c
CLNK   = $(CC)
COPT   = -O3
CDEB   =
CFLAGS = $(COPT) $(CDEB)

# Include paths
F90INCDIR = -I/usr/include
CINCDIR   = -I libLOCA

# Libraries
LIBARPACK = -larpack
LIBLAPACK = -llapack
LIBBLAS   = -lblas
LIBMUMPS  = -ldmumps -lzmumps -lpord
LIBMPI    = -lmpi
LIBLOCA   = -L./libLOCA/ -lloca
LIBS      = $(LIBMUMPS) $(LIBMPI) $(LIBLOCA) $(LIBARPACK) $(LIBLAPACK) $(LIBBLAS)

FOLDER = ./
NAME   = pasta_2d

# FOLDER = ~/BTSync/Software/bin/

###############################################################################
# Rules

.SUFFIXES:
.SUFFIXES: .c .f90 .o

.f90.o:
	$(F90CMP) $(F90FLAGS) $(F90INCDIR) $<

.c.o:
	$(CCMP) $(CFLAGS) $(CINCDIR) $<
