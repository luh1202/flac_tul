# -*- Makefile -*-
#
# Makefile for flac
#
# Author: Eh Tan <tan2@mail.utexas.edu>
#

## Execute "make" if making production run. Or "make debug=1" for debugging run.
##
## debug = 0: optimized build; 1: debugging build
## omp = 1: enable openmp support
## coverage = 1: enable code coverage analysis (used for debugging)
## gprof = 1: enable profiling (used for tuning optimization)


########################################################################
## Select default fortran 90 compiler
########################################################################

## GNU version >= 4.2, otherwise openmp won't work
F90 = gfortran

## Intel, tested with v11.1
#F90 = ifort

## PGI group
#F90 = pgf90

## Sun
#F90 = f90


########################################################################
## any customized flags
########################################################################
CFLAGS = -g -fno-align-commons
LDFLAGS =


########################################################################
## default command line options
########################################################################
debug = 0
omp = 1
coverage = 0
gprof = 0


########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

## for GNU gfortran >= 4.2
ifeq ($(F90), gfortran)
	ifeq ($(omp), 1)
		CFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif
	ifeq ($(debug), 0)
		CFLAGS += -O3 -ffast-math
	else
		CFLAGS += -O0 -Wall -fbounds-check -ftrapv -fbacktrace -fdump-core \
			-Waliasing -Walign-commons -Wampersand -Wintrinsic-shadow \
			-Wline-truncation -Wsurprising -Wunderflow \
			-finit-integer=-99 -finit-real=nan -frange-check -Werror
	endif

	ifeq ($(coverage), 1)
		CFLAGS += -fprofile-arcs -ftest-coverage
		LDFLAGS += -fprofile-arcs -ftest-coverage
	endif

	ifeq ($(gprof), 1)
		CFLAGS += -pg
		LDFLAGS += -pg
	endif
endif

## for Intel v11.1
ifeq ($(F90), ifort)
	CFLAGS += -assume byterecl
	ifeq ($(omp), 1)
		CFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif
	ifeq ($(debug), 0)
		CFLAGS += -fast
	else
	       CFLAGS += -O0 -no-opt-subscript-in-range -ftrapuv -check all -traceback
	       LDFLAGS +=
	endif
endif

## for PGI
ifeq ($(F90), pgf90)
	ifeq ($(omp), 1)
		CFLAGS += -mp
		LDFLAGS += -mp
	endif
        ifeq ($(debug), 0)
                CFLAGS += -O3 -fast -fastsse -Mnoframe -Mipa=fast
        else
               CFLAGS += -O0 -Mbounds -Mchkstk -traceback
               LDFLAGS +=
        endif
endif

## for Solaris
ifeq ($(F90), f90)
	ifeq ($(omp), 1)
		CFLAGS += -xopenmp=parallel
		LDFLAGS += -xopenmp
	endif

# Explanation of flags --
# -r8const: Promote single precision const to double precision.
# -C: Enable runtime subscript range checking
# -fpover=yes: Provide run-time overflow check during READ
# -xcheck=%all: Enable all run-time check
	ifeq ($(debug), 0)
		CFLAGS += -O3 -r8const -fast -xloopinfo -s -xipo
	else
		CFLAGS += -w3 -r8const -C -fpover=yes -xcheck=%all
	endif
endif


########################################################################
## list of source files
########################################################################

# module files have to be compiled before other f90 files
SRCS =	\
		bc_update.f90 \
		dt_adjust.f90 \
		dt_mass.f90 \
		flac.f90 \
		fl_injectheat.f90 \
		fl_move.f90 \
		fl_node.f90 \
		fl_rheol.f90 \
		fl_srate.f90 \
		fl_therm.f90 \
		functions.f90 \
		init_areas.f90 \
		init_bc.f90 \
		init_cord.f90 \
		init_phase.f90 \
		init_stress.f90 \
		init_temp.f90 \
		init_visc.f90 \
		matprops.f90 \
		outflac.f90 \
		par.f90 \
		read_more.f90 \
		read_params.f90 \
		rem_cord.f90 \
		remesh.f90 \
		rem_test.f90 \
		rh_elastic.f90 \
		rh_maxwell.f90 \
		rh_plastic.f90 \
		rh_visc.f90 \
		rmasses.f90 \
		rsflac.f90 \
		saveflac.f90 \
		setflac.f90 \
		user_ab.f90 \
		user_ito.f90

#MODULES = $(MODULE_SRCS:.f90=.mod)
#MODULE_OBJS = $(MODULE_SRCS:.f90=.o)
OBJS = $(SRCS:.f90=.o)

EXE = flac

## Is this a mercurial repository?
HAS_HG := $(shell hg log -r tip --template '{node}' 2>/dev/null)

########################################################################
## actions
########################################################################

.PHONY: all clean clean-data take-snapshot

all: $(EXE) 

#$(EXE): $(MODULE_OBJS) $(OBJS)
#	$(F90) $(LDFLAGS) $(MODULE_OBJS) $(OBJS) -o $@
$(EXE): $(OBJS)
	$(F90) $(LDFLAGS) $(MODULE_OBJS) $(OBJS) -o $@

$(OBJS): %.o : %.f90 params.inc arrays.inc precision.inc
	$(F90) $(CFLAGS) -c $<

take-snapshot:
	@# snapshot of the code for building the executable
	@echo Flags used to compile the code: > snapshot.diff
	@echo '  '  CXX=$(CXX) debug=$(debug) openmp=$(omp) >> snapshot.diff
	@echo '  '  PATH=`pwd ${PATH}` >> snapshot.diff
	@echo '  '  LD_LIBRARY_PATH=`pwd ${LD_LIBRARY_PATH}` >> snapshot.diff
ifneq (HAS_HG),)
	@echo >> snapshot.diff
	@echo >> snapshot.diff
	@echo '==== Summary of the code ====' >> snapshot.diff
	@hg summary >> snapshot.diff
	@echo >> snapshot.diff
	@echo >> snapshot.diff
	@echo '== Code modification (not checked-in) ==' >> snapshot.diff
	@hg diff >> snapshot.diff
	@echo >> snapshot.diff
	@echo >> snapshot.diff
	@echo '== Code modification (checked-in but not public) ==' >> snapshot.diff
	@hg log --patch -r "not public()" >> snapshot.diff
else
	@echo \'hg\' is not in path, cannot take code snapshot. >> snapshot.diff
endif



clean:
	@rm -f $(MODULES) $(MODULE_OBJS) $(OBJS) $(EXE)

clean-data:
	@rm -f *.0 *.rs sys.msg marker_init.asc vbc.s area.dat output.asc temp0.dat
