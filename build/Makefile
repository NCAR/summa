#========================================================================
# Makefile to compile SUMMA
#========================================================================
#
# Recommended use: Copy this file to Makefile.local, edit it to your 
# heart's content, and then run `make -f build/Makefile.local` from
# your top level SUMMA directory. Don't include the Makefile.local in
# any pull requests you make.
#
# Note that Makefile configurations that we commonly use can be found on
# the SUMMA wiki at: 
# https://github.com/NCAR/summa/wiki/SUMMA-Makefile-Part-0-configuration
# feel free to add yours to that page.
#
# To troubleshoot your paths and setup, type 'make check'
#
# At a minimum you will need to set the following:
#  * F_MASTER         - top level summa directory
#  * FC               - compiler suite
#  * FC_EXE           - compiler executable
#  * INCLUDES         - path to include files
#  * LIBRARIES        - path to and libraries to include
# 
# Some further options can be specified for OpenMP, etc. See in Part 0 of
# the Makefile. You don't need to make any changes in PART 1 and 
# following unless you are doing SUMMA development and changed what 
# needs to be compiled

#========================================================================
# PART 0: User-configurable part
#========================================================================

# The variables can be specified in one of two ways:
# * delete the '##' in front of the variable, fill out the entry, 
#   save the file and run make
# * make no changes to this file, but specify the variables in your
#   environment before you run make

# Define core directory below which everything resides. This is the
# parent directory of the 'build' directory
##F_MASTER =

# Define the Fortran Compiler. If you are using gfortran, then this needs
# to be version 6 or higher. This variable is simply used to select the right
# compiler flags in the ifeq statements in this Makefile. The compiler
# executable is set separately as FC_EXE
# Currently this is either gfortran or ifort
##FC =

# Define the path for the compiler executable. This is the actual executable
# that is invoked. For example, FC=gfortran and FC_EXE=/usr/bin-gfortran-mp-6
# FC and FC_EXE have to be consistent
##FC_EXE =

# Define the NetCDF and LAPACK libraries and path to include files. 
# INCLUDES needs to be of the form (no quotes around the string):
# INCLUDES = -I<include path 1> -I<include path 2> -I<...> -I<include path n>
# LIBRARIES needs to be of the form ( no quotes around the string):
# LIBRARIES = '-L<netcdf library path> -lnetcdff -L<blas library path> -lblas -L<etc> -l<etc>'
# If none of this makes sense, please talk to your system
# administrator.
##INCLUDES =
##LIBRARIES =
       
# Eventually we plan move to a real configure script, but for now we like
# to keep track of successful compilations of SUMMA on different platforms
# and with different compilers. If you are successful compiling SUMMA,
# please add your configuration (operating system and compiler plus
# part 0 of the Makefile) to the SUMMA wiki on github.

# Define compiler flags. If you use a different compiler,
# you will need to figure out what the equivalent flags are
# and may need to update this section

# ------------ define compiler flags ----------------------------------------

# define open MP flags
isOpenMP =
FLAGS_OMP =
LIBOPENMP =

# Define compiler flags. If you use a different compiler,
# you will need to figure out what the equivalent flags are
# and may need to update this section

# gfortran compiler flags
ifeq "$(FC)" "gfortran"

	ifeq "$(isOpenMP)" "yes"
 		FLAGS_OMP = -fopenmp
	endif

# Production runs
FLAGS_NOAH = -O3 -ffree-form -ffree-line-length-none -fmax-errors=0 $(FLAGS_OMP)
FLAGS_COMM = -O3 -ffree-line-length-none -fmax-errors=0 $(FLAGS_OMP)
FLAGS_SUMMA = -O3 -ffree-line-length-none -fmax-errors=0 $(FLAGS_OMP)

# Debug runs
#FLAGS_NOAH = -p -g -ffree-form -ffree-line-length-none -fmax-errors=0 -fbacktrace -Wno-unused -Wno-unused-dummy-argument
#FLAGS_COMM = -p -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=bounds
#FLAGS_SUMMA = -p -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=bounds

endif

# ifort compiler flags
ifeq "$(FC)" "ifort"

	ifeq "$(isOpenMP)" "yes"
 		FLAGS_OMP = -qopenmp
	endif

# Production runs
FLAGS_NOAH = -O3 -noerror_limit -FR -auto -fltconsistency $(FLAGS_OMP)
FLAGS_COMM = -O3 -FR -auto -fltconsistency -fpe0 $(FLAGS_OMP)
FLAGS_SUMMA = -O3 -FR -auto -fltconsistency -fpe0 $(FLAGS_OMP)

# Debug runs
#FLAGS_NOAH = -O0 -p -g -warn nounused -noerror_limit -FR -auto -WB -traceback -fltconsistency
#FLAGS_COMM = -O0 -p -g -debug -warn all -check all -FR -auto -WB -traceback -fltconsistency -fpe0
#FLAGS_SUMMA = -O0 -p -g -debug -warn all -check all -FR -auto -WB -traceback -fltconsistency -fpe0
endif

#========================================================================
# PART 1: Define directory paths
#========================================================================

# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)/build/source

# Location of the compiled modules
MOD_PATH = $(F_MASTER)/build

# Define the directory for the executables
EXE_PATH = $(F_MASTER)/bin

#========================================================================
# PART 2: Assemble all of the SUMMA sub-routines
#========================================================================

# Define directories
DRIVER_DIR = $(F_KORE_DIR)/driver
HOOKUP_DIR = $(F_KORE_DIR)/hookup
NETCDF_DIR = $(F_KORE_DIR)/netcdf
DSHARE_DIR = $(F_KORE_DIR)/dshare
NUMREC_DIR = $(F_KORE_DIR)/numrec
NOAHMP_DIR = $(F_KORE_DIR)/noah-mp
ENGINE_DIR = $(F_KORE_DIR)/engine

#  utilities
SUMMA_NRUTIL= \
		nrtype.f90 \
		f2008funcs.f90 \
		nr_utility.f90
NRUTIL = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_NRUTIL))

#
# Numerical recipes procedures
# NOTE: all numerical recipes procedures are now replaced with free versions
SUMMA_NRPROC= \
      expIntegral.f90 \
      spline_int.f90
NRPROC = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_NRPROC))

# Hook-up modules (set files and directory paths)
SUMMA_HOOKUP= \
      ascii_util.f90 \
      summaFileManager.f90

HOOKUP = $(patsubst %, $(HOOKUP_DIR)/%, $(SUMMA_HOOKUP))

# Data modules
SUMMA_DATAMS= \
		multiconst.f90 \
		var_lookup.f90 \
		data_types.f90 \
		globalData.f90 \
		flxMapping.f90 \
		get_ixname.f90 \
		popMetadat.f90 \
		outpt_stat.f90
DATAMS = $(patsubst %, $(DSHARE_DIR)/%, $(SUMMA_DATAMS))

# utility modules
SUMMA_UTILMS= \
		time_utils.f90 \
		mDecisions.f90 \
		snow_utils.f90 \
		soil_utils.f90 \
		updatState.f90 \
		matrixOper.f90
UTILMS = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_UTILMS))

# Model guts
SUMMA_MODGUT= \
		MODGUT = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_MODGUT))

# Solver
SUMMA_SOLVER= \
		vegPhenlgy.f90 \
		diagn_evar.f90 \
		stomResist.f90 \
		groundwatr.f90 \
		vegSWavRad.f90 \
		vegNrgFlux.f90 \
		ssdNrgFlux.f90 \
		vegLiqFlux.f90 \
		snowLiqFlx.f90 \
		soilLiqFlx.f90 \
		bigAquifer.f90 \
		computFlux.f90 \
		computHeatCap.f90 \
		computThermConduct.f90 \
		computResid.f90 \
		computJacob.f90 \
		eval8summa.f90 \
		summaSolve.f90 \
		systemSolv.f90 \
		varSubstep.f90 \
		opSplittin.f90 \
		coupled_em.f90 \
		run_oneHRU.f90 \
		run_oneGRU.f90
SOLVER = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_SOLVER))

# Define routines for SUMMA preliminaries
SUMMA_PRELIM= \
		conv_funcs.f90 \
		sunGeomtry.f90 \
		convE2Temp.f90 \
		allocspace.f90 \
		checkStruc.f90 \
		childStruc.f90 \
		ffile_info.f90 \
		read_attrb.f90 \
		read_pinit.f90 \
		pOverwrite.f90 \
		read_param.f90 \
		paramCheck.f90 \
		check_icond.f90
PRELIM = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_PRELIM))

SUMMA_NOAHMP= \
		module_model_constants.F \
		module_sf_noahutl.F \
		module_sf_noahlsm.F \
		module_sf_noahmplsm.F
NOAHMP = $(patsubst %, $(NOAHMP_DIR)/%, $(SUMMA_NOAHMP))

# Define routines for the SUMMA model runs
SUMMA_MODRUN = \
      indexState.f90 \
      getVectorz.f90 \
      t2enthalpy.f90 \
      updateVars.f90 \
		var_derive.f90 \
		read_force.f90 \
		derivforce.f90 \
		snowAlbedo.f90 \
		canopySnow.f90 \
		tempAdjust.f90 \
		snwCompact.f90 \
		layerMerge.f90 \
		layerDivide.f90 \
		volicePack.f90 \
		qTimeDelay.f90
MODRUN = $(patsubst %, $(ENGINE_DIR)/%, $(SUMMA_MODRUN))

# Define routines for the solver
SUMMA_MSOLVE = \

# Define NetCDF routines
SUMMA_NETCDF = \
		netcdf_util.f90 \
		def_output.f90 \
		modelwrite.f90 \
		read_icond.f90
NETCDF = $(patsubst %, $(NETCDF_DIR)/%, $(SUMMA_NETCDF))

# ... stitch together common programs
COMM_ALL = $(NRUTIL) $(NRPROC) $(HOOKUP) $(DATAMS) $(UTILMS)

# ... stitch together SUMMA programs
SUMMA_ALL = $(NETCDF) $(PRELIM) $(MODRUN) $(SOLVER)

# Define the driver routine
SUMMA_DRIVER= \
		summa_type.f90 \
		summa_util.f90 \
		summa_alarms.f90 \
		summa_globalData.f90 \
		summa_defineOutput.f90 \
		summa_init.f90 \
		summa_setup.f90 \
		summa_restart.f90 \
		summa_forcing.f90 \
		summa_modelRun.f90 \
		summa_writeOutput.f90 \
		summa_driver.f90
DRIVER = $(patsubst %, $(DRIVER_DIR)/%, $(SUMMA_DRIVER))

# Define the executable
DRIVER__EX = summa.exe

# Define version number
VERSIONFILE = $(DRIVER_DIR)/summaversion.inc
VERSION = $(shell git tag | tail -n 1)
BULTTIM = $(shell date)
GITBRCH = $(shell git describe --long --all --always | sed -e's/heads\///')
GITHASH = $(shell git rev-parse HEAD)

#========================================================================
# PART 3: Checks
#======================================================================
# make sure that the paths are defined. These are just some high level checks
ifndef F_MASTER
 $(error F_MASTER is undefined)
endif
ifndef FC
 $(error FC is undefined: Specify your compiler)
endif
ifndef FC_EXE
 $(error FC_EXE is undefined: Specify your compiler executable)
endif
ifndef FLAGS_SUMMA
 $(error Specify flags for your compiler: $(FC))
endif
ifndef INCLUDES
 $(error INCLUDES is undefined)
endif
ifndef LIBRARIES
 $(error LIBRARIES is undefined)
endif

#========================================================================
# PART 4: compilation
#======================================================================

# Compile
all: compile_noah compile_comm compile_summa link clean install
part: compile_noah compile_comm compile_summa

check:
	$(info)
	$(info Displaying make variables:)
	$(info   F_MASTER   : $(F_MASTER))
	$(info   EXE_PATH   : $(EXE_PATH))
	$(info   FC         : $(FC))
	$(info   INCLUDES   : $(INCLUDES))
	$(info   LIBRARIES  : $(LIBRARIES))
	$(info   FLAGS_NOAH : $(FLAGS_NOAH))
	$(info   FLAGS_COMM : $(FLAGS_COMM))
	$(info   FLAGS_SUMMA: $(FLAGS_SUMMA))
	$(info)

# update version information
update_version:
	echo "character(len=64), parameter     :: summaVersion = '${VERSION}'" > $(VERSIONFILE)
	echo "character(len=64), parameter     :: buildTime = '${BULTTIM}'" >> $(VERSIONFILE)
	echo "character(len=64), parameter     :: gitBranch = '${GITBRCH}'" >> $(VERSIONFILE)
	echo "character(len=64), parameter     :: gitHash = '${GITHASH}'" >> $(VERSIONFILE)

# compile Noah-MP routines
compile_noah:
	$(FC_EXE) $(FLAGS_NOAH) -c $(NRUTIL) $(NOAHMP)

# compile common routines

compile_comm:
	$(FC_EXE) $(FLAGS_COMM) -c $(COMM_ALL) $(INCLUDES) 

# compile SUMMA routines
compile_summa: update_version
	$(FC_EXE) $(FLAGS_SUMMA) -c $(SUMMA_ALL) $(DRIVER) \
	$(INCLUDES)

# link routines
link:
	$(FC_EXE) -g *.o $(LIBRARIES) -o $(DRIVER__EX)

# Remove object files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f soil_veg_gen_parm__genmod.f90

# Copy the executable to the bin directory
install:
	@mkdir -p $(EXE_PATH)
	@mv $(DRIVER__EX) $(EXE_PATH)
	$(info $(DRIVER__EX) successfully installed in $(EXE_PATH))
