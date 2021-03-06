#-------------------------------------------------------------------------------
# Makefile generated by PELICANS depend utility
# USER DEFINED PART : START
#-------------------------------------------------------------------------------
PELICANSHOME := /home/leo/Pelican/pelicans
LIBDIR        = $(PELICANSHOME)/lib/Linux-icc
BINDIR        = /global/home/leo/Pelican/Pel3.1_FluidSolver/lib/Linux-icc/opt2/
OPT           = $(flagsopt2)
CPPFLAGS      = -DLEVEL=2

MAKE_PEL = 0
LINK_PEL = 0
WITH_PROFILE = 0
WITH_COVERAGE = 0
# runtime consistency : MAKE_PEL=1 => LINK_PEL=1
ifeq ($(MAKE_PEL),1)
  LINK_PEL = 1
endif

SRC = 
SRC += $(wildcard /global/home/leo/Pelican/Pel3.1_FluidSolver/SecondOrderAdvectionDiffusion/src/*.cc)
SRC += $(wildcard /global/home/leo/Pelican/Pel3.1_FluidSolver/PELbase/src/*.cc)
SRC += $(wildcard /global/home/leo/Pelican/Pel3.1_FluidSolver/NavierStokes/src/*.cc)
SRC += $(wildcard /global/home/leo/Pelican/Pel3.1_FluidSolver/PDEsolver/src/*.cc)
SRC += $(wildcard /global/home/leo/Pelican/Pel3.1_FluidSolver/FrameFE/src/*.cc)
SRC += $(wildcard /global/home/leo/Pelican/Pel3.1_FluidSolver/RefSol/src/*.cc)

PRECOMP_OBJ = 

INC = 
INC += /home/leo/Pelican/pelicans/include
INC += /global/home/leo/Pelican/Pel3.1_FluidSolver/FrameFE/include
INC += /global/home/leo/Pelican/Pel3.1_FluidSolver/NavierStokes/include
INC += /global/home/leo/Pelican/Pel3.1_FluidSolver/PDEsolver/include
INC += /global/home/leo/Pelican/Pel3.1_FluidSolver/PELbase/include
INC += /global/home/leo/Pelican/Pel3.1_FluidSolver/RefSol/include
INC += /global/home/leo/Pelican/Pel3.1_FluidSolver/SecondOrderAdvectionDiffusion/include

LIBPATH= 
LIBPATH += $(LIBDIR)

LDLIBS=
LDLIBS += -lpel1
#-------------------------------------------------------------------------------
# USER DEFINED PART : END
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# EXTRA COMPONENT PART : START
#-------------------------------------------------------------------------------
# Inlining /home/leo/Pelican/pelicans/etc/extra-Linux.mak
# Common Makefile for Linux plateform runnning gcc
# Settings for Extra Components

ifeq ($(LINK_PEL),1)

WITH_MPI      = 0

##### PETSc external API
# if PETSc is enabled, uncomment the appropriate definition of PETSc_VERSION
WITH_PETSc    = 0
# PETSc_VERSION=2.3.0
# PETSc_VERSION=2.3.1
# PETSc_VERSION=2.3.3

WITH_AZTEC    = 0
WITH_SPARSKIT = 0
WITH_SILO     = 0
WITH_OPENGL   = 0
WITH_UMFPACK  = 1
WITH_METIS    = 0
WITH_MUMPS    = 0


endif

WITH_ZLIB = 1
WITH_EXTENDED_MATH = 1
WITH_JNI = 0
EXTRA_LIBS_DIR=/home/semar/pelican/ExternalPackages/etch

###############################################################################
ifeq ($(WITH_MUMPS),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/MUMPS/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/MUMPS/include
endif

MUMPS_DIR = $(EXTRA_LIBS_DIR)/mumps

SCALAP = -L$(EXTRA_LIBS_DIR)/scalapack/lib -lscalapack  -lblacs -lblacsF77 -lblacsC  -lblacs -lblacsF77 -lblacsC

#LIBBLAS = -lblas -llapack -lmpi_f77 -lmpi_f90
LIBBLAS = -lblas 

CPPFLAGS += -I$(MUMPS_DIR)/include 
LIBPATH  += $(MUMPS_DIR)/lib
LDLIBS   +=  -ldmumps -lmumps_common ${SCALAP} ${LIBBLAS}  -lpord  

endif

###############################################################################
ifeq ($(WITH_PETSc),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/PETSc_$(PETSc_VERSION)/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/PETSc_$(PETSc_VERSION)/include
endif

ifeq ($(PETSc_VERSION),2.3.0)
PETSc_DIR=/home/semar/pelican/LA_Packages/petsc/petsc-2.3.0
else
ifeq ($(PETSc_VERSION),2.3.1)
PETSc_DIR=/home/semar/pelican/LA_Packages/petsc/petsc-2.3.1-p0
else
ifeq ($(PETSc_VERSION),2.3.3)
PETSc_DIR=$(EXTRA_LIBS_DIR)/petsc
endif
endif
endif

PETSc_ARCH=etch
CPPFLAGS += -I$(PETSc_DIR)/include -I$(PETSc_DIR)/bmake/$(PETSc_ARCH)
LIBPATH  += $(PETSc_DIR)/lib/$(PETSc_ARCH)
LDLIBS   += -lpetscsnes -lpetscksp -lpetscmat -lpetscvec -lpetsc -lpetscdm 

WITH_X11 = 1
WITH_BLAS = 1
endif


###############################################################################
ifeq ($(WITH_AZTEC),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/Aztec/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/Aztec/include
endif

AZTECPATH = $(EXTRA_LIBS_DIR)/aztec
CPPFLAGS += -I$(AZTECPATH)/include
LIBPATH  += $(AZTECPATH)/lib
LDLIBS   += -laztec

WITH_SYSF77 = 1
endif


###############################################################################
ifeq ($(WITH_SPARSKIT),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/SPARSKIT/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/SPARSKIT/include
endif

SPARSKITPATH = /home/semar/pelican/LA_Packages/SPARSKIT2
CPPFLAGS += -I$(AZTECPATH)/lib
LIBPATH  += $(SPARSKITPATH)/lib/Linux
LDLIBS   += -lskit

WITH_SYSF77 = 1
endif


###############################################################################
 ifeq ($(WITH_UMFPACK),1)

     ifeq ($(MAKE_PEL),1)
     SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/UMFPACK/src/*.cc)
     CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/UMFPACK/include
     endif

     UMFPACKPATH = /home/leo/UMFPACKv4.4/UMFPACK/Lib
     CPPFLAGS += -I/home/leo/UMFPACKv4.4/UMFPACK/Include
     LIBPATH  += $(UMFPACKPATH)
     LDLIBS   += -lumfpack  -lamd
     WITH_BLAS = 1
     endif



###############################################################################
ifeq ($(WITH_MPI),1)

# Mandatory definition MPIRUN
# this variable should contain the complete path to the mpirun command
# (without using any other variable)
MPIRUN = /home/semar/pelican/ExternalPackages/etch/openmpi/bin/mpirun

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/MPI/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/MPI/include -DMPIRUN=\"$(MPIRUN)\"
endif

CPPFLAGS += -I$(MPIPATH)/include
MPIPATH   = $(EXTRA_LIBS_DIR)/openmpi
LIBPATH  += $(MPIPATH)/lib
LDLIBS   += -lmpi -lmpi_cxx
endif


###############################################################################
ifeq ($(WITH_SILO),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/Silo/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/Silo/include
endif

SILOPATH  = /home/semar/pelican/VISU_Packages/silo020506
CPPFLAGS += -I$(SILOPATH)/include
LIBPATH   += $(SILOPATH)/lib
LDLIBS   += -lsilo
endif


###############################################################################
ifeq ($(WITH_OPENGL),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/OpenGL/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/OpenGL/include
endif

OPENGLPATH = /usr/
CPPFLAGS += -I$(OPENGLPATH)/include
LIBPATH   += $(OPENGLPATH)/lib
LDLIBS   += -lglut -lGL -lGLU

WITH_X11=1
endif


###############################################################################
ifeq ($(WITH_X11),1)

X11PATH = /usr/
CPPFLAGS += -I$(X11PATH)/include
LIBPATH  += $(X11PATH)/lib
LDLIBS   += -lnsl -lXmu -lXt -lX11 
endif


###############################################################################
ifeq ($(WITH_SYSF77),1)
#LDLIBS   += -lg2c -lm
LDLIBS   += -lm
endif


###############################################################################
ifeq ($(WITH_EXTENDED_MATH),1)
CPPFLAGS += -DEXTENDED_MATH
endif


###############################################################################
ifeq ($(WITH_BLAS),1)
BLASPATH  = /global/software/lib64/intel/GotoBLAS
LIBPATH  += $(BLASPATH)
LDLIBS   += -lblas  -Wl,-rpath,/global/software/lib64/intel/GotoBLAS
endif


###############################################################################
ifeq ($(WITH_ZLIB),1)
ZLIBPATH = /usr
LIBPATH  += $(ZLIBPATH)/include
CPPFLAGS += -I$(ZLIBPATH)/include -DZLIB
LDLIBS   += -lz
endif

###############################################################################
ifeq ($(WITH_METIS),1)

METISPATH = $(EXTRA_LIBS_DIR)/metis
LDLIBS   += $(METISPATH)/libmetis.a

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/METIS_4.0.1/src/*.cc)
CPPFLAGS += -I$(METISPATH)/Lib
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/METIS_4.0.1/include
endif

endif

###############################################################################
ifeq ($(WITH_JNI),1)
#JAVAPATH=/export/opt/J2SDK/current/include
JAVAPATH=/usr/lib/jvm/java-1.5.0-sun/include
JAVAMACHINE=linux
CPPFLAGS += -I$(JAVAPATH) -I$(JAVAPATH)/$(JAVAMACHINE)
endif
#-------------------------------------------------------------------------------
# EXTRA COMPONENT PART : END
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ARCHITECTURE PART : START
#-------------------------------------------------------------------------------
# Inlining /home/leo/Pelican/pelicans/etc/Linux-icc.mak
# * Environment variables should be set prior to invoking the icc compiler.
# * In the intel documentation, it is said to source the file iccvars.sh
#       located in the root installation path (eg /opt/intel/cc/10.0.xxx/bin)
# * Sourcing this file modifies the LD_LIBRARY_PATH environment variable by
#       adding a path for reaching Intel shared libraries (which are not in 
#       the default library path).
# * But this modification of the LD_LIBRARY_PATH environment variable is not 
#       seen by mpirun, which causes an error for parallel execution.
# * Hence, we add the path for intel shared libraries to the LIBPATH variable
#
# * For the default installation:
#       LIBPATH += /opt/intel/cc/10.0.xxx/lib
# * For the installation at IRSN:
LIBPATH += /export/opt/INTEL/cc/10.0.026/lib

# Common Makefile for Linux plateform runnning icc
flagsdbg  = -g
flagsopt0 = -O3
flagsopt1 = -O3
flagsopt2 = -O3

# setting the valid profile flags when needed
ifeq ($(WITH_PROFILE),1)
OPT += -pg
endif
# setting the valid test coverage flags when needed
ifeq ($(WITH_COVERAGE),1)
OPT += -fprofile-arcs -ftest-coverage
endif

CC       = icc
CFLAGS  += $(OPT) -axSSE4.2,SSE4.1 -xSSSE3 -ip -fPIC  

CXX      = icc
CPP      = icc
CXXFLAGS += $(OPT) -axSSE4.2,SSE4.1 -xSSSE3 -ip -fPIC 

FC       = ifort
FFLAGS  += $(OPT) -axSSE4.2,SSE4.1 -xSSSE3 -ip -fPIC

LDFLAGS += $(OPT) -fPIC
LDFLAGS += -lm

LD.so       = icc -shared -o $@ 
LDFLAGS.so  = $(OPT)
DYNAMIC_LIB_EXT = .so
LDLIBSSO = 

MKDEP.c  = gcc -M $(CPPFLAGS)
MKDEP.cc = gcc -M $(CPPFLAGS)

#Generating LDFLAGS from LIBPATH
# a_path => -La_path -Xlinker -rpath -Xlinker a_path
_lpath := $(foreach path,$(LIBPATH), -L$(path) -Xlinker -rpath -Xlinker $(path))
LDFLAGS += $(_lpath)
LDFLAGS.so += $(_lpath)

#Generating CPPFLAGS from INC
# a_path => -Ia_path
CPPFLAGS += $(foreach path,$(INC), -I$(path))
FFLAGS += $(foreach path,$(INC), -I$(path))

#-------------------------------------------------------------------------------
# ARCHITECTURE PART : END
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# GENERIC PART : START
#-------------------------------------------------------------------------------
# Inlining /home/leo/Pelican/pelicans/etc/generic_targets.mak

#standard targets:
# all   : builds the target (TARGET variable) in BINDIR
# obj   : builds the objects only
# clean : clean the dependencies, the objects files and the target
# cleandepend : clean the dependencies
# cleanobj    : clean the objects files
# cleantarget : clean the target file

.PHONY : default_target all dirs
.PHONY : clean cleandepend cleanobj cleantarget
.SUFFIXES :
.SUFFIXES : .c .cc .d .f .F
default_target: all

###############################################################################
#default link mode (cc,c,f,...)
###############################################################################
LINK_MODE=cc

###############################################################################
#declare here the final target and the start sources (which may be in
#different dirs).
###############################################################################
TARGET = exe

MKDIR    = mkdir -p
RMDIR    = rm -fr
DEPFILTER = perl -e 'while(<>) {/^\s*\/usr\// && next;s{^($(*F)).o\s*:}{$(@D)/\1.o $(@D)/\1.d :}; print;}print "\n";'

###############################################################################
#object files definition
NOTDIRSRC := $(notdir $(SRC))
VPATH:=$(subst ' ',':',$(sort $(dir $(SRC))))
OBJ := $(patsubst %.c,   %.o, $(filter %.c,   $(NOTDIRSRC)))\
       $(patsubst %.cc,  %.o, $(filter %.cc,  $(NOTDIRSRC)))\
       $(patsubst %.cpp, %.o, $(filter %.cpp, $(NOTDIRSRC)))\
       $(patsubst %.F,  %.o, $(filter %.F,  $(NOTDIRSRC)))\
       $(patsubst %.f,  %.o, $(filter %.f,  $(NOTDIRSRC)))

OBJ_X = $(addprefix $(BINDIR), $(OBJ))

#target file definition
TARGET_X += $(addprefix $(BINDIR)/, $(TARGET))

#(sub)directories creation in BINDIR
ifdef OBJ_X
DIRS := $(shell $(MKDIR) $(sort $(dir $(OBJ_X))))
endif

#dependencies inclusion
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),cleandepend)
ifneq ($(MAKECMDGOALS),cleanobj)
ifneq ($(MAKECMDGOALS),cleantarget)
ifdef OBJ_X
include $(OBJ_X:.o=.d)
endif
endif
endif
endif
endif

#standard targets
all:$(TARGET_X)
obj: $(OBJ_X)
cleandepend: ; $(RM) $(OBJ_X:.o=.d)
cleanobj:    ; $(RM) $(OBJ_X)
cleantarget: ; $(RM) $(TARGET_X)
clean       : cleandepend cleanobj cleantarget

#the target cmd depends on the suffix (.a, .so or other (executable))
$(TARGET_X):$(OBJ_X)
ifeq ($(suffix $(TARGET_X)),.a)
	$(AR) $(ARFLAGS) $@ $?
else
	$(RM) $(TARGET_X)
ifeq ($(suffix $(TARGET_X)),$(DYNAMIC_LIB_EXT))
	 cd $(BINDIR) ; $(LD.so) $(LDFLAGS.so) $(OBJ) $(PRECOMP_OBJ) $(LDLIBSSO) $(LDLIBS)
else 
	$(LINK.$(LINK_MODE)) -o $@ $(LDFLAGS) $(OBJ_X) $(PRECOMP_OBJ) $(LDLIBS)
endif
endif
###############################################################################
# Production rules ... Don't look further !
###############################################################################

$(BINDIR)%.d:%.F   ; touch $@
$(BINDIR)%.d:%.f   ; touch $@
$(BINDIR)%.d:%.c 
	@echo "Building dependencies for $<"
	@$(MKDEP.c) $< | $(DEPFILTER) > $@
$(BINDIR)%.d:%.cc 
	@echo "Building dependencies for $<"
	@$(MKDEP.c) $< | $(DEPFILTER) > $@
$(BINDIR)%.d:%.cpp 
	@echo "Building dependencies for $<"
	@$(MKDEP.c) $< | $(DEPFILTER) > $@

$(BINDIR)%.o:%.f  ;	$(COMPILE.f)   $< -o $@
$(BINDIR)%.o:%.F  ;	$(COMPILE.f)   $< -o $@
$(BINDIR)%.o:%.c  ;	$(COMPILE.c)   $< -o $@
$(BINDIR)%.o:%.cc ;	$(COMPILE.cc)  $< -o $@
$(BINDIR)%.o:%.cpp ;	$(COMPILE.cc)  $< -o $@
#-------------------------------------------------------------------------------
# GENERIC PART : END
#-------------------------------------------------------------------------------
