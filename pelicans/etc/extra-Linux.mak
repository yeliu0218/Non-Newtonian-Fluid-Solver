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
