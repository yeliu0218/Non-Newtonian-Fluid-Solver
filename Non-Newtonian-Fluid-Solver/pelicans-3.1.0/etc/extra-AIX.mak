# Common Makefile for AIX plateform runnning native xlC
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
WITH_UMFPACK  = 0
WITH_METIS    = 0

endif

WITH_JNI  = 0

###############################################################################
ifeq ($(WITH_PETSC),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/PETSc/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/PETSc/include
endif

PETSCARCH = 
PETSCPATH = 
CPPFLAGS += -I$(PETSCPATH)/include -I$(PETSCPATH)/bmake/$(PETSCARCH)
LIBPATH   += $(PETSCPATH)/lib/libO/$(PETSCARCH)
LDLIBS   += -lpetscgsolver -lpetscgrid -lpetscmesh -lpetscts -lpetscsnes -lpetscsles -lpetscdm -lpetscmat -lpetscvec -lpetsc -lmpe

WITH_X11 = 1
WITH_BLAS = 1
endif


###############################################################################
ifeq ($(WITH_AZTEC),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/Aztec/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/Aztec/include
endif

AZTECPATH = /home/semar/pelican/LA_Packages/Aztec
CPPFLAGS += -I$(AZTECPATH)/lib
LIBPATH   += $(AZTECPATH)/lib/AIX
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
CPPFLAGS += -I$(SILOPATH)/lib
LIBPATH   += $(SPARSKITPATH)/lib/AIX
LDLIBS   += -lskit

WITH_SYSF77 = 1
endif


###############################################################################
ifeq ($(WITH_UMFPACK),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/UMFPACK/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/UMFPACK/include
endif

UMFARCH =
UMFPACKPATH = 
CPPFLAGS += -I$(UMFPACKPATH)/UMFPACK/Include
LIBPATH  += $(UMFPACKPATH)/UMFPACK/Lib/$(UMFARCH) $(UMFPACKPATH)/AMD/Lib/$(UMFARCH)
LDLIBS   += -lumfpack -lamd

WITH_BLAS = 1
endif


###############################################################################
ifeq ($(WITH_MPI),1)

# Mandatory definition MPIRUN
# this variable should contain the complete path to the mpirun command
# (without using any other variable)
MPIRUN = 

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/MPI/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/MPI/include -DMPIRUN=\"$(MPIRUN)\"
endif

CPPFLAGS += -I$(MPIPATH)/include
MPIPATH   =
LIBPATH  += $(MPIPATH)/lib
LDLIBS   += -lmpi -lmpi_cxx

endif


###############################################################################
ifeq ($(WITH_SILO),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/Silo/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/Silo/include
endif

SILOPATH = /usr/local
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

OPENGLPATH = /usr/local
CPPFLAGS += -I$(OPENGLPATH)/include
LIBPATH   += $(OPENGLPATH)/lib
LDLIBS   += -lglut -lGL -lGLU

WITH_X11=1
endif


###############################################################################
ifeq ($(WITH_X11),1)
LDLIBS   += -lnsl -lXmu -lXt -lX11 
endif


###############################################################################
ifeq ($(WITH_SYSF77),1)
# LDLIBS   +=
endif


###############################################################################
ifeq ($(WITH_FLEX),1)
LIBPATH += /usr/local/lib
LDLIBS  += -lfl
endif
ifeq ($(WITH_LEX),1)
LDLIBS   += -ll
endif


###############################################################################
ifeq ($(WITH_EXTENDED_MATH),1)
CPPFLAGS += -DEXTENDED_MATH
endif


###############################################################################
ifeq ($(WITH_BLAS),1)
BLASPATH =
LIBPATH  += $(BLASPATH)
LDLIBS   += -lblas -llapack
endif

###############################################################################
ifeq ($(WITH_METIS),1)

ifeq ($(MAKE_PEL),1)
METISPATH = /home/semar/pelican/METIS/AIX/metis-4.0
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/METIS_4.0.1/src/*.cc)
CPPFLAGS += -I$(METISPATH)/Lib
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/METIS_4.0.1/include
LDLIBS   += $(METISPATH)/libmetis.a
endif

endif

###############################################################################
ifeq ($(WITH_JNI),1)
JAVAPATH=/usr/java/include
CPPFLAGS += -I$(JAVAPATH)
endif
