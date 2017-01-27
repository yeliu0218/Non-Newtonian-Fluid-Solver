# Common Makefile for Linux plateform runnning gcc
# Settings for Extra Components

ifeq ($(LINK_PEL),1)

WITH_MPI      = 0
##### PETSc external API
# if PETSc is enabled, uncomment the appropriate definition of PETSc_VERSION
WITH_PETSc    = 0
# PETSc_VERSION=2.3.0
# PETSc_VERSION=2.3.1
PETSc_VERSION=2.3.3

WITH_AZTEC    = 0
WITH_SPARSKIT = 0
WITH_SILO     = 0
WITH_OPENGL   = 1
WITH_UMFPACK  = 0
WITH_METIS    = 0

endif

WITH_ZLIB  = 1
WITH_EXTENDED_MATH = 1
WITH_JNI = 0

###############################################################################
ifeq ($(WITH_PETSc),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/PETSc_$(PETSc_VERSION)/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/PETSc_$(PETSc_VERSION)/include
endif

PETSc_ARCH=Darwin
PETSc_DIR=/Users/piar/soft/PETSc/petsc-2.3.3-p15
CPPFLAGS += -I$(PETSc_DIR)/include -I$(PETSc_DIR)/bmake/$(PETSc_ARCH)
LIBPATH  += $(PETSc_DIR)/lib/$(PETSc_ARCH)
# LDLIBS   += -L$(PETSc_DIR)/lib/$(PETSc_ARCH) -lpetscsnes -lpetscksp -lpetscmat -lpetscvec -lpetsc -lpetscdm
LDLIBS   += $(PETSc_DIR)/lib/$(PETSc_ARCH)/libpetscsnes.a $(PETSc_DIR)/lib/$(PETSc_ARCH)/libpetscksp.a $(PETSc_DIR)/lib/$(PETSc_ARCH)/libpetscmat.a $(PETSc_DIR)/lib/$(PETSc_ARCH)/libpetscvec.a $(PETSc_DIR)/lib/$(PETSc_ARCH)/libpetsc.a $(PETSc_DIR)/lib/$(PETSc_ARCH)/libpetscdm.a

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
CPPFLAGS += -I$(SILOPATH)/lib
LIBPATH   += $(AZTECPATH)/lib/Linux
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
LIBPATH   += $(SPARSKITPATH)/lib/Linux
LDLIBS   += -lskit

WITH_SYSF77 = 1
endif


###############################################################################
ifeq ($(WITH_UMFPACK),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/UMFPACK/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/UMFPACK/include
endif

UMFARCH = Linux
UMFPACKPATH = /home/semar/pelican/LA_Packages/UMFPACKv4.3
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
MPIRUN = /Users/piar/soft/openmpi-1.2.8/bin/mpirun

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/MPI/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/MPI/include -DMPIRUN=\"$(MPIRUN)\"
endif

# with openmpi, the options may be obtained with the commands
#    openmpicxx -showme:compile
#    openmpicxx -showme:link
CPPFLAGS += -D_REENTRANT -I/Users/piar/soft/openmpi-1.2.8/include
LDFLAGS  += -Wl,-u,_munmap -Wl,-multiply_defined,suppress 
# ????? où est utilisé LIBPATH ????
# LIBPATH  += /opt/local/lib
LDLIBS   += -L/Users/piar/soft/openmpi-1.2.8/lib -lmpi_cxx -lmpi -lopen-rte -lopen-pal 

endif


###############################################################################
ifeq ($(WITH_SILO),1)

ifeq ($(MAKE_PEL),1)
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/Silo/src/*.cc)
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/Silo/include
endif

SILOPATH = /home/semar/pelican/share/Linux
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

LDLIBS += -framework OpenGL -framework GLUT

WITH_X11=0
endif


###############################################################################
ifeq ($(WITH_X11),1)

X11PATH = /usr/X11R6
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
BLASPATH = /usr/lib/atlas/3dnow
BLASPATH2 = /usr/lib/3dnow/atlas
LIBPATH  += $(BLASPATH)
LIBPATH  += $(BLASPATH2)
LDLIBS   += -lblas -llapack
endif


###############################################################################
ifeq ($(WITH_METIS),1)

ifeq ($(MAKE_PEL),1)
METISPATH = /Users/piar/soft/metis-4.0
SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/METIS_4.0.1/src/*.cc)
CPPFLAGS += -I$(METISPATH)/Lib
CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/METIS_4.0.1/include
LDLIBS   += $(METISPATH)/libmetis.a
endif
endif

###############################################################################
ifeq ($(WITH_JNI),1)
LDLIBS += -framework JavaVM
JAVAPATH=/System/Library/Frameworks/JavaVM.framework/Headers
CPPFLAGS += -I$(JAVAPATH)
endif

