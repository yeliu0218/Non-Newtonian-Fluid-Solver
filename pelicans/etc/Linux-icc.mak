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

