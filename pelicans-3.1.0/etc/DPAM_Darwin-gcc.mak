# Common Makefile for Linux plateform runnning gcc
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

CC       = /opt/local/bin/gcc-mp-4.3
CFLAGS  += $(OPT) -fPIC

CXX      = /opt/local/bin/g++-mp-4.3
CPP      = /opt/local/bin/gcc-mp-4.3
CXXFLAGS += $(OPT) -fPIC -Wall -Wno-long-long -Wno-ctor-dtor-privacy -pedantic -W -Wcast-qual -Wwrite-strings -Wconversion -Wshadow -Wno-unused-parameter

FC       = /opt/local/bin/gfortran-mp-4.3
FFLAGS  += $(OPT) -fPIC

LDFLAGS += $(OPT) -fPIC
LDFLAGS += -lm

LD.so       = /opt/local/bin/g++-mp-4.3 -dynamiclib -o $@ 
LDFLAGS.so  = $(OPT)
DYNAMIC_LIB_EXT = .dylib
LDLIBSSO = 

MKDEP.c  = /opt/local/bin/gcc-mp-4.3 -M $(CPPFLAGS)
MKDEP.cc = /opt/local/bin/gcc-mp-4.3 -M $(CPPFLAGS)

#Generating LDFLAGS from LIBPATH
# a_path => -La_path -Xlinker -rpath -Xlinker a_path
_lpath := $(foreach path,$(LIBPATH), -L$(path))
LDFLAGS += $(_lpath)

#Generating CPPFLAGS from INC
# a_path => -Ia_path
CPPFLAGS += $(foreach path,$(INC), -I$(path))
FFLAGS += $(foreach path,$(INC), -I$(path))

