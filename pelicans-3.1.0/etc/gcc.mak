# Common Makefile for plateform runnning gcc
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

CC       = gcc
CFLAGS  += $(OPT) -fPIC

CXX      = g++
CPP      = $(CXX)
CXXFLAGS+= $(OPT) -fPIC -Wall -Wno-ctor-dtor-privacy -pedantic -W -Wcast-qual -Wwrite-strings -Wconversion -Winline -Wshadow -Wno-unused-parameter

FC       = g77
FFLAGS  += $(OPT) -fPIC

LDFLAGS += $(OPT) -fPIC
LDLIBS  += -lnsl -lg2c -lm

LD.so       = $(CXX) -shared -o $@ 
LDFLAGS.so  = $(OPT)
DYNAMIC_LIB_EXT = .so
LDLIBSSO = 

MKDEP.c  = $(CC) -M $(CPPFLAGS)
MKDEP.cc = $(CXX) -M $(CPPFLAGS)

#Generating LDFLAGS from LIBPATH
# a_path => -La_path -Xlinker -rpath -Xlinker a_path
_lpath := $(foreach path,$(LIBPATH), -L$(path) -Xlinker -rpath -Xlinker $(path))
LDFLAGS += $(_lpath)
LDFLAGS.so += $(_lpath)

#Generating CPPFLAGS from INC
# a_path => -Ia_path
CPPFLAGS += $(foreach path,$(INC), -I$(path))
FFLAGS += $(foreach path,$(INC), -I$(path))

