# Common Makefile for Linux plateform runnning gcc
flagsdbg  = -g
flagsopt0 = -O
flagsopt1 = -O
flagsopt2 = -O

# setting the valid profile flags when needed
ifeq ($(WITH_PROFILE),1)
OPT += -pg
endif
# setting the valid test coverage flags when needed
ifeq ($(WITH_COVERAGE),1)
$(error I dont know how to set "coverage" compilation flags for this architecture.)
endif

CC       = cc
CXX      = cxx
CPP      = cxx
CXXFLAGS += $(OPT) -fPIC -Wall -Wno-ctor-dtor-privacy -pedantic -W -Wcast-qual -Wwrite-strings -Wconversion -Winline -Wshadow -Wno-unused-parameter
CPPFLAGS += -D__USE_STD_IOSTREAM


FC       = f77
FFLAGS  += $(OPT) -fPIC

LD.cc    = cxx
LDFLAGS += $(OPT) -fPIC
LDLIBS  += -lUfor -lFutil -lfor -lm
LIBPATH += "/usr/shlib"

LD.so       = cxx -shared -o \$@
LDFLAGS.so  = $(OPT)
DYNAMIC_LIB_EXT = .so
LDLIBSSO = 

MKDEP.c  = cxx -M -noimplicit_include $(CPPFLAGS)
MKDEP.cc = cxx -M -noimplicit_include $(CPPFLAGS)

#Generating LDFLAGS from LIBPATH
# a_path => -La_path -Xlinker -rpath -Xlinker a_path
_lpath := $(foreach path,$(LIBPATH), -L$(path) -rpath $(path))
LDFLAGS += $(_lpath)
LDFLAGS.so += $(_lpath)

#Generating CPPFLAGS from INC
# a_path => -Ia_path
CPPFLAGS += $(foreach path,$(INC), -I$(path))
FFLAGS += $(foreach path,$(INC), -I$(path))

