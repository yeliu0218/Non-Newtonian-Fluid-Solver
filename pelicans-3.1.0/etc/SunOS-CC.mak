# Common Makefile for SunOS plateform running native CC
flagsdbg  = -g
flagsopt0 = -fast -xtarget=native 
flagsopt1 = -fast -xtarget=native 
flagsopt2 = -fast -xtarget=native 

# setting the valid profile flags when needed
ifeq ($(WITH_PROFILE),1)
OPT += -xpg
endif
# setting the valid test coverage flags when needed
ifeq ($(WITH_COVERAGE),1)
OPT += -xprofile=tcov
endif

CC       = cc
CFLAGS  += $(OPT) -KPIC

CXX      = CC
CPP      = cpp
CXXFLAGS += $(OPT) -KPIC

FC       = f77
FFLAGS  += $(OPT) -KPIC

LDFLAGS += $(OPT) -KPIC
LDLIBS  += -lm

LD.so       = CC -G -o $@ 
LDFLAGS.so  = $(OPT)
DYNAMIC_LIB_EXT = .so
LDLIBSSO =  -lCstd -lCrun

MKDEP.c  = gcc -M $(CPPFLAGS)
MKDEP.cc = gcc -M $(CPPFLAGS)

#Generating LDFLAGS from LIBPATH
# a_path => -La_path -Ra_path
_lpath := $(foreach path,$(LIBPATH), -L$(path) -R$(path))
LDFLAGS += $(_lpath)
LDFLAGS.so += $(_lpath)

#Generating CPPFLAGS from INC
# a_path => -Ia_path
CPPFLAGS += $(foreach path,$(INC), -I$(path))
FFLAGS += $(foreach path,$(INC), -I$(path))

