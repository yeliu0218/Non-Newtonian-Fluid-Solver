# Common Makefile for AIX plateform runnning xlC
flagsdbg  = -g
flagsopt0 = -O3
flagsopt1 = -O3
flagsopt2 = -O3
#
# # setting the valid profile flags when needed
ifeq ($(WITH_PROFILE),1)
OPT += -pg
endif
# setting the valid test coverage flags when needed
ifeq ($(WITH_COVERAGE),1)
$(error I dont know how to set "coverage" compilation flags for this architecture.)
endif

#
CC       = xlC
CFLAGS  += $(OPT) 

CXX      = xlC
CPP      = xlC
CXXFLAGS += -qrtti=all $(OPT) 
#
FC       = f77
FFLAGS  += $(OPT) -qextname
#
LDFLAGS += $(OPT) -brtl
LDLIBS  += -lm -lxlf -lxlf90
#
LD.so       = xlC -G -o $@ 
LDFLAGS.so  = $(OPT)
DYNAMIC_LIB_EXT = .so
LDLIBSSO = 
#
MKDEP.c  = gcc -M $(CPPFLAGS)
MKDEP.cc = gcc -M $(CPPFLAGS)
#
#Generating LDFLAGS from LIBPATH
# a_path => -La_path -Xlinker -R -Xlinker a_path
_lpath := $(foreach path,$(LIBPATH), -L$(path) )
LDFLAGS += $(_lpath)
LDFLAGS.so += $(_lpath)
#
#Generating CPPFLAGS from INC
# a_path => -Ia_path
CPPFLAGS += $(foreach path,$(INC), -I$(path))
FFLAGS += $(foreach path,$(INC), -I$(path))
#
