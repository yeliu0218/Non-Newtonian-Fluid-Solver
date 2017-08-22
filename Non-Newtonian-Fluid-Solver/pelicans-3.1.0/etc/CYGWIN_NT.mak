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
$(error I dont know how to set "coverage" compilation flags for this architecture.)
endif

CC       = gcc
CFLAGS  += $(OPT)

CXX      = g++
CPP      = gcc
CXXFLAGS += $(OPT) -Wall -Wno-ctor-dtor-privacy -pedantic -W -Wcast-qual -Wwrite-strings -Wconversion -Winline -Wshadow -Wno-unused-parameter

FC       = g77
FFLAGS  += $(OPT)

LDFLAGS += $(OPT) -Wl,--dll-search-prefix=lib
LDLIBS  += -lg2c -lm

LD.so       = g++ -shared -Wl,--out-implib=$@.a -o $@
#$(@D)/$(subst lib,cyg,$(@F))
LDFLAGS.so  = $(OPT) -Wl,--export-all-symbols,--enable-auto-import,--whole-archive,--dll-search-prefix=lib
DYNAMIC_LIB_EXT = .dll
LDLIBSSO = -Wl,--no-whole-archive -lfl

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

