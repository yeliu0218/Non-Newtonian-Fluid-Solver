# Common Makefile for Windows plateform runnning gcc
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
CFLAGS  += $(OPT)

CXX      = g++
CPP      = gcc
CXXFLAGS += $(OPT) -DYY_NO_UNISTD_H -Wall -Wno-long-long -Wno-ctor-dtor-privacy -pedantic -W -Wcast-qual -Wwrite-strings -Wconversion -Wshadow -Wno-unused-parameter

FC       = g77
FFLAGS  += $(OPT)

LDFLAGS += $(OPT)
LDFLAGS += -lm -Wl,--enable-auto-import

LD.so       = g++ -shared -o $@ -Wl,--out-implib,$@.a -Wl,--export-all-symbols -Wl,--enable-auto-import -Wl,--add-needed
LDFLAGS.so  = $(OPT)
DYNAMIC_LIB_EXT = .dll
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

