# Common Makefile for Linux plateform running CC
flagsdbg  = -g
flagsopt0 = -fast -xtarget=native
flagsopt1 = -fast -xtarget=native
flagsopt2 = -fast -xtarget=native

CC       = cc
CFLAGS  += $(OPT) -KPIC

CXX      = CC
CPP      = cpp
CXXFLAGS += $(OPT) -KPIC 

FC       = f77
FFLAGS  += $(OPT) -KPIC

LDFLAGS += $(OPT) -KPIC
LDFLAGS += -lm

LD.so       = CC -G -o $@ 
LDFLAGS.so  = $(OPT)
DYNAMIC_LIB_EXT = .so
LDLIBSSO = -lCstd -lCrun 

MKDEP.c  = gcc -M $(CPPFLAGS)
MKDEP.cc = gcc -M $(CPPFLAGS)

#Generating LDFLAGS from LIBPATH
# a_path => -La_path -Xlinker -rpath -Xlinker a_path
_lpath := $(foreach path,$(LIBPATH), -L$(path) -R$(path))
LDFLAGS += $(_lpath)
LDFLAGS.so += $(_lpath)

#Generating CPPFLAGS from INC
# a_path => -Ia_path
CPPFLAGS += $(foreach path,$(INC), -I$(path))
FFLAGS += $(foreach path,$(INC), -I$(path))

