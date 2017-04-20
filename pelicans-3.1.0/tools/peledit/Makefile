#
#  Recompilation de l'interface JAVA/C++ pour pelguis.
#

.SUFFIXES: .c .h

CCC := gcc
OPT_TEST := opt0

# PELICANS
PELEDITHOME = $(shell pwd)
PELICANSHOME = $(PELEDITHOME)/../..
pel_arch   = cd $(PELICANSHOME) && perl -I $(PELICANSHOME)/tools/pel $(PELICANSHOME)/tools/pel/arch.pl
pel_depend = perl -I $(PELICANSHOME)/tools/pel $(PELICANSHOME)/tools/pel/depend.pl

SYS := $(shell $(pel_arch) $(CCC))
ifndef SYS
   $(shell $(pel_arch) -v $(CCC))
   $(error Unable to find a valid architecture for this computer)
endif

WITH_JNI := $(shell $(pel_arch) -getvariable_extra=WITH_JNI $(CCC))

CPP_SRC = $(PELICANSHOME)/PELbase $(PELICANSHOME)/RefSol $(PELICANSHOME)/Geometry/src/GE_RoundoffExp.cc $(PELICANSHOME)/Geometry/src/GE_ColorExp.cc $(PELICANSHOME)/Geometry/src/GE_Color.cc $(PELICANSHOME)/Geometry/src/GE_PerturbatedMeshingExp.cc
CPP_LIB = $(PELICANSHOME)/lib/$(SYS)/$(OPT_TEST)
DYNAMIC_LIB_EXT = $(shell $(pel_arch) -getvariable DYNAMIC_LIB_EXT)

# Directory:
jardir = jars
libdir = $(jardir)/lib
tmpdir = $(libdir)/$(SYS)

#JNI Environment
LIBJNI = $(tmpdir)/PELjni$(DYNAMIC_LIB_EXT)
JNIDIR = jnisrc
JNISRC = $(wildcard $(JNIDIR)/*.cc)
JNIOBJECTS = $(patsubst %.cc, %.o, $(JNISRC))
JNIHEADERS = $(wildcard $(JNIDIR)/*.h)
RENAME = java -cp $(jardir)/pelguis.jar RenameJni

all:
ifeq ($(WITH_JNI),1)
	@ $(MAKE) jni
endif

info:
	@echo "Making peledit JNI:"
	@echo "   ARCH = $(SYS)"
	@echo "   PELICANSHOME = $(PELICANSHOME)"
	@echo "   LEVEL = $(OPT_TEST)"

clean:
	@$(RM) $(libdir)
	@$(RM) $(jardir)/lib*

jni : info $(LIBJNI)
	@chmod -R  +r $(jardir)/*
	@chmod -R g+r $(jardir)/*
	@chmod -R o+r $(jardir)/*

doc :

$(LIBJNI) : $(tmpdir)/Makefile
ifneq ($(WITH_JNI),1)
	$(error "JNI path not defined in configuration file (cf. etc directory)")
endif
	@cd $(tmpdir) ; $(MAKE) TARGET=$(@F)
	@$(RM) $(jardir)/$(@F)
	@$(MV) $(tmpdir)/$(@F) $(jardir)/$(@F)
	@$(RENAME) $(jardir)/$(@F)
	@$(RM) $(tmpdir)

$(tmpdir):
	@$(MKDIR) $@

$(tmpdir)/Makefile : $(tmpdir)
ifeq ($(OPT_TEST),optpg)
	@$(pel_depend) -compiler $(CCC) opt0 $(tmpdir) $(JNIDIR) $(CPP_SRC) $(patsubst %, -precomp %, $(CPP_LIB))
else
	@$(pel_depend) -compiler $(CCC) $(OPT_TEST) $(tmpdir) $(JNIDIR) $(CPP_SRC) $(patsubst %, -precomp %, $(CPP_LIB))
endif

help :
	@echo "Usage : $(MAKE) <targets ...> with targets names in :"
	@echo "   all         : builds jnis"
	@echo "   clean       : clean temporary files"

MKDIR = mkdir -p
RM = rm -rf
MV = mv
