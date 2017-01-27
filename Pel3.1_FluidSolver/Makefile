##################################################################
# Standard GNU Makefile for PELICANS Applications                #
##################################################################

# Set here the PELICANS environment configuration file
CONFIG = $(PELICANSHOME)/bin/init

# Set here the directories where the pel scripts will recursively
# find the source files of the application
SRC = .

# default compilation system (used to find the architecture)
CCC := icc

# the name of the application (used in the documentation)
APPLINAME := pel-application

# Test if the current architecture is valid for PELICANS
SYS    := $(shell . $(CONFIG).sh ; pel arch $(CCC))
ifndef SYS
$(shell . $(CONFIG).sh ; pel arch  -v $(CCC))
$(error Unable to find a valid architecture for this computer)
endif

# directory where objects files and binaries are stored
libdir := lib/$(SYS)

# targets that are not the name of a file. There are just names
# for some commands to be executed on an explicit request
.PHONY : exe all clean exe% FORCE

# Tells make not to remove intermediate target such as Makefile(s) in $(libdir)
.SECONDARY :

# the default target
all : help

clean:
	$(BEGIN_TARGET) clean
	$(RM) lib 	
	$(RM) doc
	$(RM) tests
	$(END_TARGET) clean

cleanlib:
	$(BEGIN_TARGET) cleanlib
	$(RM) lib
	$(END_TARGET) cleanlib

cleandoc:
	$(BEGIN_TARGET) cleandoc
	$(RM) doc
	$(END_TARGET) cleandoc


##################################################################
# Executables                                                    #
# the user can make predefined targets :                         #
#       make exe2 => application opt2, pelicans opt1             #
#       make exe0 => application opt0, pelicans opt0             #
# the user may also build non predefined targets :               #
#       make exe2 PELLIB=pel0 => application opt2, pelicans opt0 #
#       make exeg PELLIB=pel2 => application dbg, pelicans opt2  #
# so more generally,                                             #
#       make exeOOO PELLIB=PPP => application OOO, pelicans PPP  #
# some very special targets are :                                #
#       make exegg => application dbg, pelicans dbg              #
#       make exepg => application optpg, pelicans pel0           #
#       make exepgpg => application optpg, pelicans optpg        #
#            optpg option is used to profile an application      #
##################################################################

# the application is linked with the PELLIB version of the PELICANS librarie
PELLIB = pel0

# default target generates 3 executables
exe: exe0 exe2 exeg

# overrides PELLIB for some known targets
exe2    : PELLIB=pel1
exe22   : PELLIB=pel2
exeg    : PELLIB=pel1
exegg   : PELLIB=pelg
exepgpg : PELLIB=$(PELICANSHOME)/$(libdir)/optpg

$(libdir)/opt2%/Makefile: 
	$(BEGIN_TARGET) $@
	$(MKDIR) $(@D)
	. $(CONFIG).sh ; pel depend -l pel$* -compiler $(CCC) dbg $(@D) $(SRC)

$(libdir)/optg/Makefile: 
	$(BEGIN_TARGET) $@
	$(MKDIR) $(@D)
	. $(CONFIG).sh ; pel depend -l pel1 -compiler $(CCC) dbg $(@D) $(SRC)

$(libdir)/optg%/Makefile: 
	$(BEGIN_TARGET) $@
	$(MKDIR) $(@D)
	. $(CONFIG).sh ; pel depend -l pel$* -compiler $(CCC) dbg $(@D) $(SRC)

$(libdir)/optpg%/Makefile: 
	$(BEGIN_TARGET) $@
	$(MKDIR) $(@D)
	. $(CONFIG).sh ; pel depend -l $(PELLIB) -compiler $(CCC) optpg $(@D) $(SRC)

$(libdir)/%/Makefile:
	$(BEGIN_TARGET) $@
	$(MKDIR) $(@D)
	. $(CONFIG).sh ; pel depend -l $(PELLIB) -compiler $(CCC) $* $(@D) $(SRC)
	$(END_TARGET) $@

exe% : $(libdir)/exe%
	$(END_TARGET) $@

$(libdir)/exe% : $(libdir)/opt%/Makefile FORCE
	$(BEGIN_TARGET) $@
	. $(CONFIG).sh ; pel build -exe $(libdir)/opt$*
	if [ ! -f $@ ] ; then cd $(@D) ; ln -s opt$*/exe $(@F) ; fi
	$(END_TARGET) $@

##################################################################
# Test                                                           #
# Builds the test executable and runs the application            #
# non regression tests. Stores results and reports in the tests  #
# subdirectory.                                                  #
# By default, exe0 is used to run the tests but this may be      #
# changed like this : make tests exe_test=exe1                   #
##################################################################

exe_test = exe0
libtest  = tests/$(SYS)

test: $(libtest) FORCE
	$(BEGIN_TARGET) $@
	cd $(libtest) ; . $(CONFIG).sh ; pel test -build_pattern pattern.pel ../../$(libdir)/exe0 ../../RegressionTests
	cd $(libtest) ; . $(CONFIG).sh ; pel test -verify_pattern pattern.pel ../../$(libdir)/exe0 ../../RegressionTests
	$(END_TARGET) $@

$(libtest):
	$(MKDIR) $@


##################################################################
# Documentation                                                  #
##################################################################

doc: $(exe_test) FORCE
	$(BEGIN_TARGET) $@
	$(RM) doc/Html
	$(MKDIR) doc/Html
	. $(CONFIG).sh ; pel predoc -Wno_unresolved doc/Html/description.txt doc/Html/data_doc.pel $(APPLINAME) .
	cd doc/Html ; . $(CONFIG).sh ; pel run -noverb -noautocheck ../../$(libdir)/exe0 data_doc.pel peldoc.out
	$(END_TARGET) $@

##################################################################
# Help                                                           #
##################################################################

ifdef CALLER
COMMANDE := $(CALLER)
else
COMMANDE := $(MAKE)
endif

help:
	@echo
	@echo "Usage: $(COMMANDE) CCC=<arg> <target>"
	@echo "<target>"
	@echo "    exe0  : build executable : opt0  ; libpel0.so"
	@echo "    exe2  : build executable : opt2  ; libpel1.so"
	@echo "    exe22 : build executable : opt2  ; libpel2.so"
	@echo "    exeg  : build executable : dbg   ; libpel1.so"
	@echo "    exegg : build executable : dbg   ; libpelg.so"
	@echo "    exepg : build executable : optpg ; libpel0.so"
	@echo "    exepgpg : build full profiled executable"
	@echo "    exe   : same as exe0, exe2 and exeg"
	@echo "    test  : regression testing"
	@echo "    clean : remove all generated files (lib, tests, doc)"
	@echo "    cleanlib : remove generated files (lib)"
	@echo "    cleandoc : remove generated files (doc)"
	@echo "<arg> : compiler system (default: gcc)"
	@echo "detected architecture:" $(SYS)
	@echo "*******************************************"
	@echo "Reminder: Create documentation with doxygen"
	@echo "*******************************************"
	@echo

FORCE:

##################################################################
# internal commands                                              #
##################################################################

RM    := rm -fr
MKDIR := mkdir -p
BEGIN_TARGET := @printf "****************************************\n* \
Begin target : %s \n****\n" 
END_TARGET := @printf "****\n* \
target %s completed\n****************************************\n"

