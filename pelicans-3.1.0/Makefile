#
# Makefile for PELICANS platform
#
LICENSE := "accept"
CCC := gcc
OPT_TEST := opt0
RM    := rm -fr
MKDIR := mkdir -p
DATE  := date

SYS:=$(shell perl -I tools/pel tools/pel/arch.pl $(CCC))
ifndef SYS
$(shell perl -I tools/pel tools/pel/arch.pl -v $(CCC))
$(error Unable to find a valid architecture for this computer)
endif

CONFIG=bin/init
BEGIN_TARGET=@printf "****************************************\n* \
Begin target : %s \n****\n" 
END_TARGET=@printf "****\n* \
target %s completed\n****************************************\n"

# determination de l'extension des DLL par perl
DYNAMIC_LIB_EXT := $(shell perl -I tools/pel tools/pel/arch.pl $(CCC) -getvariable DYNAMIC_LIB_EXT)

libdir = lib/$(SYS)
abslibdir := $(shell perl -MFile::Spec -e 'print File::Spec->rel2abs("$(libdir)")')
lib0 = $(libdir)/libpel0$(DYNAMIC_LIB_EXT)
lib1 = $(libdir)/libpel1$(DYNAMIC_LIB_EXT)
lib2 = $(libdir)/libpel2$(DYNAMIC_LIB_EXT)
libg = $(libdir)/libpelg$(DYNAMIC_LIB_EXT)


LIB_TEST = $(subst dbg,pelg,$(subst opt,pel,$(OPT_TEST)))
lib_for_test = $(libdir)/lib$(LIB_TEST)$(DYNAMIC_LIB_EXT)
export LIB_TEST

PELICANS_VCPROJECT_2003 = admin/win32/vc2003/pelicans.vcproj
PELICANS_VCPROJECT_2008 = admin/win32/vc2008/pelicans.vcproj
STD_VCPROJECT_2008 = admin/win32/vc2008/examplesofapplication.vcproj

##########################
# MACRO TARGETS          #
##########################

default:help

##################################
# GNU-like MACRO TARGETS         #
##################################

all:libraries

clean:clean_no_doc

#####################################
# MAINTAINER Installing/Updating    #
#####################################

maintainer-update: libraries extra-libraries  tools testing 
	$(END_TARGET) $@

maintainer-install:cleangenerated maintainer-clean generated doc vcproj dist maintainer-update
	$(END_TARGET) $@

maintainer-build_and_smoke_tests: cleangenerated maintainer-clean generated doc vcproj dist tools testing
	$(END_TARGET) maintainer-build_and_smoke_tests

#######################
# cleaner             #
#######################
#DIRTOOLS:=$(shell ls tools/*/Makefile 2>/dev/null)
DIRTOOLS:=$(wildcard tools/*/Makefile)
DIRTOOLS.clean:=$(addsuffix .clean,$(DIRTOOLS))
cleanenv:
	$(RM) $(CONFIG).csh $(CONFIG).sh
cleangenerated:
	cd PELbase/grammar ; $(MAKE) clean
	cd DocApplication/grammar ; $(MAKE) clean
cleanlib:
	$(RM) lib
cleantests:
	$(RM) tests
	$(RM) testing
cleantools:$(DIRTOOLS.clean)

cleaninclude:FORCE
	$(RM) include

cleandoc: cleantools cleaninclude
	$(RM) doc/Html doc/Site/APPLICATION_Documentations doc/Site/pel_CLI doc/man

clean_no_doc:FORCE cleanenv cleanlib cleantests cleandist cleanpeledit cleanvcproj
	$(RM) doc/Html/src

maintainer-clean:clean_no_doc cleandoc
	$(RM) $(STAMP)
	$(END_TARGET) $@

LIBSDIR:=$(wildcard lib/*/*/Makefile)
LIBSDIR.strip:=$(addsuffix .strip,$(LIBSDIR))
LIBTESTSDIR:=$(wildcard tests/lib/*/*/Makefile)
LIBTESTSDIR.strip:=$(addsuffix .strip,$(LIBSDIR))
strip:FORCE  $(LIBSDIR.strip) $(LIBTESTSDIR.strip)
	$(BEGIN_TARGET) $@
	$(RM) tests/SingleApplication tests/doc tests/pCheckTime.bd
	$(RM) $(STAMP)
	$(END_TARGET) $@

%.strip:$(CONFIG).sh
	$(BEGIN_TARGET) $*
	. $(CONFIG).sh ; $(MAKE) -f $* cleanobj ; $(MAKE) -f $* cleandepend
	$(END_TARGET) $*

cleanmakefiles:FORCE
	@find lib -name "Makefile" -exec $(RM) {} \; -print
	@find tests/lib -name "Makefile" -exec $(RM) {} \; -print
	@find tools/peledit/jars -name "Makefile" -exec $(RM) {} \; -print

#######################
# Environment setting #
#######################

environment:$(CONFIG).sh 

$(CONFIG).sh:
	admin/setEnvironment $(CONFIG) $(LICENSE)
	$(END_TARGET) environment 

#######################
# Generated source    #
#######################

generated:FORCE
	cd PELbase/grammar ; $(MAKE) 
	cd DocApplication/grammar ; $(MAKE) 
	$(END_TARGET)  $@

#######################
# Tools builder       #
#######################
DIRTOOLS.all:=$(addsuffix .all,$(DIRTOOLS))
tools:$(DIRTOOLS.all)

%.all:$(CONFIG).sh
	$(BEGIN_TARGET) $*
	. $(CONFIG).sh ; cd $(@D) ; $(MAKE) libdir=$(abslibdir) CCC=$(CCC) LIB_TEST=$(LIB_TEST) OPT_TEST=$(OPT_TEST)
	$(END_TARGET) $*

%.clean:
	$(BEGIN_TARGET) $*
	cd $(@D) ; $(MAKE) clean
	$(END_TARGET) $*

#######################
# Include directory   #
#######################

include:$(CONFIG).sh
	$(BEGIN_TARGET) $@
	. $(CONFIG).sh ; cd admin ; ./setDoc
	$(END_TARGET) $@

$(PELICANS_VCPROJECT_2003) $(PELICANS_VCPROJECT_2008) $(STD_VCPROJECT_2008) vcproj:FORCE $(CONFIG).sh
	$(BEGIN_TARGET) $@
	cp $(PELICANS_VCPROJECT_2003)_template $(PELICANS_VCPROJECT_2003)  
	. $(CONFIG).sh ; pel depend -vcproj $(PELICANS_VCPROJECT_2003) -mSTD -mPELICANS
	cp $(PELICANS_VCPROJECT_2008)_template $(PELICANS_VCPROJECT_2008)
	. $(CONFIG).sh ; pel depend -vcproj $(PELICANS_VCPROJECT_2008) -mPELICANS
	cp  $(STD_VCPROJECT_2008)_template  $(STD_VCPROJECT_2008)
	. $(CONFIG).sh ; pel depend -vcproj $(STD_VCPROJECT_2008) -mSTD
	$(END_TARGET) $@

cleanvcproj:FORCE
	$(RM) $(PELICANS_VCPROJECT_2003) $(PELICANS_VCPROJECT_2008) $(STD_VCPROJECT_2008)

#######################
# Libraries builder   #
#######################

lib0 : $(lib0)
lib1 : $(lib1)
lib2 : $(lib2)
libg : $(libg)
libpg : $(libdir)/libpelpg$(DYNAMIC_LIB_EXT)

libraries : lib1 lib0
	$(END_TARGET) $@

extra-libraries : lib2 libg
	$(END_TARGET) $@

.PRECIOUS:$(libdir)/%/lib$(DYNAMIC_LIB_EXT) $(libdir)/%/Makefile

$(libdir)/libpelg$(DYNAMIC_LIB_EXT): $(libdir)/dbg/Makefile $(CONFIG).sh FORCE
	. $(CONFIG).sh ; pel build -make $(MAKE) -archive ../$(@F) $(<D)

$(libdir)/libpel%$(DYNAMIC_LIB_EXT): $(libdir)/opt%/Makefile $(CONFIG).sh FORCE
	. $(CONFIG).sh ; pel build -make $(MAKE) -archive ../$(@F) $(<D)

$(libdir)/%/Makefile: $(CONFIG).sh
	$(MKDIR) $(@D)
	. $(CONFIG).sh ; pel depend -mPELICANS -compiler $(CCC) $* $(@D)


#######################
# PelEditor builder   #
#######################

peledit: $(CONFIG).sh $(lib_for_test)
	$(BEGIN_TARGET) $@
	. $(CONFIG).sh ; cd tools/peledit ; $(MAKE) jni
	$(END_TARGET) $@

cleanpeledit:
	cd tools/peledit ; $(MAKE) clean


############################
# Examples of applications #
############################

libdirtests=tests/$(libdir)/$(OPT_TEST)

exetests=$(libdirtests)/exe
makefile_lib_exe=$(libdirtests)/Makefile

exe_test:$(lib_for_test) $(exetests)
	$(END_TARGET) $@

$(exetests):$(libdirtests) $(makefile_lib_exe) FORCE
	. $(CONFIG).sh ;  pel build -make $(MAKE) -exe $(libdirtests)

$(libdirtests):
	$(MKDIR) $(libdirtests)

$(makefile_lib_exe):include 
ifeq ($(OPT_TEST),optpg)
	. $(CONFIG).sh ;  pel depend -l $(libdir)/optpg -compiler $(CCC) optpg $(libdirtests) ExamplesOfApplication
else
	. $(CONFIG).sh ;  pel depend -l $(LIB_TEST) -compiler $(CCC) -path $(abslibdir) $(OPT_TEST) $(libdirtests) ExamplesOfApplication
endif


#########################
# Single application    #
#########################

sa_libdirtests=tests/SingleApplication/$(SYS)/$(OPT_TEST)
sa_exetests=$(sa_libdirtests)/exe
sa_makefile_lib_exe=$(sa_libdirtests)/Makefile

sa_exe_test:$(lib_for_test) $(sa_exetests)
	$(END_TARGET) $@

$(sa_exetests):$(sa_libdirtests) $(sa_makefile_lib_exe) FORCE
	. $(CONFIG).sh ;  pel build -make $(MAKE) -exe $(sa_libdirtests)

$(sa_libdirtests):
	$(MKDIR) $(sa_libdirtests)

$(sa_makefile_lib_exe):include 
ifeq ($(OPT_TEST),optpg)
	. $(CONFIG).sh ;  pel depend -l $(libdir)/optpg -compiler $(CCC) optpg $(sa_libdirtests) SingleApplication
else
	. $(CONFIG).sh ;  pel depend -l $(LIB_TEST) -compiler $(CCC) -path $(abslibdir) $(OPT_TEST) $(sa_libdirtests) SingleApplication
endif

#########
# Tests #
#########

libtest=tests/$(SYS)
test_stamp=$(libtest)/testing

checkfull:$(test_stamp) 

testing:checkfull timing pelgenerated_doc
	. $(CONFIG).sh ; admin/verify_building.sh $(lib_for_test) $(exetests) > testing
	$(END_TARGET) $@


$(test_stamp):$(libtest) exe_test sa_exe_test
	. $(CONFIG).sh ; cd $(libtest) ; pel test -exact -build_then_verify_pattern pattern.pel ../../$(exetests) ../../UnitTests ../../ExamplesOfApplication 
	. $(CONFIG).sh ; cd $(libtest) ; pel test -exact -peltest_exe ../../$(exetests) ../../$(sa_exetests) ../../SingleApplication
	$(DATE) > $(test_stamp)

WITH_MPI:=$(shell perl -I tools/pel tools/pel/arch.pl $(CCC) -getvariable_extra=WITH_MPI)
ifeq ($(WITH_MPI),1)
 MPIRUN:=$(shell perl -I tools/pel tools/pel/arch.pl $(CCC) -getvariable_extra=MPIRUN)
 RES_MPI_RUN=$(shell if test -x $(MPIRUN) ; then echo "1" ; fi )
 ifneq ($(RES_MPI_RUN),1)
  $(error "No valid MPIRUN variable defined : "$(MPIRUN))
 endif
endif

$(libtest):
	$(MKDIR) $(libtest)

check:$(libtest) exe_test sa_exe_test
	. $(CONFIG).sh ; cd $(libtest) ; pel test -peltest_exe ../../$(exetests) ../../$(sa_exetests) ../../SingleApplication
	. $(CONFIG).sh ; cd $(libtest) ; pel run -noverb -noautocheck ../../$(exetests) ../../admin/check.pel resu_checks

#########################
# External API          #
#########################

PETSc_VERSION:=$(shell perl -I tools/pel tools/pel/arch.pl $(CCC) -getvariable_extra=PETSc_VERSION)

PETSc_TESTS =
ifdef PETSc_VERSION
   PETSc_TESTS += ../../ExternalAPI/PETSc_$(PETSc_VERSION)
   PETSc_TESTS += ../../UnitTests/ExternalAPI/PETSc_$(PETSc_VERSION)
endif

check_PETSc: $(libtest) exe_test
   ifndef PETSc_TESTS
	@echo ; echo "--------------------------------------------------" 
	@echo "| Test failed: ExternalAPI PETSc not enabled"
	@echo "--------------------------------------------------"; echo
   else
	. $(CONFIG).sh ; cd $(libtest) ; pel test -exact ../../$(exetests) $(PETSc_TESTS)
   endif

check_MPI: $(libtest) exe_test
	. $(CONFIG).sh ; cd $(libtest) ; pel test -exact ../../$(exetests) ../../ExternalAPI/MPI ../../UnitTests/ExternalAPI/MPI 

check_METIS: $(libtest) exe_test
	. $(CONFIG).sh ; cd $(libtest) ; pel test -exact ../../$(exetests) ../../ExternalAPI/METIS_4.0.1 ../../UnitTests/ExternalAPI/METIS 

check_UMFPACK  : tttt_UMFPACK
check_Aztec    : tttt_Aztec
check_SPARSKIT : tttt_SPARSKIT

tttt_%: $(libtest) exe_test
	. $(CONFIG).sh ; cd $(libtest) ; pel test -exact ../../$(exetests) ../../ExternalAPI/$* ../../UnitTests/ExternalAPI/$*

#######
# Doc #
#######

DIRTOOLS.doc:=$(addsuffix .doc,$(DIRTOOLS))

%.doc:$(CONFIG).sh
	$(BEGIN_TARGET) $*
	. $(CONFIG).sh ; cd $(@D) ; $(MAKE) doc
	$(END_TARGET) $*

doc:pelgenerated_doc doc/Site/API_Documentation/general.htm $(DIRTOOLS.doc)
	$(END_TARGET) $@

pelgenerated_doc:classdoc exappdoc

DIR_DOC_HTML:=doc/Html

$(DIR_DOC_HTML):
	$(MKDIR) $@

classdoc:exe_test onlyclassdoc

onlyclassdoc:$(CONFIG).sh $(DIR_DOC_HTML)
	$(BEGIN_TARGET) 'CLASS DOCUMENTATION'
	. $(CONFIG).sh ; cd $(DIR_DOC_HTML) ; pel run $(CURDIR)/$(exetests) ../share/data.pel peldoc.out -noverb -noautocheck
	$(END_TARGET) 'CLASS DOCUMENTATION'

DIR_STD_DOC_HTML:=tests/doc/Html

$(DIR_STD_DOC_HTML):
	$(MKDIR) $@

exappdoc:exe_test $(CONFIG).sh $(DIR_STD_DOC_HTML) $(libtest)
	$(BEGIN_TARGET) 'DOCUMENTATION: EXAMPLES OF APPLICATION'
	. $(CONFIG).sh ; pel predoc tests/doc/description.txt tests/doc/data_doc.pel StdApp ExamplesOfApplication -Wno_unresolved
	. $(CONFIG).sh ; cd $(DIR_STD_DOC_HTML) ; pel run $(CURDIR)/$(exetests) ../data_doc.pel ../peldoc.out -noverb -noautocheck
	$(END_TARGET) 'DOCUMENTATION: EXAMPLES OF APPLICATION'

doc/Site/API_Documentation/general.htm:exe_test $(CONFIG).sh $(DIRTOOLS.all)
	$(BEGIN_TARGET) 'HTML GENERATED DOCUMENTATION'
	. $(CONFIG).sh ; cd admin ; perl ./build_html_docs ../$(exetests)
	$(END_TARGET) 'HTML GENERATED DOCUMENTATION'

##########
# Timing #
##########

timing:$(CONFIG).sh
	$(BEGIN_TARGET) 'Timing'
	. $(CONFIG).sh ; cd tests ; pel time -build_and_add pCheckTime.bd
	$(END_TARGET) 'Timing'

###########
# Dist #
###########

distfile=$(shell $(DATE) +pelicans.%d_%m_%Y.tar.gz)

dist:
	$(BEGIN_TARGET) $(distfile)
	admin/archive_version . pelicans
	$(END_TARGET) $(distfile)

cleandist:
	$(RM) $(distfile)

########
# Help #
########

help:
	@echo
	@echo "Usage: $(MAKE) CCC=<arg> <target>"
	@echo "<target>"
	@echo "    environment : build environment configuration files"
	@echo "    all         : build libraries"
	@echo "    check       : check install of built-in functionalities"
	@echo "    check_<API> : check activation of an External API with"
	@echo "                  <API> = MPI, PETSc, UMFPACK, METIS, Aztec" 
	@echo "                          or SPARSKIT"
	@echo "    clean       : remove generated files"
	@echo "    maintainer-help : help for advanced maintainers"
	@echo "<arg> : compiler system (default: gcc)"
	@echo

maintainer-help:
	@echo
	@echo "Usage: $(MAKE) CCC=<arg1> OPT_TEST=<arg2> <target>"
	@echo "<target> for standard users:"
	@echo "    environment : build environment configuration files"
	@echo "    all         : build libraries (libpel0.so, libpel1.so)"
	@echo "    check       : check install of built-in functionalities"
	@echo "    check_<API> : check activation of an External API with"
	@echo "                  <API> = MPI, PETSc, UMFPACK, METIS, Aztec" 
	@echo "                          or SPARSKIT"
	@echo "    peledit     : build PELICANS editor jni (interface C++/java)"
	@echo "    clean       : remove generated files"
	@echo "<target> for advanced users:"
	@echo "    generated   : build bison and flex generated files"
	@echo "    cleangenerated : remove bison and flex generated files"
	@echo "    extra-libraries : build libpel2.so and libpelg.so"
	@echo "    exe_test    : build executable of ExamplesOfApplication"
	@echo "    doc         : build the whole documentation"
	@echo "    classdoc    : build documentation of libraries"
	@echo "    onlyclassdoc: build documentation of libraries but not executable"
	@echo "    exappdoc   : build documentation of ExamplesOfApplication"
	@echo "    dist        : build the tarball of the current version"
	@echo "    testing     : run all self-tests"
	@echo "    timing      : build the CPU database"
	@echo "    maintainer-clean"
	@echo "      same as: clean all except bison and flex generated files"
	@echo "    maintainer-update"
	@echo "      same as: libraries extra-libraries testing"
	@echo "    maintainer-install"
	@echo "      same as: cleangenerated maintainer-clean doc dist maintainer-update"
	@echo "<arg1> : compiler system    (default: gcc)"
	@echo "<arg2> : optimization level (default: opt1)"
	@echo

FORCE:
# DO NOT DELETE
