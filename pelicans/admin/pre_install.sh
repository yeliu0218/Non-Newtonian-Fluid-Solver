#!/bin/sh -e
set -x

if [ -z "$PELICANSHOME" ]
then
  echo "No PELICANSHOME defined"
  exit
fi

##############################
# PELICANS user library extraction
# Extract part of PELICANS library and prepare it
#  to be installed on a remote system in user mode
##################################################


##################################################
#
# Configurable items
##
RELEASE=HEAD
INSTALLDIR=pelicans
DATE=`date`
VERSIONFILE=PELbase/src/PEL_Version.cc
WASTE=$HOME/Waste
PACKAGES="Geometry LinearAlgebra PDEsolver PELbase"
UNWANTEDDIR="doc/Site/HowTo doc/Site/local doc/Site/man doc/Site/pel_CLI doc/Html/src tools bin Makefile"
STD_API="HelloWorld FiniteElementInterpolation Helmholtz_Galerkin Poisson_Mortar Stokes_StabilizedGalerkin Stefan_ALE FSflow_CharacteristicGalerkin ConvectionDiffusion_QUICK NavierStokes_StaggeredGrid ShallowWater_Godunov visu "

##################################################
#
# And now, ladies an gentlemen, GO !
#
MAKE=make
if gmake -version > /dev/null
then
  MAKE=gmake
fi
# Export pelicans from CVS
##
if [ ! -d $INSTALLDIR ]
then
  cvs -d $CVSROOT export -r $RELEASE -d $INSTALLDIR pelicans
  \rm -fr $WASTE/sem_*
fi

# go to pelicans
cd $INSTALLDIR
if [ ! -d $WASTE ]
then
  mkdir $WASTE
fi

# Set date
sem=$WASTE/sem_date
if [ ! -f $sem ] 
then
  mv $VERSIONFILE $VERSIONFILE.old
  subst_date="s/DATE/$DATE/"
  subst_rel="s/RELEASE/$RELEASE/"
  cat $VERSIONFILE.old | sed "$subst_date" | sed "$subst_rel" > $VERSIONFILE
  date > $sem
fi

# make include generated and doc
sem=$WASTE/sem_include
if [ ! -f $sem ]
then
 $MAKE generated
 $MAKE doc
 cp -r -p include include.new
 \rm -fr include
 mv include.new include  
 date > $sem
fi

# make source directory
sem=$WASTE/sem_source
libsrcdir=library
if [ ! -f $sem ]
then
 mkdir $libsrcdir
 mv $PACKAGES $libsrcdir
 date > $sem
fi

# recover examples of applications 
sem=$WASTE/sem_stdapp
if [ ! -f $sem ]
then
 mv ExamplesOfApplication ExamplesOfApplication.old
 mkdir ExamplesOfApplication
 cd ExamplesOfApplication.old
 mv $STD_API ../ExamplesOfApplication
 cd ..
 \rm -fr ExamplesOfApplication.old
 date > $sem
fi

# remove unwanted packages
sem=$WASTE/sem_del
if [ ! -f $sem ]
then
  for p in $UNWANTEDDIR
  do
    if [ -d $p -o -f $p ]
    then
      \rm -fr $p
    fi
  done 
  find doc -name "*.ps" -exec \rm {} \;
  date > $sem
fi



