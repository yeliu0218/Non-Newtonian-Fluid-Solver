#! /bin/sh
# Compile all PELICANS to be used with tcov tool.
. ../../bin/init.sh

if [ ! -d src ]
then
    mkdir src
    pel --mirror src 
fi

if [ ! -d lib ]
then
    mkdir lib
fi
pel --depend --profile opt2 lib src
pel --build --exe lib

pel --test --build_pattern pattern.pel lib/exe
pel --test --verify_pattern pattern.pel
tmp=explore_tcov.tmp_dir

if [ -d $tmp ]
then
    \rm -fr $tmp
fi
mkdir $tmp
set -x
MYSYSTEM=`uname`
if [ "$MYSYSTEM" = "SunOS" ]
then
 here=`pwd`
 prof=exe.profile
 cat `find ./ -name tcovd` > toto

 mv toto $tmp/tcovd

 cd $tmp

 files=`ls ../src/*.cc`
 tcov -a -x . $files

 grep "Percent of the file executed" *.tcov | sort -n -k 2,2 > $here/tcov_sheet
fi
if [ "$MYSYSTEM" = "Linux" ]
then
  echo "coucou"
  LIST=`find ./ -name gmon.out`
  gprof -A lib/exe $LIST > gprof_sheet
fi
