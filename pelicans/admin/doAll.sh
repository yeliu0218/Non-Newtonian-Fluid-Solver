#! /bin/sh -e
# Run all tests.

if [ -z "$PELICANSHOME" ]
then
    echo "Environement variable PELICANSHOME not defined"
    exit
fi

EXE="$PELICANSHOME"/tests/lib/opt2SunOSCC/exe
PRUN="$PELICANSHOME"/bin/pRun
if [ ! -x "$EXE" ]
then
    echo "Executable $EXE not available !"
    exit
fi

DIR=$PELICANSHOME"/ExamplesOfApplication"

for oneTest in `find $DIR -name data.pel -print`
do
   cd `dirname $oneTest`
   extra=""
   if [ -f "pTests.cfg" ]
   then
    extra=`grep "CHECKING LEVEL" "pTests.cfg"|cut -d"=" -f2`
   fi
   $PRUN $EXE data.pel resu $extra
done

