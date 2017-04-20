#! /bin/sh

RM="rm -fr"

if [ ! -d "$PELICANSHOME" ]
then
  echo "No valid PELICANSHOME variable defined !"
  exit
fi

$RM $PELICANSHOME/tests 

# Clean up objects and depedencies.
# Two passes to not remove SunOSCC/dbg objects used in debugger
find "$PELICANSHOME" -name dbg -prune -o \( -name "*.o" -o -name "*.d" -o -name SunWS_cache \) -print -exec $RM {} \;
find "$PELICANSHOME" -name SunOSCC -prune -o \( -name "*.o" -o -name "*.d" -o -name SunWS_cache \) -print -exec $RM {} \;

#LaTex temporary and generated postscipts
find "$PELICANSHOME" \( -name "*.log" -o -name "*.toc" -o -name "*.aux" -o -name "*.dvi" \) -print  -exec $RM {} \;
