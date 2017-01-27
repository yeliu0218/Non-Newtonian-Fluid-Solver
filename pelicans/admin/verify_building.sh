#! /bin/sh 
# This script is aimed to verify that all targets implied in building
# are successful

verify_peldoc() {

if [ ! -f $peldoc_output ]
then
  error_message=$error_message"No file $peldoc_output found.\n"
fi
if [ -s $peldoc_output ]
then
  error_message=$error_message"$peldoc_output : "
  error_message=$error_message`cat $peldoc_output |wc -l`
  error_message=$error_message" error(s) found by peldoc."
fi
}

#---------------------------------------------------------------------
#   MAIN
#---------------------------------------------------------------------
test_string="failed"
test_goal="Test results are strictly compared with data base"
created_files=""
while( [ $# -ge 1 ] )
do
   case "$1" in
      "-h")
         echo "PELICANS : `basename $0` utility"
         echo "           located in `dirname $0`"
         echo "usage `basename $0` [-test_difference_allowed] [list_of_files]"
         exit 1
         ;;
      "-test_difference_allowed")
         test_string="Test failed"
         test_goal="Normal execution of the tests are checked"
         shift
         ;;
      *)
         created_files="$created_files $1"
         shift
   esac
done

cd $PELICANSHOME
error_message=""
ECHO="/usr/bin/echo"
if [ ! -x $ECHO ]
then
  ECHO="echo -e"
fi

# Test checking
echo " "
report=`ls $PELICANSHOME/tests/*/report*`
if [ -z "$report" ]
then
  error_message=$error_message"No report file generated.\n"
else
  echo $test_goal
  for report in $report
  do
    echo "   checking file : $report"
    if ( grep "$test_string" $report 2>&1 ) > /dev/null 
    then
      error_message=$error_message"Some tests failed in :\n  "$report"\n"
    fi
  done
fi

# Timing
timing_output="tests/pCheckTime.bd"
timing_ref="admin/pCheckTime_ref.bd"
if [ ! -f $timing_output ]
then
  error_message=$error_message"No timing report "$timing_output"\n"
else
  if( pel time -compare $timing_ref $timing_output )
  then
    echo " "
    echo "Test execution times are checked"
  else
    error_message=$error_message"Time comparison failed\n"
    error_message=$error_message
  fi
fi

# peldoc report checking
echo " "
echo "Documentation is checked"
peldoc_output="doc/Html/peldoc.out"
echo "   checking file : $peldoc_output"
verify_peldoc
peldoc_output="tests/doc/peldoc.out"
echo "   checking file : $peldoc_output"
verify_peldoc

# Verify files
echo " "
for file in $created_files
do
   if [ -f $file ]
   then
      echo "file $file : ok"
   else
      echo "file $file : not built"
      error_message=$error_message"file $file : not built.\n"
   fi
done

$ECHO "\n\n******************************"
$ECHO "* verify_building.sh reports : "
$ECHO "******"
if [ -n "$error_message" ]
then
  $ECHO $error_message
  $ECHO "Building failed : let correct previous errors and try it again !"
  exit 1
else
  $ECHO "Congratulations : Building were successful !"
  exit 0
fi
