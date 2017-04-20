#!/bin/csh -f


set dir = $cwd
cd `dirname $0`
set dirbase = $cwd
set dirsrc = $dirbase/sources
set PELGUISHOME = $dirbase
cd $dir

set main = peleditor.Peleditor
set CPFLAG = -classpath

    set CP = ${dirbase}:${dirbase}/jars/peleditor.jar
    set ld_path = $dirbase/jars

set CLASSARGS=(--noframe)
set JAVAARGS
set JAVAARGS = ($JAVAARGS -DCWD=$dir)
set JAVAARGS = ($JAVAARGS -Djava.library.path=$ld_path)
set JAVAARGS = ($JAVAARGS -Dpelguis-homedir=${PELGUISHOME})
set JAVAARGS = ($JAVAARGS -Duser.language=en)
set JAVAARGS = ($JAVAARGS -Dverbosity=SEVERE)

# Default application name:
set APPLINAME = "PelEditor"
set CONTROL
set EXE
set APPLIDIR

if ( "$1" == "-help" || "$1" == "-h" ) then
    echo 
    echo 'Usage : ' `basename $0` -h
    echo '        ' `basename $0` '[OPTIONS] [<data_file>]'

    sed -n 's/^###//p' $0
   exit 1
endif

while ("$1" != "")
    switch("$1")
	case "-lang":
	    shift
	    set JAVAARGS = ($JAVAARGS -Duser.language=$1)
	    breaksw
	case "-vv":
	    set JAVAARGS = ($JAVAARGS -Dverbosity=ALL)
	    breaksw
	case "-v":
	    set JAVAARGS = ($JAVAARGS -Dverbosity=FINE)
	    breaksw
	case "-i":
	    set JAVAARGS = ($JAVAARGS -Dverbosity=INFO)
	    breaksw
	case "-q":
	    set JAVAARGS = ($JAVAARGS -Dverbosity=SEVERE)
	    set CLASSARGS = ($CLASSARGS --noframe)
	    breaksw
	case "-control":
	    shift
	    set p = "$1"
	    set car = `echo $p | cut -c-1`
	    if ("$car" != '/') set p = $dir/$p
            if ( ! -f $p ) then
               echo "   error: control file $p not found"
               exit 1
            endif
            set CONTROL = $p
	    set JAVAARGS=($JAVAARGS -Dcontrol-filename=$p)
	    breaksw
	case "-exe":
	    shift
	    set p = "$1"
	    set car = `echo $p | cut -c-1`
	    if ("$car" != '/') set p = $dir/$p
            if ( ! -x $p ) then
               echo "   error: executable $p not found"
               exit 1
            endif
            set EXE = $p
	    set JAVAARGS = ($JAVAARGS -Dexecutable-filename=$p)
	    breaksw
	case "-application-datadir":
	    shift
	    set p = "$1"
	    set car = `echo $p | cut -c-1`
	    if ("$car" != '/') set p = $dir/$p
            if ( ! -d $p ) then
               echo "   error: application directory $p not found"
               exit 1
            endif
            set APPLIDIR = $p
	    set JAVAARGS = ($JAVAARGS -Dapplication-datadir=$1)
	    breaksw
	case "-application-name":
	    shift
            set APPLINAME = $1
	    set JAVAARGS = ($JAVAARGS -Dapplication-name=$1)
	    breaksw
	default:
	    set CLASSARGS=($CLASSARGS "$1")
    endsw
shift
end

echo "Launching $APPLINAME"
echo "     PELGUISHOME: $PELGUISHOME"
if( "$APPLIDIR" != "" ) echo "     APPLIHOME: $APPLIDIR"
if( "$CONTROL" != "" )  echo "     control file: $CONTROL"
if( "$EXE" != "" )      echo "     executable: $EXE"
java -Xmx400m -Djava.library.path=$ld_path $JAVAARGS $CPFLAG ${APPLIDIR}:$CP $main $CLASSARGS


### 
###  Options:
###     -exe <exe>
###        set executable for running
###     -control <controler_file>
###        set controler file
###     -application-datadir <directory>
###        set application directory
###     -application-name <name>
###        set application name
###     -lang <lang>
###        change default language (ex. -lang fr)
###
###  Verbosity options:
###     -q
###        quiet mode: only severe errors are displayed
###     -i
###        active verbose flag
###     -v
###        active some debug flags
###     -vv
###        active all debug flags
### 
