#
#  Copyright 1995-2010 by IRSN
#
#  This software is an application framework, with a set of integrated
#  reusable components, whose purpose is to simplify the task of developing
#  softwares of numerical mathematics and scientific computing.
#
#  This software is governed by the CeCILL-C license under French law and
#  abiding by the rules of distribution of free software. You can use, modify
#  and/or redistribute the software under the terms of the CeCILL-C license
#  as circulated by CEA, CNRS and INRIA at the following URL
#  "http://www.cecill.info".
#
#  As a counterpart to the access to the source code and rights to copy,
#  modify and redistribute granted by the license, users are provided only
#  with a limited warranty and the software's author, the holder of the
#  economic rights, and the successive licensors have only limited liability.
#
#  In this respect, the user's attention is drawn to the risks associated
#  with loading, using, modifying and/or developing or reproducing the
#  software by the user in light of its specific status of free software,
#  that may mean that it is complicated to manipulate, and that also
#  therefore means that it is reserved for developers and experienced
#  professionals having in-depth computer knowledge. Users are therefore
#  encouraged to load and test the software's suitability as regards their
#  requirements in conditions enabling the security of their systems and/or
#  data to be ensured and, more generally, to use and operate it in the same
#  conditions as regards security.
#
#  The fact that you are presently reading this means that you have had
#  knowledge of the CeCILL-C license and that you accept its terms.
#

##
# MAIN

package main ;

use strict ;

use Getopt::Long ;
use Pod::Usage ;
use File::Find ;
use Text::Wrap ;

use Util ;
use Defs ;

#-------------------------------------------------------------------------
# Global variables
#-------------------------------------------------------------------------
my $verbose = 0 ;
my $successful_test = 1 ;
my $default_options = '' ;
my $test_dir = '' ;


#-------------------------------------------------------------------------
#  MAIN starts here
#-------------------------------------------------------------------------

my $dir = Util::absolute_pathname( File::Spec->curdir() ) ;

my $help = 0 ;
my $man  = 0 ;
my @extra=() ;
my $build_pattern  = '' ;
my $verify_pattern = '' ;
my $build_then_verify_pattern = '' ;
my $exe_test = '' ;
my $exact = 0 ;
my $dbl_eps = undef ;
my $dbl_min = undef ;
my $result = GetOptions (
     'help'       => \$help
   , 'verbose'    => \$verbose
   , 'man'        => \$man
   , 'build_pattern=s'  => \$build_pattern
   , 'verify_pattern=s' => \$verify_pattern
   , 'build_then_verify_pattern=s' => \$build_then_verify_pattern
   , 'test_directory=s' => \$test_dir
   , 'Cpost'      => sub { $default_options .= " -post"    ; }
   , 'Call'       => sub { $default_options .= " -all"     ; }
   , 'peltest_exe=s'  => \$exe_test
   , 'exact'      => \$exact
   , 'dbl_eps=s'  => \$dbl_eps
   , 'dbl_min=s'  => \$dbl_min
   ) ;

pod2usage( -verbose => 1 ) if $help  ;
pod2usage( -verbose => 2 ) if $man   ;
if( scalar @ARGV < 2 || !$result ) {
  pod2usage( -verbose => 0 ) ;
}
my( $exe, @dirs_of_args ) = @ARGV ;

if( ! $exe_test ) { $exe_test = $exe ; }

if( $build_pattern && $verify_pattern ||
    $build_pattern && $build_then_verify_pattern ||
    $verify_pattern && $build_then_verify_pattern ) {
  my $mesg = "options --build_pattern, --verify_pattern and --build_then_verify_pattern" ;
  $mesg .= "are mutually exclusive" ;
  pod2usage( $mesg ) ;
}

if( $build_pattern ) {
  $default_options .= " -pattern_build " ;
  $default_options .= Util::absolute_pathname( $build_pattern ) ;
}

if( $verify_pattern ) {
  Defs::Error( "Invalid file: $verify_pattern" ) unless ( -r $verify_pattern );
  $default_options .= " -pattern_verify " ;
  $default_options .= Util::absolute_pathname( $verify_pattern ) ;
}

if( $build_then_verify_pattern ) {
  $default_options .= " -pattern_build_then_verify " ;
  $default_options .= Util::absolute_pathname( $build_then_verify_pattern ) ;
}

if( $dbl_eps && $dbl_min ) {
  if( $exact ) {
     pod2usage('options --exact and --dbl_eps,--dbl_min are incompatible');
  }
  $default_options .= " -dbl_eps ".$dbl_eps ;
  $default_options .= " -dbl_min ".$dbl_min ;
} else {
  if( !( !$dbl_eps && !$dbl_min ) ) {
     pod2usage('options --dbl_eps and --dbl_min should be used together');
  }
}

if( $exact ) {
  $default_options .= " -exact " ;
}

$exe = Util::absolute_pathname( $exe ) ;
Defs::Error( "Invalid executable: $exe" ) unless ( -X $exe ) ;

#--- search for list of tests to run
my $tdirs = "" ;
for ( @dirs_of_args ) {
  if( ! -d $_ ) { Defs::Error( "Unknown directory : $_" ) ; }
  $tdirs .= " $_" ;
}

my $cmd = "$exe_test -A peltest -executable $exe $default_options -output_file resu -test_directories $tdirs" ;

#--- resu directories
if( $test_dir ) {
   my $odirs = "" ;
   my @adirs = @dirs_of_args ;
   for( @adirs ) { $_ = Util::absolute_pathname( $_ ) ; }
   my $common_root = Util::common_root( @adirs ) ;
   for ( @adirs ) {
      my $d = $_ ;
      $d =~ s/$common_root\/// ;
      $d = File::Spec->join( $test_dir, $d ) ;
      $odirs .= " $d" ;
   }
   $cmd .= " -output_directories $odirs" ;
}

#--- mpirun
my @script = Defs::script( "arch" ) ;
my $mpirun_line = "@script -getvariable_extra MPIRUN" ;
my $mpirun = `$mpirun_line` ;
if( -x "$mpirun" ) { $cmd .= " -mpirun $mpirun" ; }

#--- create report file
my $suffixe=`date +.%d_%m_%y.%Hh%M` ;
chomp $suffixe ;
my $report = "report$suffixe" ;
my $i = 1 ;
while( -e $report ) {
   $i++ ;
   $report = "report$suffixe.$i" ;
}
$cmd .= " -verbose" if( $verbose ) ;
$cmd .= " | tee $report" ;

# print "COMMANDE : $cmd \n" ;
if( $verbose ) {
   print "running: $cmd\"\n" ;
}
exec $cmd ;
exit(0) ;


##
# POD Documentation
#
__END__

=head1 NAME

test - comparison between runs of a PELICANS-based application (for regression testing)

=head1 SYNOPSIS

pel test [-help|-man]

pel test [options...] F<exe> F<dirs>

pel test -build_pattern F<filename> F<exe> F<dirs>

pel test -verify_pattern F<filename> F<exe> F<dirs>

pel test -build_then_verify_pattern F<filename> F<exe> F<dirs>

=head1 DESCRIPTION

Comparing in details two runs of an application is important for
the sake of non regression or installation testing.

Given a hierarchy of directories containing reference runs
of a given application, C<pel test> will perform the following actions:

=over

=item 1.

create, in the
working directory, a subdirectory
in which that hierarchy is duplicated;

=item 2.

in all subdirectories of the duplicate hierarchy, run the application
with the associated reference data file (called F<data.pel>);

=item 3.

compare the results of this run with the reference results.

=back

On completion, the conclusions of all these comparisons are recorded
in a report file located in the current directory.

The essentials of C<pel test> tasks are forwarded to the PELICANS-based
application "peltest" (contained in the executable specified by the
B<-peltest_exe> option if any, or by the argument F<exe>).

Further details are given below.

=over

=item .

If a test failure is pronounced, the character sequence
"Test failed" will appear in the report file.

=item .

Runs are performed with the commands like:

sequential: F<exe> F<data.pel> -v [opts...] > F<resu>

parallel:   F<mpirun> [mpi_opts...] F<exe> F<data.pel> -v -o F<resu> [opts...]

where the options C<opts...> and C<mpi_opts...> are determined
for all runs by the calling options of C<pel test> (see below)
and for one specific run by a possible file F<config.pel> in
the reference directory of that run.

A file F<config.pel> containing:

  MODULE test_config
    run_options = vector( "...", "..." )
  END MODULE test_config

leads to the addition in C<opts...> of all items in the StringVector
(only for the run associated to the directory
containing the considered file F<config.pel>).

A file F<config.pel> containing:

  MODULE test_config
    mpi_options = vector( "...", "..." )
  END MODULE test_config

switches to parallel execution and leads to the addition in C<mpi_opts...> of
all items in the StringVector (only for the run associated to the directory
containing the considered file F<config.pel>).
An additional optional StringVector data of keyword C<mpi_machinefile>
can be used to specify the list of possible machines to run on
(this data is written on a temporary machine file transmitted to F<mpirun>
via the option C<-machinefile>).

=item .

The exit code is tested. The test failure is pronounced
if it is non zero, unless it exists a file F<config.pel>, stored in
the reference directory, containing:

  MODULE test_config
    failure_expected = true
  END MODULE test_config

in which case success might be pronounced if exit code is non zero and
one or more produced files called F<expected.err*> are identical
(same name and same content) to those present in the reference directory.

=item .

The test failure is pronounced if the F<resu> file has not been produced.

=item .

All files that have been produced during the run (other than F<resu>) 
are compared (as described below) with the reference ones (that must exist).
This comparison might be avoided from some particular file of a particular
run if the associated reference directory stores a file F<config.pel>
containing:

  MODULE test_config
    files_to_ignore = vector( "...", "..." )
  END MODULE test_config

The items of the StringVector correspond to files produced during the run for
which no comparison will be performed.

=item .

The comparison method between the reference and produced files
depends on the format of these files.

There are three "native" formats understood by PELICANS:
the format called GENE, for C<TIC> postprocessing; the
format called PEL, for Hierarchical Data Structures with the PELICANS format;
the format called CSV, for comma separated values.

The files with format GENE, PEL or CSV are compared to the reference ones
with the PELICANS-based application "pelcmp" (contained in the
executable specified by the B<-peltest_exe> option if any, or by the
argument F<exe>). If they are not identical, the comparison results
are recorded in the report file (the test failure is not pronounced since
differences may be acceptable, depending on the use case).

The other files are compared line by line with the reference ones.
If they are not the same, the test failure is pronounced.

=item .

The format of a file, say F<save.zzz>, is determined as follows.
It can be specified via a configuration file F<config.pel> stored
in the reference directory:

  MODULE test_config
    MODULE PEL_Comparator
      MODULE xxx               // xxx is a non significant name
        filename = "save.zzz"
        format = "CSV"         // either "GENE", "PEL" of "CSV"
      END MODULE xxx
    END MODULE PEL_Comparator
  END MODULE test_config

If such a specification is absent, the format is identified on the
basis of a motif appearing in the file name: F<.gene> gives the GENE
format, F<.pel> gives the PEL format, F<.csv> gives the CSV format.

=item .

The files with format GENE or PEL contain data identified by keywords.
The comparison might ignore some of these data if the associated
reference directory stores a file F<config.pel>
containing:

  MODULE test_config
    MODULE PEL_Comparator
      MODULE xxx               // xxx is a non significant name
        filename = "save.zzz"  // file with format GENE or PEL
        ignore_data = vector( "...", "..." )
      END MODULE xxx
    END MODULE PEL_Comparator
  END MODULE test_config

The items of the StringVector correspond to keywords of data that
should be ignored during the comparison.

=item .

The floating point values contained in the reference and produced files
(with format GENE, PEL or CSV) are compared with PEL::double_equality.
The last two arguments of
this member function are respectively called a_dbl_eps (a
kind of tolerance on relative errors) and a_dbl_min (a lower bound under which
values are undistinguishable from zero).

By default, a_dbl_eps and a_dbl_min are equal to zero (which means that
comparisons without any tolerance are performed).
They can be given other values either globally (for all runs)
using the C<-dbl_eps> and C<-dbl_min> options, or for a specific
run via a file F<config.pel> in the reference directory of that
run.

For instance, a file F<config.pel> containing:

  MODULE test_config
     MODULE PEL_Comparator
        MODULE xxx                       // xxx is a non significant name
           filename = "save.csv"
           MODULE double_comparison
              dbl_min = 1.e-10
              dbl_eps = 1.e-8
           END MODULE double_comparison
        END MODULE xxx
     END MODULE PEL_Comparator
  END MODULE test_config

will set a_dbl_min=1.e-10 and a_dbl_eps=1.e-8 for comparisons between
the floating point values of the files F<save.csv>.

Note that the command line options C<-dbl_eps> and C<-dbl_min>
always overread the options stated in the files F<config.pel>.
Moreover the line option C<-exact> can be used to ignore any
setting of a_dbl_min and a_dbl_max in the F<config.pel> files.

=back

When the C<-verify_pattern> option is activated, the behavior
of C<pel test> is slightly different: the only test performed
is the conformance of the reference data file F<data.pel> with
the given pattern file.

=head1 ARGUMENTS

=over

=item B<exe>

Name of the executable of the PELICANS-based application to run.

=item B<dirs>

List of the directories defining the reference runs. Any subdirectory
of an item of F<dirs> containing a file F<data.pel> is considered
by C<pel test> as a definition of a reference run whose data file
is F<data.pel>. This subdirectory must contain the reference
version of all the files produced when calling F<exe> with that
data file. It might also contain (see above) a file called F<config.pel>
and, more rarely, files called F<expected.err*>.

=back

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-Cpost>

Call C<pel run> with this option for all runs.

=item B<-Call>

Call C<pel run> with this option for all runs.

=item B<-build_pattern> F<filename>

Call C<pel run> with this option for all runs.

=item B<-verify_pattern> F<filename>

Do not perform the runs, but instead use the PELICANS-based application
"check" (contained in the argument F<exe>) to check the conformity
of all reference data file F<data.pel> with the pattern file F<filename>.

=item B<-build_then_verify_pattern> F<filename>

Call C<pel run> with the option B<-build_pattern> F<filename> 
for all runs and then check the conformity
of all reference data file F<data.pel> with the created pattern file
F<filename> (equivalent to two calls of C<pel test>
with successively the B<-build_pattern> and the B<-verify_pattern> options).

=item B<-test_directory> F<dirname>

Duplicate the hierarchy of directories containing the reference runs
in the subdirectory F<dirname> of the working directory,
and run the application in the subdirectories of F<dirname> for further
result comparison with the reference runs (default: F<PELICANS_TEST>).

=item B<-peltest_exe> F<texe>

Specify the executable containing the "peltest" and "pelcmp" applications.
Default is the argument F<exe> itself.

=item B<-dbl_eps> F<eps>

Specify the a_dbl_eps argument in calls to PEL::double_equality
when comparing floating point values. This option is only
significant for files with format PEL, CSV or GENE.

=item B<-dbl_min> F<min>

Specify the a_dbl_min argument in calls to PEL::double_equality
when comparing floating point values. This option is only
significant for files with format PEL, CSV or GENE.

=item B<-exact>

Always perform comparisons between floating point values
without any tolerance, whatever settings of a_dbl_eps
and a_dbl_eps in files F<config.pel>.

=back

=head1 EXAMPLES

=over

=item C<pel test ../bin/exe ../RegressionTests>

Run executable F<exe> located in the directory F<../bin>
with all data files F<data.pel> contained in the subdirectories
of F<../RegressionTests> and compare the results with the reference ones.
Create a report file in the current directory recording the conclusions
of all comparisons.


=item C<pel test -build_pattern pat.pel ../bin/exe ../Tests>

Same as before, with in addition the learning and storage of the requested
structure of the data files in F<pat.pel>.

=item C<pel test -verify_pattern pat.pel ../bin/exe ../Appli>

Check the conformance with F<pat.pel> of all files F<data.pel> contained 
in a subdirectory of F<../Appli>, and record the conclusions in a report 
file in the current directory.

=back

=cut
