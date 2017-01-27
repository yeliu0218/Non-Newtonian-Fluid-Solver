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
my $default_options = '-noautocheck' ;
my $test_dir = "PELICANS_TESTS" ;

#-------------------------------------------------------------------------
sub create_dir_then_go_and_clean {
  my @xxx = File::Spec->splitdir( shift @_ ) ;
  for (@xxx) {
     if ( ! -d $_ ) {
         mkdir( $_ ) || Defs::Error( "cannot mkdir $_ ($!)" ) ; 
     }
     chdir( $_ ) || Defs::Error( "cannot chdir $_ ($!)" ) ;
  }
  #-- delete files
  opendir( DIR, File::Spec->curdir() ) ;
  foreach ( readdir( DIR ) ) {
     next if( $_ eq '.' ) ;   # to avoid the '.' directory
     next if( $_ eq '..' ) ;  # to avoid the '.' directory
     unlink( $_ ) ||  Defs::Error( "cannot delete $_ ($!)" ) ;
  }
}

#-------------------------------------------------------------------------
# IO procedures
#-------------------------------------------------------------------------

sub begin {
  my $test_name = shift @_ ;
  print "--------------------------------------------\n" ;
  print "| $test_name \n" ;
  print "--------------------------------------------\n" ;
  $successful_test = 1 ;
}

sub output {
  my $arg = shift ;
  my @row = split /\n/, $arg ;
  chomp @row ;
  map { print "| $_\n" ; } @row ;
  1;
}

sub decide_failure {
  $successful_test = 0 ;
}

sub decide_undetermined {
  if( $successful_test == 1 ) {
    $successful_test = 2 ;
  }
}

sub end {
  if( $successful_test == 1 ) {
     print ( "| Test is successful\n" ) ;
  }
  elsif ( $successful_test == 0 ) {
     print ( "| Test failed\n" ) ;
  }
  elsif ( $successful_test == 2 ) {
     print ( "| Test success ?... to be analysed... \n" ) ;
  }
  print "--------------------------------------------\n" ;
}

#-------------------------------------------------------------------------
# Check for remaining objects
#-------------------------------------------------------------------------
sub check_objects {
  my $file = shift @_ ;
  open HANDLE, "$file" ;
  my $remain = grep /number of remaining P-objects/, <HANDLE> ;
  if( $remain ) { 
    output( "Object lifetime management failed" ) ;
    decide_failure() ; 
  }
}

#-------------------------------------------------------------------------
# comparison between two files using PEL_Application "pelcmp"
#-------------------------------------------------------------------------
sub compare_with_pelcmp {
  my $file_ref = shift @_ ;
  my $file_new = shift @_ ;
  my $exec = shift @_ ;

  my $cmd = "$exec -A pelcmp $file_ref $file_new" ;
  open( PROC, "$cmd|" ) ;
  my $report ;
  while( <PROC> ) { $report .= $_ ; }
  close( PROC ) ;
  if( $report ) {
    output( "strict comparison of $file_new with reference failed" ) ;
    output( $report ) ;
    decide_undetermined() ;
  }
  $report ;
}

#-------------------------------------------------------------------------
# comparison between two files using UNIX "cmp"
#-------------------------------------------------------------------------
sub compare_with_cmp_of_unix {
  my $file_ref = shift @_ ;
  my $file_new = shift @_ ;

  my $cmd = "cmp -s $file_new $file_ref" ;
  my $report = system( "$cmd" ) ;
  if( $report ) {
    output( "$file_new comparison failed" ) ;
    decide_undetermined() ;
  }
  $report ;
}

#-------------------------------------------------------------------------
# Run of a single test
#-------------------------------------------------------------------------
sub run_one_test {
  my $exe = shift @_ ;
  my $test = shift @_ ;
  my $common_root = shift @_ ;

  #--- Data filename recovering
  my $data = File::Spec->join( $test, "data.pel" ) ;
  Defs::Error( "Bad directory $test" ) unless( -d $test && -f $data ) ;

  #--- Test directory building
  my $subfix = $test ;
  $subfix =~ s/$common_root\/// ;
  begin( $subfix ) ;

  my $output = File::Spec->join( $test_dir, $subfix ) ;
  create_dir_then_go_and_clean( $output ) ;

  #-- Check for additional test configuration
  my @files_to_ignore = () ;
  my $pel_run_options = '' ;
  my $objects_checking = 1 ;
  if( open CFG, "$test/pel_test.pl" ) {
    if( $verbose ) { print STDOUT "Found a configuration file\n" }
    while ( my $line = <CFG> ) {
      chomp $line ;
      if( $verbose ) { print STDOUT "Evaluating : $line\n" }
      eval $line ;
      if( $@ ) { Defs::Error( "Invalid line in \"pel_test.pl\":\n--->$line" ) }
    }
  }
  else {
    if( $verbose ) { print STDOUT "no configuration file\n" }
  }

  #--- Run command
  my $level = $default_options." ".$pel_run_options ;
  my $resu = "resu" ;
  my @script = Defs::script( "run" ) ;
  my $cmd = "@script $level $exe $data $resu" ;
  if( $verbose ) { print STDOUT wrap ("Executing : ", "\t", "$cmd\n" ) } ;
  my $ret = system( $cmd ) >> 8 ;
  if( $verbose ) { print STDOUT "Execution result : $ret\n" } ;

  #--- Post checking
  my @g = glob "expected.err*" ;
  if ( $ret && ! ( scalar @g ) ) {
    output( "Test execution failed : (exit code $ret)" ) ;
    decide_failure() ;
  }

  # if any file of @save_list exists in the reference directory,
  # it has to be produced during the test execution
  my @save_list = ( 'save.pel', 'save.gene' ) ;
  for( @save_list ) {
    if( ( -f File::Spec->join( $test, $_ ) ) && ( ! -f $_ ) ) {
       output( "No $_ produced" ) ;
       decide_failure() ;
    }
  }

  my %ignore_file = () ;
  for( @files_to_ignore ) { $ignore_file{$_} = 1 } ;
# Test for parallel execution
  if( ( ! -f $resu ) && ( -f "$resu.0" ) )  { $resu = "$resu.0" ; }
  if( ! -f $resu ) {
    # the test execution should have produced a $resu file
    output( "No result file $resu" ) ;
    decide_failure() ;
  }
  else {
    if( $objects_checking ) { check_objects( $resu ) ; }

    # all produced files should compare successfully with the reference ones
    opendir( DIR, File::Spec->curdir() ) || die "Unable to opendir ($!)\n" ;
    foreach ( readdir( DIR ) ) {
      next if( $_ eq '.' || $_ eq '..' ) ;   # to avoid '.' and '..'
      next if( $_ =~ /resu\.*/ ) ;
      if( $ignore_file{$_} ) {
         $ignore_file{$_} = 28 ;
         next ;
      }
      my $file_new = $_ ;
      my $file_ref = File::Spec->join( $test, $_ ) ;
      if( ! -f $file_ref ) {
         output( "Reference file $file_ref does not exist" ) ;
         decide_failure() ;
      }
      else {
        if( $file_new =~ /.*\.pel$/ || $file_new =~ /.*\.gene$/ ) {
           compare_with_pelcmp( $file_ref, $file_new, $exe ) ;
        }
        else {
          compare_with_cmp_of_unix( $file_ref, $file_new ) ;
        }
      }
    }
    closedir( DIR ) || die "Unable to closedir ($!)\n" ;

    while( my ($key,$val) = each %ignore_file ) {
      if( $val != 28 ) {
        output( "the file \"$key\", quoted in \"pel_test.pl\", has not been produced" ) ;
        decide_failure() ;
      }
    }
  }
  end ;
}

#---------------------------------------------------------------------------
sub verify_one_data {
  my $pattern = shift @_ ;
  my $exe = shift @_ ;
  my $test = shift @_ ;
  my $common_root = shift @_ ;

  #--- Data filename recovering
  my $data = File::Spec->join( $test, "data.pel" ) ;
  Defs::Error( "Bad directory $test" ) unless( -d $test && -f $data ) ;
  chdir( $test ) || Defs::Error( "cannot chdir : $test" ) ;

  #--- Test directory building
  my $subfix = $test ;
  $subfix =~ s/$common_root\/// ;
  $subfix = "verifying data.pel of : $subfix" ;
  begin( $subfix ) ;

  #--- Run command
  my $cmd = "$exe -A check -s $pattern data.pel" ;
  open( PROC, "$cmd|" ) ;
  my $report ;
  while( <PROC> ) { $report .= $_ ; }
  close( PROC ) ;
  if( $report ) {
     output( $report ) ;
     decide_failure() ;
  }
  end ;
}

#-------------------------------------------------------------------------
#  MAIN starts here
#-------------------------------------------------------------------------

my $dir = Util::absolute_pathname( File::Spec->curdir() ) ;

my $help = 0 ;
my $man  = 0 ;
my @extra=() ;
my $build_pattern  = '' ;
my $verify_pattern = '' ;
my $result = GetOptions (
     'help'     => \$help
   , 'verbose'  => \$verbose
   , 'man'      => \$man
   , 'build_pattern=s'  => \$build_pattern
   , 'verify_pattern=s' => \$verify_pattern
   , 'test_directory=s' => \$test_dir
   , 'Cpost'    => sub { $default_options .= " -Cpost"    ; }
   , 'Call'     => sub { $default_options .= " -Call"     ; }
   ) ;

pod2usage( -verbose => 1 ) if $help  ;
pod2usage( -verbose => 2 ) if $man   ;
if( scalar @ARGV < 2 || !$result ) {
  pod2usage( -verbose => 0 ) ;
}
my( $exe, @dirs_of_args ) = @ARGV ;

if( $build_pattern && $verify_pattern ) {
  my $mesg = "options --build_pattern and --verify_pattern " ;
  $mesg .= "are mutually exclusive" ;
  pod2usage( $mesg ) ;
}

if( $build_pattern ) {
  $default_options .= " -build_pattern " ;
  $default_options .= Util::absolute_pathname( $build_pattern ) ;
}

if( $verify_pattern ) {
  Defs::Error( "Invalid file: $verify_pattern" ) unless ( -r $verify_pattern );
  $verify_pattern = Util::absolute_pathname( $verify_pattern ) ;
}

$exe = Util::absolute_pathname( $exe ) ;
Defs::Error( "Invalid executable: $exe" ) unless ( -X $exe ) ;

#--- search for list of tests to run
for ( @dirs_of_args ) {
  if( ! -d $_ ) { Defs::Error( "Unknown directory : $_" ) ; }
}
my @tests = () ;

File::Find::find({ "follow" => 1, 
                   "wanted" => sub { /^data\.pel$/ &&
				     push @tests, $File::Find::dir }},
	          @dirs_of_args ) ;

if(!scalar(@tests)) { Defs::Error( "No test found" ) ;}
for( @tests ) {
  $_ = Util::absolute_pathname( $_ ) ;
}
@tests = sort @tests ;
my $common_root = Util::common_root( @tests ) ;
if( $verbose ) {
  map { print "$_ : has \"data.pel\"\n" ; } @tests ;
  print "$common_root : common root to all directories with \"data.pel\"\n" ;
}

#--- create report file
my $suffixe=`date +.%d_%m_%y.%Hh%M` ;
chomp $suffixe ;
my $report = "report$suffixe" ;
my $i = 1 ;
while( -e $report ) {
   $i++ ;
   $report = "report$suffixe.$i" ;
}
open( OUT, ">$report" ) or Defs::Error( "cannot create : $report" )  ;
select OUT ;

#--- do the tests
print "executable : $exe\n" ;
print "tests      : $common_root\n" ;
for ( @tests ) {
  chdir( $dir ) || Defs::Error( "cannot chdir : $dir" ) ;
  if( $verify_pattern ) {
    verify_one_data( $verify_pattern, $exe, $_, $common_root ) ;
  }
  else {
    run_one_test( $exe, $_, $common_root ) ;
  }
}

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

Further details are given below.

=over

=item .

If a test failure is pronounced, the character sequence
"Test failed" will appear in the report file.

=item .

Runs are performed with the command:
 pel run -noautocheck [opts...] F<exe> F<data.pel> F<resu>
where the options C<[opts...]> are determined
for all runs by the calling options of C<pel test> (see below)
and for one specific run by a possible file F<pel_test.pl> in
the reference directory of that run. A line, in a file F<pel_test.pl>,
of type: 
 $pel_run_options = "..." leads to the addition
in C<[opts...]> of anything appearing between the quotes at
the right of the = sign
(only for the run associated to the directory containing 
the considered file F<pel_test.pl>).

=item .

The exit code is tested. The test failure is pronounced
if it is non zero, unless the reference directory contains
one or more files called F<expected.err*>, in which case success
might be pronounced if identical files (same name and same content)
have been produced.

=item .

The produced F<resu> file is used to check if all instances
of subclasses of C<PEL_Object> have been destroyed when the
run terminates. If not, the test failure is pronouced.
This control can by bypassed by adding in the file F<pel_test.pl> 
a line of type:
 $objects_checking = 0

=item .

If the files called F<save.gene> (for C<TIC> postprocessing) or F<save.pel>
have been produced, they are compared with the reference ones (that
must exit) with the PELICANS-based application "pelcmp" (contained
in the argument F<exe>). If they are not identical, the comparison results
are recorded in the report file (the test failure is not pronounced since
differences may be acceptable, depending on the use case).

=item .

The files other that F<save.gene> or F<save.pel> that are produced
during the run are compared with the reference ones (that must exist) 
with the unix "cmp" utility. 
If they are not the same, the test failure is pronounced.
This comparison might be avoided from some particular file of a particular
run if the associated reference directory contains a file F<pel_test.pl>
with a line of type:
 @files_to_ignore = ( "...", "..." )
where the right of the = sign is a parentized comma separated list
of names. These names correspond to files produced during the run for
which no comparison will be performed.

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
data file. It might also contain (see above) a file called F<pel_test.pl>
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

=item B<-test_directory> F<dirname>

Duplicate the hierarchy of directories containing the reference runs
in the subdirectory F<dirname> of the working directory,
and run the application in the subdirectories of F<dirname> for further
result comparison with the reference runs (default: F<PELICANS_TEST>).

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
