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

package main ;

use diagnostics ;
use strict ;
use Pod::Usage;
use File::Spec ;
use Defs ;

## Command recovering

pod2usage( -verbose => 1 ) unless $ARGV[0] ;
if ( $ARGV[0] eq "-h" || $ARGV[0] eq "-help" ) {
  pod2usage( -verbose => 1 ) ;
  shift @ARGV ;
}
if ( $ARGV[0] eq "-m" || $ARGV[0] eq "-man" ) {
  pod2usage( -verbose => 2 ) ;
  shift @ARGV ;
}
my $verbose = 0 ;
if ( $ARGV[0] eq "-v" || $ARGV[0] eq "-verbose" ) {
  $verbose=1;
  shift @ARGV;
}

my $cmd = shift @ARGV ;
pod2usage( -verbose => 0 ) if !$cmd ;

## --Run script

my @script_cmd = Defs::script( "$cmd" ) ;
print "Running " . "@script_cmd" . " ". "@ARGV" ."\n" if $verbose;
#exec ( @script_cmd, @ARGV ) or Defs::Error "Unable to execute $cmd" ;
# Peut-être... à voir avec le boss...
my $exit_statut = system ( @script_cmd, @ARGV ) ;
exit $exit_statut

##
# POD Documentation
#
__END__

=head1 NAME

pel - management of PELICANS-based applications

=head1 SYNOPSIS

pel [-help|-man]

pel [options...] command [command-options-and-arguments]

command is either :

  depend         generation of makefiles
  build          generation of executables, objects and libraries
  run            execution of an application
  test           comparison between runs (for regression testing)
  time           management of a CPU time database
  arch           display of the available compiler architecture
  predoc         possible preparation before using "peldoc"
  newclass       text files creation for a new class
  cmp            requests forward to the PELICANS-based application "pelcmp"
  cvgce          extraction of error norms from "save_*.gene" files


[command-options-and-arguments] : depends on the chosen command.

=head1 DESCRIPTION  

C<pel> is a collection of utilities devoted to the management
of PELICANS-based applications.

=head1 OPTIONS

The following options are those of C<pel> itself. See the
manpages of the various commands for their specific options.

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=back

=cut
