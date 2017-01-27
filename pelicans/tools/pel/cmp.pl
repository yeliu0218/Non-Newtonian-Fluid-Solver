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
use File::Spec ;
use File::Find ;
use Pod::Usage ;
use Arch ;
use Defs ;

my $verbose = 0 ;

#----------------------------------------------------------------------------
sub pel_executable {
#----------------------------------------------------------------------------
  my $archdir = File::Spec->join( $Defs::pelicanshome, "etc" ) ;
  Arch::set_mak_dir( $archdir ) ;
  my $arch = Arch::query_best_arch( "gcc" ) ;
  my $dir = File::Spec->join( $Defs::pelicanshome,"tests","lib", $arch ) ;
  if( $verbose ) {
    print "Looking for a PELICANS executable in:\n  $dir\n" ;
  }
  my $exe = '' ;
  File::Find::find( sub { /exe/ && ($exe = $File::Find::name) }, $dir ) ;
  if( ! -x $exe ) {
    Defs::Error( "Unable to find a valid executable file" ) ;
  }
  if( $verbose ) {
    print "Found : $exe\n" ;
  }
  return $exe ;
}

#----------------------------------------------------------------------------
# MAIN start here
#----------------------------------------------------------------------------

my $help = 0 ;
my $man  = 0 ;
my $outfile = '' ;
my @modules = () ;
my @data = () ;
my $result = GetOptions( 'help'      => \$help
                       , 'man'       => \$man
                       , 'verbose'   => \$verbose
		       , 'output'    => \$outfile 
		       , 'modules=s' => \@modules
		       , 'data=s'    => \@data
		       ) ;

pod2usage( -verbose => 1 ) if $help  ;
pod2usage( -verbose => 2 ) if $man   ;
if( scalar(@ARGV) != 2 || !$result ) {
  pod2usage( -verbose => 0 ) ;  
}
my ( $file1, $file2 ) = @ARGV ;

Defs::Error( "Invalid file : $file1" ) unless ( -r $file1 ) ;
Defs::Error( "Invalid file : $file2" ) unless ( -r $file2 ) ;

# my $exe = Defs::bin_tool( "pel" ) ;
my $exe = pel_executable() ;

my @cmd = ( "$exe", "-PEL_Application", "pelcmp", 
	  "left_file=\"$file1\"", "right_file=\"$file2\"" ) ;

if( scalar @modules ) {
  my $arg=Defs::string_vector( "valid_module", @modules ) ;
  push @cmd, $arg ;
}
if( scalar @data ) {
  my $arg=Defs::string_vector( "valid_data", @data ) ;
  push @cmd, $arg ;
}

if( $outfile ) { push( @cmd, " output_file=\"$outfile\"" ) ; }
push( @cmd, "verbose=true" ) if( $verbose ) ;

if( $verbose ) {
  print "<pel -cmp> Command line : @cmd\n" ; }
exec @cmd ;

##
# POD Documentation
#
__END__

=head1 NAME

cmp - requests forward to the PELICANS-based application "pelcmp"

=head1 SYNOPSIS

pel cmp [-help|-man] 

pel cmp [options...] F<file1> F<file2>

=head1 DESCRIPTION  

The PELICANS-based application "pelcmp" is called with the given
options and arguments.

=head1 OPTIONS

=over 8

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-m, -modules> C<mod>

Limit the comparison to C<mod>.

=item B<-d, -data> C<datum>

Limit the comparison to C<datum>

=back

=head1 ARGUMENTS

=over

=item F<file1>, F<file2> 

Name of the files containing containing the 
data to be compared.

=back

=head1 EXAMPLE

=over

=item C<pel cmp file_ref.pel file.pel>

Compare the data of the files F<file_ref.pel> and F<file.pel>.

=back

=cut
