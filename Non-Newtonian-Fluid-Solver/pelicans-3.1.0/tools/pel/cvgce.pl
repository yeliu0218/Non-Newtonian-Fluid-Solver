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

#-------------------------------------------------------------------------
sub nbsave{

  my $nbrsave = 0 ;
  foreach my $name ( <save_*.gene> ) {
    if( $name=~ /^save_(\d*).gene/ ) {
      my $tmp = $1 ;
      if( $tmp > $nbrsave){ $nbrsave = $tmp }
    }
  }
  return $nbrsave ;
}

#-------------------------------------------------------------------------
sub extract_scalar_variables{
  my ($nbsaves) = @_ ;

  if( $verbose ) { print STDOUT "extraction of scalar variables...\n" }
  my( $val, $it ) ;
  my %hash ;
  for( $it=0 ; $it<=$nbsaves ; $it++ ) {
    open( HANDLE, "save_${it}.gene" ) ||
                  Defs::Error( "cannot open save_${it}.gene : $!" ) ;
    if( $verbose ) { print STDOUT "   save_$it.gene\n" }
    while( <HANDLE> ) {
      chomp( $_ ) ;
      if( $_ =~ /^\'(....)\'\s0/ ) {
	my $nomvar = $1 ;
	$val = <HANDLE> ;
	chomp( $val ) ;
	push( @{$hash{$nomvar}}, $val ) ;
      }
    }
    close( HANDLE ) ;
  }
  %hash ;
}

#-------------------------------------------------------------------------
sub last_time{
  my ($hashname) = @_ ;

  my %tabletmp = %${hashname} ;
  my $value ;
  my $oldvalue = -123456789 ;
  foreach( @{$tabletmp{'TIME'}} ) {
    $value = $_ ;
    if( $value < $oldvalue ){ return $oldvalue ; }
    $oldvalue = $value ;
  }
}

#-------------------------------------------------------------------------
sub check_consistency{
  my ($hashname) = @_ ;

  my %tabletmp = %${hashname} ;
  my $nbe = undef ;
  my $nn = undef ;
  if( $verbose ) { print STDOUT "found " }
  while( my ($key,$nm) = each %tabletmp ) {
    if( $verbose ) { print STDOUT "\'$key\' " }
    my $nb_vals = scalar( @{$nm} ) ;
    if( ! defined $nbe ) {
      $nbe = $nb_vals ;
      $nn = $key
    } elsif( $nbe != $nb_vals ) {
      my $mesg = "unconsistency detected in files save*.gene :\n   " ;
      $mesg .= "variable \'$nn\' has $nbe values whereas\n   " ;
      $mesg .= "variable \'$key\' has $nb_vals values" ;
      Defs::Error "$mesg" ;
    }
  }
  if( $verbose ) { print STDOUT "\nwith $nbe values each\n" }
}

#-------------------------------------------------------------------------
sub decompose_tablehash{
  my ($hashname) = @_ ;

  my %tabletmp = %${hashname};
  my %tablenew ;
  my ( $timestep, $mesh, $time, $value, $oldtime, $exist ) = () ;
  my $i = 0 ;
  foreach( @{$tabletmp{'TIME'}} ) {
    $time = $_ ;

    $timestep = shift(@{$tabletmp{'TIST'}}) ;
    $mesh = shift(@{$tabletmp{'XH  '}}) ;

    my $key_timestep = sprintf( "%2.6f", $time )."TIST" ;
    $exist = 0 ;
    foreach (@{$tablenew{$key_timestep}}){
      if( $_ == $timestep){ $exist = 1 }
      last if( $exist == 1 )
    }
    if( $exist == 0 ) {
      push( @{$tablenew{$key_timestep}}, $timestep ) ;
    }

    my $key_mesh = sprintf( "%2.6f", $time )."XH" ;
    $exist = 0 ;
    foreach (@{$tablenew{$key_mesh}}){
      if( $_ == $mesh){ $exist = 1;}
      last if( $exist == 1 )
    }
    if( $exist == 0 ) {
      push( @{$tablenew{$key_mesh}},$mesh) ;
    }

    foreach (keys %tabletmp){
      if( $_ ne 'TIME' && $_ ne 'TIST' && $_ ne 'XH  ' && $_ ne 'XD  ' ) {
        $value =  shift( @{$tabletmp{$_}} ) ;
        my $key = sprintf("%2.6f",$time) ;
        $key .= 'mesh'.sprintf("%2.6f",$mesh) ;
        $key .= 'timestep'.sprintf("%2.6f",$timestep) ;
        $key .= "$_" ;
        push( @{$tablenew{$key}}, $value );
      }
    }
    $oldtime = $time ;
  }
  %tablenew ;
}

#-------------------------------------------------------------------------
sub generate_output_files{
  my $hashname = shift @_ ;
  my $timekeys = sprintf( "%2.6f", shift @_ ) ;
  my $graph    = shift @_ ;

  my ( $fname, $timestep, $mesh, $key ) = () ;
  my $i =0 ;
  my @abscissetimestep = sort { $a <=> $b } @{$${hashname}{$timekeys."TIST"}} ;
  my @abscissemesh = sort { $a <=> $b } @{$${hashname}{$timekeys."XH"}} ;

  # generation of file : variable as a function of the time step
  # ------------------------------------------------------------
  $fname = "$graph"."_dt.txt" ;
  open( OUT,">$fname") || Defs::Error( "cannot open $fname : $!" ) ;
  if( $verbose ) {
     print STDOUT "   $fname : norm $graph function of dt (for all h)\n" ;
  }
  my $icol = 1 ;
  print OUT "\#-------\n" ;
  print OUT "\# column $icol : time step\n" ;
  foreach( @abscissemesh ) {
     $icol = $icol+1 ;
     print OUT "\# column $icol : $graph for h = $_\n" ;
  }
  print OUT "\#-------\n" ;
  foreach( @abscissetimestep ) {
    $timestep = sprintf("%2.6f",$_) ;
    print OUT $_ ;
    foreach( @abscissemesh ) {
       $key  = $timekeys ;
       $key .= 'mesh'.sprintf("%2.6f",$_) ;
       $key .= 'timestep'.$timestep ;
       $key .= $graph ;
       print OUT " ".${$${hashname}{$key}}[0] ;
    }
    print OUT "\n";
  }
  close( OUT ) ;

  # generation of file : variable as a function of the mesh size
  # ------------------------------------------------------------
  $fname = "$graph"."_h.txt" ;
  open( OUT,">$fname") || Defs::Error( "cannot open $fname : $!" ) ;
  if( $verbose ) {
     print STDOUT "   $fname  : norm $graph function of h  (for all dt)\n" ;
  }
  $icol = 1 ;
  print OUT "\#-------\n" ;
  print OUT "\# column $icol : mesh size\n" ;
  foreach( @abscissetimestep ) {
     $icol = $icol+1 ;
     print OUT "\# column $icol : $graph for dt = $_\n" ;
  }
  print OUT "\#-------\n" ;
  foreach( @abscissemesh ) {
    $mesh = sprintf("%2.6f",$_) ;
    print OUT $_ ;
    foreach( @abscissetimestep ) {
      $key  = $timekeys ;
      $key .= 'mesh'.$mesh ;
      $key .= 'timestep'.sprintf( "%2.6f", $_ ) ;
      $key .= $graph ;
      print OUT " ".${$${hashname}{$key}}[0] ;
    }
    print OUT "\n";
  }
  close( OUT ) ;
}

#-------------------------------------------------------------------------
#  MAIN starts here
#-------------------------------------------------------------------------

my $help = 0 ;
my $man  = 0 ;
my $result = GetOptions (
     'help'     => \$help
   , 'verbose'  => \$verbose
   , 'man'      => \$man
   ) ;

pod2usage( -verbose => 1 ) if $help  ;
pod2usage( -verbose => 2 ) if $man   ;
if( scalar @ARGV < 1 || !$result ) {
  pod2usage( -verbose => 0 ) ;
}
my( @variables ) = @ARGV ;

my %tbhash = extract_scalar_variables( nbsave() ) ;

check_consistency(\%tbhash) ;

my %finalhash = decompose_tablehash(\%tbhash) ;

my $lasttime = last_time(\%tbhash) ;

my %is_ok = ();
for (@variables) { $is_ok{$_} = 1 }
if( $verbose ) { print STDOUT "generation of output files...\n" ; }
foreach my $key (keys(%tbhash))
{
  if( $is_ok{$key} ) {
    generate_output_files( \%finalhash, $lasttime, $key, $key ) ;
  }
}
if( $verbose ) {
  print STDOUT "all data are considered at time : $lasttime\n" ;
}

##
# POD Documentation
#
__END__

=head1 NAME

cvgce - extraction of error norms from F<save_*.gene> text files

=head1 SYNOPSIS

pel cvgce [-help|-man]

pel cvgce [-v] <norm_saving_names>

=head1 DESCRIPTION

The convergence properties of numerical schemes when refining
the time step and the mesh size may be tested using the classes
F<FE_Launcher> and F<FE_ComparatorWithAnalytic>. In this context and
when the calculation
results are saved with F<PEL_TICwriter> (in the "text" mode), a sequence
of F<save_*.gene> files are generated.
Each F<save_*.gene> corresponds to a given
time step, to a given mesh size (specified on creation of F<FE_Launcher>),
and stores at all cycles the error norms that where requested
from F<FE_ComparatorWithAnalytic>.

Given such a sequence of F<save_*.gene> files, C<pel cvgce> extracts
all the values of the error norms, reorganizes these data to stress the
dependencies with the time step and the mesh size (for a given mesh size,
dependence with the time step ; for a given time step, dependence with
the mesh size), and finally produces multicolumn ascii files that may be
used for graphic postprocessing by usual tools (eg gnuplot or xmgrace).

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely (B<recommended>).

=back

=head1 ARGUMENTS

=over

=item B<norm_saving_names>

Items occuring in the entry of keyword C<norm_saving_names>
in the Hierarchical Data Structure used to create the
F<FE_ComparatorWithAnalytic> object.

=back

=head1 EXAMPLE

=over

=item C<pel cvgce -v XLDU XHDU>

Extract from all the F<save_*.gene> ascii files of the current directory
the entries C<XLDU> and C<XHDU>, and produce the following files:

F<XLDU_h.txt> : for all time steps, XLDU as a function of the mesh size

F<XLDU_dt.txt> : for all mesh sizes, XLDU as a function of the time step

F<XHDU_h.txt> : for all time steps, XHDU as a function of the mesh size

F<XHDU_dt.txt> : for all mesh sizes, XHDU as a function of the time step

The F<save_*.gene> files have been produced with an application using
a F<FE_ComparatorWithAnalytic> which has been created using a PELICANS
Hierarchical Data Structure. This latter must be such that C<"XLDU">
and C<"XLDU"> occur in its entry of keyword C<norm_saving_names>.

=back

=cut
