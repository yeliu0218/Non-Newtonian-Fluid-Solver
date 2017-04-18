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
# PELICANS Global definitions

package Defs ;

use strict ;
use diagnostics ;
use File::Spec ;
use Text::Wrap ;
use POSIX ;

BEGIN {
  require Exporter;
  
  use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
  $VERSION = 1.0;
  @ISA = qw(Exporter);
  @EXPORT = qw( package UnitTests script_cmd bin_tool Error string_vector );
  @EXPORT_OK = qw( $pelicanshome $pelpath );

  use vars qw( $pelicanshome $pelpath ) ;

  $pelicanshome = $ENV{PELICANSHOME} or Defs::Error( "PELICANSHOME not set" ) ;
  $pelpath = File::Spec->catfile( $pelicanshome, "tools", "pel" ) ;

}

sub search{
  my ($rep,$underRep ) = @_ ;
  my $dirstd = File::Spec->catfile( $pelicanshome, $rep ) ;
  opendir DIRHAND, $dirstd ;
  my @LIST;
  map { push @LIST, "$dirstd/$_" if -d "$dirstd/$_/$underRep" } File::Spec->no_upwards( readdir DIRHAND ) ;
  @LIST ;
}


sub package{
  my $rootdir = shift @_ ;

  my @STD = search( "ExamplesOfApplication", "src" ) ;

  my %package = ( "PEL" => [ "$pelicanshome/PELbase", "$pelicanshome/UnitTests/PELbase" ],
		  "LA"  => [ "$pelicanshome/LinearAlgebra", "$pelicanshome/UnitTests/LinearAlgebra" ], 
		  "GE"  => [ "$pelicanshome/Geometry", "$pelicanshome/UnitTests/Geometry" ],
		  "PDE" => [ "$pelicanshome/PDEsolver", "$pelicanshome/UnitTests/PDEsolver" ],
		  "FE" => [ "$pelicanshome/FrameFE", "$pelicanshome/UnitTests/FrameFE" ],
		  "RS" => [ "$pelicanshome/RefSol" ],
		  "LIB" => [ "$pelicanshome/DocApplication" ],
		  "STD" => [ @STD ],
		) ;
  %package ;
#  map {print "$_ @{$package{$_}}\n" ; } keys(%package) ;
}

sub UnitTests{
  my @res = search( "UnitTests", "tests" ) ;
  @res ;
}

## --Search for PELICANSHOME/tools/pel/<command>.pl script
sub script{
  my $tool = shift @_ ; 
  my $pelperl=File::Spec->catfile( $Defs::pelicanshome, 'tools', 'pel' ) ;
  my $abstool=File::Spec->catfile( $pelperl, "$tool.pl" ) ;
  Defs::Error( "Unknown command : $tool" ) unless -f $abstool ;
  my @res = ( "perl", "-w", "-I", "$pelperl", "$abstool" ) ;
  @res ;
}

sub bin_tool{
  my $tool = shift @_ ;
  my $machine = POSIX::uname() ; chomp $machine ;
  my $exe = 
    File::Spec->join( $Defs::pelicanshome, "tools", $tool, "lib", $machine, "exe" ) ;
  Defs::Error( "Bad tool $tool \n ($exe not found)" ) unless -x $exe ;
  $exe ;
}

sub Error{
print STDOUT "\n" ;
print STDOUT "-------------------------------------------------\n" ;
print STDOUT "|        Fatal error in \"pel\" utility          \n" ;
print STDOUT "-------------------------------------------------\n" ;
print STDOUT wrap ( "| ", "| ", "@_\n" ) ;
print STDOUT "-------------------------------------------------\n" ;
exit(1) ;
}

sub string_vector {
  my $key = shift @_ ;
  my $result = ' ' ;
  if( scalar @_ ) {
    $result = " $key=<" ;
    map { $result.= " \"$_\" " } @_ ;
    $result.="> " ;
  }
  $result ;
}

1;
