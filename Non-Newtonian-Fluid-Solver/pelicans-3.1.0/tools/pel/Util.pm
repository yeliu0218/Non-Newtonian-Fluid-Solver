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

package Util ;

use strict ;
use diagnostics ;
use File::Spec ;

BEGIN {
  require Exporter;
  
  use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
  $VERSION = 1.0;
  @ISA = qw(Exporter);
  @EXPORT = qw( get_field );
}

sub get_field{ 
  my ( $file, $pattern ) = @_ ;
  open HAND, $file or Defs::Error( "Unable to find $file" ) ;
  my $result ;
  foreach (<HAND>) { if(/$pattern/) { $result = $1 ; last ; } }
  $result ;
}

sub absolute_pathname {
  my( $dir ) = @_ ;
  my @dirs = File::Spec->splitdir( File::Spec->rel2abs( $dir ) ) ;
  my $i = 1 ;
  do {
    if( $dirs[$i] eq '..' ) {
      splice( @dirs, $i-1, 2 ) ;
      --$i ;
    } 
    elsif( $dirs[$i] eq '.' ) {
      splice( @dirs, $i, 1 ) ;
    }
    else {
      ++$i ;
    }
  } while( $i <= $#dirs ) ;
  my $result = File::Spec->join( @dirs ) ;
  return $result ;
}

sub common_root {
  my( @dir_elms ) = @_ ;
  my( @resu ) = () ;

  #---debug map{ print"dir_elms $_\n" } @dir_elms ;
  my $keep_going = 1 ;
  while( $keep_going ) {

    my $suivant = undef ;
    for( my $i = 0 ; $i <= $#dir_elms ; $i++ ) {
      my( @dir ) = File::Spec->splitdir( $dir_elms[$i] ) ;
      #---debug print "splitted dir : " ;
      #---debug for ( @dir ) { print "$_;" }
      my( $nn ) = shift( @dir ) ;
      #---debug print " ... tentative : $nn\n" ;
      $dir_elms[$i] = File::Spec->join( @dir ) ;
      if( $i == 0 ) {
        $suivant = $nn ;
      } elsif( defined( $suivant ) && ( $nn ne $suivant ) ) {
        $keep_going = 0 ;
        $suivant = undef ;
      }
      if( $#dir == 0 ) {
        $keep_going = 0 ;
      }
    }
    if( defined( $suivant ) ) {
       push @resu, $suivant ;
       #---debug map{ print" $_\n" } @resu ;
    }

  }
  return( File::Spec->join( @resu ) ) ;
}


sub win_relative_path {
  my ($file,$project) = @_ ;
  my $result = "" ;
  if( $file =~ /^\$\(PELICANSHOME\)/ ) {
    $result = $file ;
	$result =~ s/\//\\/g ;
  } else {
    my @filetab = () ;
    if( Arch::is_posix() == 1 ) {
      @filetab = split /\//, absolute_pathname( $file ) ;
    } else {
      @filetab = split /\\/, absolute_pathname( $file ) ;
    }  
    my @projecttab = () ;
    if( Arch::is_posix() == 1 ) {
      @projecttab = split /\//, absolute_pathname($project) ;
    } else {
      @projecttab = split /\\/, absolute_pathname($project) ;
    }
    shift @filetab ;
    shift @projecttab ;
    my $ok = 1 ;
    while( $ok ) {
      if( scalar( @filetab ) && scalar( @projecttab ) ) {
        if( "$filetab[0]" eq "$projecttab[0]" ) {
	      shift @filetab ;
	      shift @projecttab ;
        } else { last ; }
      } else { last ; }
    } 
    shift  @projecttab ;
    map { $result .= "..\\" ; } ( @projecttab ) ;
    my $least = pop @filetab ;
    map { $result .= "$_\\" ; } ( @filetab ) ;
    if ( $least ) {
      $result .= $least ;
    } else {
      $result = "." ;
    }
  }
  $result ;
}

1;
