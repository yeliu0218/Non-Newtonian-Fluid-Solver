#!/usr/local/bin/perl 

use File::Basename;
use Time::localtime;
use Time::Local;
use strict ;

my $sep = ";" ;

package Db;

sub new {
    my $class = shift;
    my $self = {};
    bless $self, $class;
    $self->{list} = () ;
    $self ;
}

sub get_enreg {
  my ( $self, $comp, $machine, $lib, $data, $t ) = @_ ;
  my $enreg = 0 ;
  my $item = {} ;
  foreach $item (@{$self->{"list"}}) {
    if( $item->{"compiler"} eq $comp &&
        $item->{"machine"} eq $machine &&
	$item->{"lib"} eq $lib &&
	$item->{"data"} eq $data ) {
      $enreg = $item ;
      last ; 
    }
  }
  $enreg ;
}

sub add_enreg {
  my ( $self, $comp, $machine, $lib, $data, $t ) = @_ ;
  my $enreg = $self->get_enreg( $comp, $machine, $lib, $data, $t ) ;
  if( $enreg ) {
    $enreg->{"time"}=$t ;
  } 
  else {
    $enreg = {} ;
    $enreg->{"compiler"} = $comp ;
    $enreg->{"machine"} = $machine ;
    $enreg->{"lib"} = $lib ;
    $enreg->{"data"} = $data ;
    $enreg->{"time"} = $t;
    push @{$self->{list}}, $enreg ;
  }  
  $enreg ;
}

sub build {
  my ($self) = @_;
  my @results=split( "\n", `find ./ -name resu -print` ) ;
  my $pel = $ENV{ "PELICANSHOME" } ;
  while( scalar @results ) {
    my $resu = shift @results ;
    my $lib ;
    my $data ;
    open HAND, "<$resu" ;
    my $line ;
    my $machine ;
    my $t ;
    my $comp ;
    foreach $line ( <HAND> ) {
      if( $line =~ /executing on (.*)$/ ) {
	$machine = $1 ; 
        chomp $machine ; }
      elsif( $line =~ /optimization level : (opt[0..2])/ ) {
	$lib = $1 ; }
      elsif( $line =~ /Data file is read from file : $pel\/(.*)$/ ) {
	$data = $1 ; }
      elsif( $line =~ /with compiler : (.*)$/ ) {
	$comp = $1 ; }
      elsif( $line =~ /^user[ \t]*(.*)$/ ) {
	$t = $1 ;
        chomp $t ; }
    }
    if( $comp eq "" || $machine eq "" || $lib eq "" || $data eq "" || $t eq "" ) {
      print "I cant read $resu \n" ;
      }
    else {
      if( $t =~ /(\d*)m(.*)s/ ) {
	  $t = 60.0*$1+$2 ;
	}
      $self->add_enreg( $comp, $machine, $lib, $data, $t ) ;
    }
    close HAND ;
  }
  $self->build_sorted_list ;
}

sub build_sorted_list {
  my ($self) = @_;
  @{$self->{list}} = sort { $a->{"compiler"}.$a->{"machine"}.$a->{"lib"} cmp $b->{"compiler"}.$b->{"machine"}.$b->{"lib"} 
                            || $b->{"data"} cmp $a->{data} } @{$self->{list}}  ;
}

sub echo {
  my ($self) = @_;
  my $k ;
  my $ref ;
  my $test ;
  my $time ;
  my $ma ;
  my $l ;
  my $d ;
  my $newer ;
  my $diff ;
format STDOUT_TOP =

 Machine ^|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
 $ma
~^|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 $ma
 Library ^<<<<<<<<<< 
 $l
                     Test name                            Time(s)  new     %
-------------------------------------------------------------------------------
.

format STDOUT =
 ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ^|||||| @|||||| @||||
 $d,                                                   $time,    $newer,  $diff
~^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
 $d
.
  my $prem = 1 ;
  foreach $k (@{$self->{"list"}}) {
    $ma=$k->{"compiler"}.$k->{"machine"} ;
    $l=$k->{"lib"} ;
    my @dpath=split( "/", $k->{"data"} ) ;
    shift @dpath ; pop @dpath ;
    my $name = shift @dpath ;
    shift @dpath ;
    unshift @dpath, $name ;
    $d = join( "/", @dpath ) ;
    $time=$k->{"time"} ;
    $newer=$k->{"newer_time"} ;
    $diff=$k->{"diff"} ;
    if( $ma ne $ref ) {
      $ref = $ma ;
      if( !$prem ) { write STDOUT_TOP ; }
      $prem = 0 ;
    }
    write STDOUT ;
  }
}

sub save {
  my ($self,$filename) = @_;
  open HAND, ">$filename" or die "Cant open file $filename" ;
  print HAND "compiler;machine;library;Test;Time (s)\n" ;
  my $k ;
  foreach $k (@{$self->{"list"}}) {
    print HAND $k->{"compiler"}.$sep.$k->{"machine"}.$sep.$k->{"lib"}.$sep.$k->{"data"}.$sep.$k->{"time"}."\n";
  }
  close HAND ;
}

sub restore {
  my ($self,$filename) = @_;
  open HAND, "<$filename" or die "Cant open file $filename" ;
  my $line = <HAND> ;
  foreach $line (<HAND>) {
    my ( $comp, $ma, $l, $d, $t ) = split( "$sep", $line ) ;
    chomp ( $comp, $ma, $l, $d, $t ) ;
    if( $ma eq "" || $l eq "" || $d eq ""|| $t eq ""  ) {
      print "Last record read in line : $line \n$comp \n$ma \n$l \n$d \n$t\n" ;
      die "Bad data file $filename\n" ;
    }
    $self->add_enreg( $comp, $ma, $l, $d, $t ) ;
  }
  close HAND ;
  $self->build_sorted_list ;
}

sub restore_old {
  my ($self,$filename) = @_;
  open HAND, "<$filename" or die "Cant open file $filename" ;
  my $line = <HAND> ;
  foreach $line (<HAND>) {
    my ( $compma, $l, $d, $t ) = split( "$sep", $line ) ;
    chomp ( $compma, $l, $d, $t ) ;
    my ( $comp, $ma ) = split( " on ", $compma ) ;
    if( $ma eq "" || $l eq "" || $d eq ""|| $t eq ""  ) {
      die "Bad data file $filename\n" ;
    }
    $self->add_enreg( $comp, $ma, $l, $d, $t ) ;
  }
  close HAND ;
  $self->build_sorted_list ;
}

sub compare {
  my ($self,$other) = @_;
  my $k = {} ;
  my $older=0.0 ;
  my $newer=0.0 ;
  my $res_db = new Db() ;

  foreach $k (@{$self->{"list"}} ) {
    my $comp = $k->{"compiler"} ;
    my $ma = $k->{"machine"} ;
    my $lib = $k->{"lib"} ;
    my $data = $k->{"data"} ;
    my $time = $k->{"time"} ;
    my $other_enreg = $other->get_enreg( $comp, $ma, $lib, $data, $time ) ;
    if( $other_enreg ) {
      $older += $time ;
      my $n = $other_enreg->{"time"} ;
      $newer += $n ;
      my $diff = ( ( $n-$time)/$time )*100.0 ; 
      $k->{"diff"} = $diff ;
      $k->{"newer_time"} = $n ;
      push @{$res_db->{"list"}}, $k ;
    }
  }
  my $percent = 0.0 ;
  if( $older > 0.0 ) { 
    $percent = ($newer-$older)*100.0/$older ;
    $res_db->echo() ;
    print "New time versus older : $newer / $older ( $percent % )\n" ;
  }
  if( $percent > 1 ) {
    exit 1 ;
  } 
  else {
    exit 0 ;
  }
}

package main ;

if( scalar @ARGV < 2 ) {
  print "$0 [-build_and_add filename] [-compare older newer] [-print filename ] [-merge older_filename added_filename]\n"
}
my $cmd = $ARGV[0]  ;
my $file = $ARGV[1] ; 

if( $cmd eq "-build_and_add" ) {
  my $db = new Db() ;
  if( -f $file ) {
    $db->restore($file) ; }
  $db->build() ;
  $db->save($file) ; 
} elsif(  $cmd eq "-print" ) {
  my $db = new Db() ;  
  $db->restore($file) ;
  $db->echo() ;
} elsif(  $cmd eq "-compare" ) {    
  my $db = new Db() ;
  $db->restore($file) ;
  my $other_db = new Db() ;
  $other_db->restore($ARGV[2]) ;
  $db->compare($other_db) ;
}elsif(  $cmd eq "-merge" ) {    
  my $db = new Db() ;
  $db->restore($file) ;
  $db->restore($ARGV[2]) ;
  $db->save($file) ; 
}elsif(  $cmd eq "-translate" ) {    
  my $db = new Db() ;
  $db->restore_old($file) ;
  $db->save($file) ; 
}
