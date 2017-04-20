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

#!/usr/local/bin/perl 

use Defs;
use File::Basename;
use File::Find;
use Time::localtime;
use Time::Local;
use strict ;
use Getopt::Long;
use Pod::Usage;
use Cwd;



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
  File::Find::find( sub { /^resu$/ && $self->read_resu( $_ ) }, "." ) ;
  $self->build_sorted_list ;
}


sub read_resu {
    my ( $self, $resu ) = @_ ;
    my $pel = $Defs::pelicanshome ;
    my $lib ;
    my $data ;
    my $dir =  Cwd::getcwd;
    if ( ! -f "$resu" ) { print "file |$resu| doesn't exist\n" ;}
    open HAND, "<$resu" or Defs::Error "Cannot open file $resu" ;
    my $line = '' ;
    my $machine = '' ;
    my $t = 'undef' ;
    my $comp = '' ;
    foreach $line ( <HAND> ) {
      if( $line =~ /Operating system: (.*)$/ ) {
	$machine = $1 ; 
        chomp $machine ; }
      elsif( $line =~ /compilation level : (opt[0..2])/ ) {
	$lib = $1 ; }
      elsif( $line =~ /\*\*\* Data file: (.*)$/ ) {
	$data = $1 ;
        $data =~ s/$pel\///g ; }
      elsif( $line =~ /compiler          : (.*)$/ ) {
	$comp = $1 ; }
      elsif( $line =~ /^user[ \t]*(.*)$/ ) {
	$t = $1 ;
        chomp $t ; }
    }
    warn( "pel::time - No compiler read in $dir/$resu" ) unless $comp ;
    warn( "pel::time - No machine read in $dir/$resu" ) unless $machine ;
    warn( "pel::time - No library read in $dir/$resu" ) unless $lib ;
    warn( "pel::time - No datafile read in $dir/$resu" ) unless $data ;
    warn( "pel::time - No user time read in $dir/$resu" ) if ($t eq 'undef');
    if( !$comp || !$machine || !$lib || !$data || ($t eq 'undef') ) {
      warn( "*** $resu can't be scanned\n" ) ;
      }
    else {
      if( $t =~ /(\d*)m(.*)s/ ) {
	  $t = 60.0*$1+$2 ;
	}
      $self->add_enreg( $comp, $machine, $lib, $data, $t ) ;
    }
    close HAND ;
}

sub build_sorted_list {
  my ($self) = @_;
  if( $self->{list} ) {
    @{$self->{list}} = sort { $a->{"machine"}.$a->{"lib"} cmp $b->{"machine"}.$b->{"lib"} 
			      || $b->{"data"} cmp $a->{data} } @{$self->{list}}  ;
  }
}

sub echo {
  my ($self) = @_;
  my $k ;
  my $ref = '' ;
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
    $ma=$k->{"machine"} ;
    $l=$k->{"lib"} ;
    my @dpath=split( "/", $k->{"data"} ) ;
    shift @dpath ; pop @dpath ;
    my $name = shift @dpath ;
    shift @dpath ;
    unshift @dpath, $name ;
    $d = join( "/", @dpath ) ;
    $time=$k->{"time"} ;
    $newer=$k->{"newer_time"} or $newer = '---' ;
    $diff=$k->{"diff"} or $diff = '---' ;
    my $newref = "$l on $ma" ; 
    if( $newref ne $ref ) {
      $ref = $newref ;
      if( !$prem ) { write STDOUT_TOP ; }
      $prem = 0 ;
    }
    write STDOUT ;
  }
}

sub save {
  my ($self,$filename) = @_;
  open HAND, ">$filename" or Defs::Error "Cannot open file $filename" ;
  print HAND "compiler;machine;library;Test;Time (s)\n" ;
  my $k ;
  foreach $k (@{$self->{"list"}}) {
    print HAND $k->{"compiler"}.$sep.$k->{"machine"}.$sep.$k->{"lib"}.$sep.$k->{"data"}.$sep.$k->{"time"}."\n";
  }
  close HAND ;
}

sub restore {
  my ($self,$filename) = @_;
  open HAND, "<$filename" or Defs::Error "Cannot open file $filename" ;
  my $line = <HAND> ;
  foreach $line (<HAND>) {
    my ( $comp, $ma, $l, $d, $t ) = split( "$sep", $line ) ;
    chomp ( $comp, $ma, $l, $d, $t ) ;
    if( $ma eq "" || $l eq "" || $d eq ""|| $t eq ""  ) {
      print "Last record read in line : $line \n$comp \n$ma \n$l \n$d \n$t\n" ;
      Defs::Error "Bad data file $filename\n" ;
    }
    $self->add_enreg( $comp, $ma, $l, $d, $t ) ;
  }
  close HAND ;
  $self->build_sorted_list ;
}

sub compare_db {
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

sub build_and_add {
  my ( $name, $file ) = @_ ;
  my $db = new Db() ;
  print "Building $file\n" ;
  $db->restore($file) if( -r $file ) ;
  $db->build() ;
  $db->save($file) ; 
}  

sub print_db {
  my ( $name, $file ) = @_ ;
  my $db = new Db() ;
  $db->restore($file) ;
  $db->echo() ;
}  

sub compare {
  my ( $name, $file ) = @_ ;
  my $newer = shift @ARGV or pod2usage( -verbose => 0 ) ;
  my $db = new Db() ;
  $db->restore($file) ;
  my $other_db = new Db() ;
  $other_db->restore($newer) ;
  $db->compare_db($other_db) ;
} 
 
sub merge {
  my ( $name, $file ) = @_ ;
  my $newer = shift @ARGV or pod2usage( -verbose => 0 ) ;
  my $db = new Db() ;
  $db->restore($file) ;
  $db->restore($newer) ;
  $db->save($file) ; 
}  


package main ;

my $verbosity = '' ;
my $H = '' ;
my $man = 0 ;
my $result = GetOptions ( 
			 'help' => \$H,
			 'man'  => \$man,
			 'build_and_add=s' => \&Db::build_and_add,
			 'compare=s' => \&Db::compare,
			 'print=s' => \&Db::print_db,
			 'merge=s' => \&Db::merge,
			 'verbose' =>  \$verbosity ) ;

pod2usage( -verbose => 1 ) if $H   ;
pod2usage( -verbose => 2 ) if $man ;
if( scalar @ARGV || !$result ) {
  pod2usage( -verbose => 0 ) ;  
}


##
# POD Documentation
#
__END__

=head1 NAME

time - management of a CPU time database for PELICANS-based applications

=head1 SYNOPSIS

pel time [-help|-man]

pel time [-verbose] -build_and_add F<file.db>

pel time [-verbose] -merge F<target.db> F<source.db>

pel time [-verbose] -compare F<file1.db> F<file2.db>

pel time [-verbose] -print F<file.db>

=head1 DESCRIPTION  

C<pel time> is devoted to the management of
data bases collecting the execution time of PELICANS-based applications.
This tool is used to manage execution timing database.
These databases are built from result files produced by PELICANS-based applications.

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-build_and_add> F<file.db>

Search recursively from the current directory all
files called F<resu>, which are supposed to contain the messages
associated to the run of a PELICANS-based application (collected
with the C<pel -run> utility). For each file F<resu>, the execution time
is read and added into F<file.db> with the directory name as key (entries
that already existed are replaced).

=item B<-merge> F<target.db> F<source.db>

Add all entries of F<source.db> into F<target.db>
(the entries that already exist in F<target.db> are replaced).

=item B<-compare> F<file1.db> F<file2.db>

Compare F<file1.db> and F<file2.db>two database and print 
the differences.

=item B<-print> F<file.db>

Display the content of F<file.db>.

=back

=head1 EXAMPLES

=over

=item C<pel time -build_and_add cpu.db>

Build the database F<cpu.db> from all the F<resu> files of the
subdirectories of the current directory.

=item C<pel time -print cpu.db>

Display the content of F<cpu.db>

=item C<pel time -compare cur.db ref.db>

Compare F<cur.db> and F<ref.db>

=back

=cut
