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
use Getopt::Long;
use Pod::Usage;
use File::Find ;
use Defs;

##
# Command line Parsing

my $help = 0 ;
my $man  = 0 ;
my $verbose = 0 ;
my @searchpath = () ;
my $result = GetOptions ( 
     'help'     => \$help
   , 'man'      => \$man
   , 'verbose'  =>  \$verbose,
   , 'pelicans' =>  sub { push @searchpath, $Defs::pelicanshome}
   , 'path=s'   =>  \@searchpath,
   ) ;

pod2usage( -verbose => 1 ) if $help  ;
pod2usage( -verbose => 2 ) if $man   ;
if( scalar @ARGV != 2 || !$result ) {
  pod2usage( -verbose => 0 ) ;  
}
my( $model, $class ) = @ARGV ;

if( !@searchpath ) { @searchpath = ( "./"  ) ; }

my %srclist ;
if( $verbose ) {
  print "Searching for model $model in :\n" ;
  map { print "  $_\n" } @searchpath ;
}
File::Find::find( \&wanted, @searchpath ) ;
sub wanted {
  /^$model\.(cc|hh|icc)\z/s &&
  ( ! -l $_ ) &&
  ( ( $srclist{$1} && Defs::Error( "Found more than one model for $_" ) ) ||
    ( $srclist{$1} = $File::Find::name ) ) ;
}

while( my ($ext, $src) = each %srclist ) {
  my $dest = "$class.$ext" ;
  if( -f $dest ) { Defs::Error ( "$dest already exists" ) }
  print "Converting $src to $dest\n" if( $verbose ) ;
  open IN, $src  or Defs::Error ( "Unable to open $src" ) ;
  open OUT, ">$dest"  or Defs::Error ( "Unable to open $dest" ) ;
  while ( my $line = <IN> ) {
    if( $line =~ /^\s*\#ifndef/ ) { 
      print OUT "\#ifndef $class"."_HH\n" ;
    } elsif( $line =~ /^\s*\#define/ ) { 
      print OUT "\#define $class"."_HH\n" ;
    } else {
      $line =~ s/$model/$class/gi ;
      print OUT $line ;
    }
  }
} 

##
# POD Documentation
#
__END__

=head1 NAME

newclass - text files creation for a new class

=head1 SYNOPSIS

pel newclass [-help|-man]

pel newclass [-path F<searchdir>] model newclass

pel newclass [-pelicans] model newclass

=head1 DESCRIPTION

C<pel newclass> creates the text files associated to
a new class by copying those of a model class
and replacing the name of the model class by the name
of the new class.

The text files associated to the model class are assumed
to have one of the following extensions : ".hh", ".icc" or ".cc" .
They are searched in a set of paths defined by the calling
options, or in the current directory if none.

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-path> F<searchdir> 

Add directory F<searchdir> to the list of paths that
C<pel newclass> will search for files
F<model.hh>, F<model.cc> or F<model.icc> .
This option may be used any number of times.

=item B<-pelicans> 

Add all directories of the PELICANS repository
defined by the environment variable C<PELICANSHOME>
to the list of paths that C<pel newclass> will search for files
F<model.hh>, F<model.cc> or F<model.icc> .

=back

=head1 ARGUMENTS

=over

=item B<model>

Name of the model class. The files F<model.hh>, F<model.cc> and F<model.icc>
will be searched, copied in the current directory if found,
and any occurence of C<model> will be replaced by C<newclass>.

=item B<newclass>

Name of the new class. The produced files will be F<newclass.hh>,
F<newclass.icc> and F<newclass.cc>, prodived that related files
for the model class were found.

=back

=head1 EXAMPLES

=over

 pel newclass -pelicans GE_QuadratureRule_TEST MY_Class_TEST>

 pel   newclass -path ../src Model NewClass>

=back

=cut
