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

package main;

use Arch;
use Getopt::Long;
use Pod::Usage;

my $help    = 0 ;
my $man     = 0 ;
my $verbose = 0 ;
my $var     = undef ;
my $extra   = undef ;
my $result = GetOptions ( man     =>  \$man
			, help    =>  \$help
			, verbose =>  \$verbose
                        , 'getvariable_extra=s' => \$extra
			, 'getvariable=s' => \$var ) ;
pod2usage( -verbose => 1 ) if $help  ;
pod2usage( -verbose => 2 ) if $man;
if( scalar @ARGV > 1 || !$result ) {
  pod2usage( -verbose => 0 ) ;  
}

if ($verbose) {
  Arch::set_verbose(1);
}
if ( $ENV{PELICANSHOME} ) {
  Arch::set_mak_dir( File::Spec->join( $ENV{PELICANSHOME}, "etc" ) ) ;
}

if( defined $extra ) { $var = $extra ; }

if( ! $var ) {
  my $arch_fin = Arch::architecture_name( $ARGV[0] ) ;
  if( Arch::get_verbose() ) { print "Compiler Architecture name: " ; }
  print $arch_fin ;
  if( Arch::get_verbose() ) { print "\n\n" ; }
}
elsif( $var ) {
  if( $extra ) {
    Arch::set_prefix('extra-') ;
  }
  my @arch = Arch::query_best_arch( $ARGV[0] ) ;
  exit(1) unless defined $arch[0] ;
  my $makefile = Arch::get_makefile( $arch[1], $arch[2] ) ;
  open IN, "<$makefile" ;
  my $line;
  my $resu="";
  foreach $line (<IN>) {
    next if( $line =~ /^\s*\#/ ) ; # skip lines starting with comment char
    if ( $line =~ /^\s*$var\s*=\s*(.*)$/ ) {
      $resu = $1;
    }
  }
  print $resu;
}

exit(0);


##
# POD Documentation
#
__END__

=head1 NAME

arch - discover and name the compiler architecture.

=head1 SYNOPSIS

pel arch [-help|-man]

pel arch [-verbose] compiler

pel arch -getvariable <var> compiler

pel arch -getvariable_extra <var> compiler

=head1 DESCRIPTION  

C<pel arch> discovers the compiler architecture by selecting the
architecture-makefile and the extra-makefile that are the more closely
related to the current machine, the chosen compiler and external APIs, and
returns a string caracterizing the matching compiler architecture (its name).

C<pel arch> searches successively two files:

=over

=item *

the architecture-makefile, called:  F<xxx.mak>,  which essentially formalizes
the usage of the current compiler on the current machine;

=item *

the extra-makefile, called F<extra-xxx.mak>, which essentially describes the
linkage of the enabled external APIs with PELICANS on the current machine.

=back

In both cases, F<xxx> denotes symbolically a character sequence which matches
one of the following patterns tried out in sequence:

=over

=over

=item 1.

<hostname>-<compiler>

=item 2.

<hostname>

=item 3.

<sysname>-<release>-<compiler>

=item 4.

<sysname>-<compiler>

=item 5.

<compiler>

=item 6.

<sysname>-<release>

=item 7.

<sysname>

=back

=back

Where :

=over 12

=item hostname

is the name of the current host. It may be substitued if a file named
'arch_file.cfg' exists in the searched paths.

=item sysname

is the name of the curent operating system name given by uname(1).

=item release

is the release of the curent operating system name given by uname(1).

=back

=head2 The Searched Paths

First, C<pel arch> searches in the directory given by the environment
variable F<PELARCHDIR> (if defined), then in the F<$PELICANSHOME/etc>
directory (if it is not defined, the subdirectory
C<etc> of the current directory is searched instead).

=head2 Hostname Substitution

When a file named 'arch_file.cfg' is encountered in the searched paths,
C<pel arch> tries to substitute the current hostname by an alias name found
in this file. The first match found returns. When no match is found, the
current hostname is used.

This file is a two columns file. Comments starts with '#' (sharp).

=over

=over

=item Column 1 :

contains a perl (perlre(1)) regular expression matching hostnames.

=item Column 2 :

contains the alias for the regular expression.

=back

=back

'arch_file.cfg' file example:
 sinux1 pinux    # sinux1 vers pinux
 sinux\d+ sinux  # les autres noeuds sinux vers sinux
 pinux\d+ pinux  # les noeuds pinux vers pinux

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-getvariable> var

Return the value of the variable var
defined in the achitecture-makefile of the current
compiler architecture.

=item B<-getvariable_extra> var

Return the value of the variable var
defined in the extra-makefile of the current
compiler architecture.

=back

=head1 ARGUMENTS

=over

=item B<compiler>

Name of the compiler for which a compiler architecture
for the current hardware platform
will be searched.

=back

=head1 EXAMPLE

=over

=item C<pel arch -verbose CC>

Find the available compiler architecture for compiler CC
and return a string caracterizing it.

=item C<pel arch -getvariable DYNAMIC_LIB_EXT CC>

Find the available compiler architecture for compiler CC
and return the value of the variable DYNAMIC_LIB_EXT that it defines
in the achitecture-makefile.

=back

=head1 ENVIRONMENT

=over

=item PELICANSHOME

The PELICANS root directory.

=item PELARCHDIR

A user directory where user's architecture GNU Makefiles are stored.

=back

=cut
