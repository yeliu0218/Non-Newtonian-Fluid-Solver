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
use Text::Wrap ;

use Util ;
use Defs ;

#-------------------------------------------------------------------------
# MAIN starts here
#-------------------------------------------------------------------------
my $help    = 0 ;
my $man     = 0 ;
my $verbose = 0 ;
my $Wno_unresolved = 0 ;
my $result = GetOptions ( 'help'     => \$help
                        , 'man'      => \$man
                        , 'verbose'  => \$verbose 
			, 'Wno_unresolved' => \$Wno_unresolved ) ;
pod2usage( -verbose => 1 ) if $help  ;
pod2usage( -verbose => 2 ) if $man   ;
if( scalar @ARGV < 3 || !$result ) {
  pod2usage( -verbose => 0 ) ;
}

my ( $descfile, $datfile, $application_name, @scanned_dirs ) = @ARGV ;

my %found_dirs = () ;
File::Find::find(
   sub { ! -l          &&
         /^(.*)\.cc$/  &&
         ( $found_dirs{ $File::Find::dir } .= " $1" ) 
       },
   @scanned_dirs ) ;

my %classes_of_dir = () ;
while ( my ( $reldir, $classes ) = each %found_dirs ) {
  my $dir = Util::absolute_pathname( $reldir ) ;
  $classes =~ s/^\s+(.*)$/$1/ ;
  $classes_of_dir{ $dir } = $classes ;
}

if( $verbose ) {
  for ( keys %classes_of_dir ) {
    print "$_ : has \"\.cc\" files\n" ;
    print wrap ( "\t", "\t", "$classes_of_dir{$_}\n" ) ;
  }
}

%found_dirs = () ;
File::Find::find( sub { ! -l && /.hh$/ &&
                        ( $found_dirs{ $File::Find::dir } = 1 ) }
                , @scanned_dirs ) ;

my @include_dirs = () ;
for ( keys %found_dirs ) {
   push @include_dirs, Util::absolute_pathname( $_ ) ;
}

if( $verbose ) {
  for ( @include_dirs ) {
    print "$_ : has \"\.hh\" files\n" ;
  }
}

open( DAT, ">$datfile" ) or Defs::Error("Can't open $datfile : $!") ;;
print DAT "MODULE PEL_Application\n" ;
print DAT "   concrete_name = \"peldoc\"\n" ;
print DAT "   ext_filter = < \".cc\" >\n" ;
print DAT "   scanned_dir = < \n";
for ( keys %classes_of_dir ) {
   print DAT "\t\"$_\"\n" ;
}
print DAT "                 > \n";
print DAT "   include_directories = < \n";
for ( @include_dirs ) {
   print DAT "\t\"$_\"\n" ;
}
print DAT "                         > \n";
print DAT "   format = \"html\"\n" ;
my @lst = File::Spec->splitdir( $descfile ) ;
print DAT "   packaging = join( this_file_dir(), \"$lst[$#lst]\" )\n" ;
print DAT "   Wno_unresolved = true\n" if( $Wno_unresolved ) ;
print DAT "END MODULE PEL_Application\n" ;
close DAT ;

my $common_root = Util::common_root( keys %classes_of_dir ) ;

print "$common_root : common root to all directories with \"\.cc\" files\n"
   if $verbose ;

my %classes_of_pack = () ;
my %printed = () ;
my %parent = () ;
while ( my ( $dir, $classes ) = each %classes_of_dir ) {
  my $pack_of_dir = undef ;
  print "$dir\n" if $verbose ;
  my @lst = File::Spec->splitdir( $dir ) ;
  my $sonpack = undef ;
  while( File::Spec->catdir( @lst ) ne $common_root ) {
    my $pack = pop @lst ;
    if( $pack eq "src" ) {
      if( File::Spec->catdir( @lst ) ne $common_root ) {
        $pack = pop @lst ;
      }
      else {
        # there is a package of name the application itself
        $pack = $application_name ;
      }
    }
    if( !defined( $pack_of_dir ) ) {
      $pack_of_dir = $pack ;
      print "   package : $pack_of_dir\n" if $verbose ;
    }
    $parent{ $pack } = "" ;
    $printed{ $pack } = 0 ;
    if( defined( $sonpack ) ) {
      $parent{ $sonpack} = $pack ;
      print "$pack parent of $sonpack\n" if $verbose ;
    }
    $sonpack = $pack ;
  }
  $classes_of_pack{ $pack_of_dir } = $classes
}

open( DESC, ">$descfile" ) or Defs::Error("Can't open $descfile : $!") ;
print DESC "$application_name\n" ;
print DESC "Package inventory :\n" ;
print DESC "Name;Comment;Package\n" ;
my $keep_going = 0 ;
do {
  while ( my ($pack, $father) = each %parent ) {
    if( ( $father eq "" ) || $printed{ $father } ) {
      my $comment = $pack ;
      print DESC "$pack\_PACK;$comment;$father;\n";
      $printed{ $pack } = 1 ;
    }
  }
  $keep_going = 0 ;
  while ( my ($pack, $done) = each %printed ) {
    if( ! $done ) { $keep_going = 1 ; last }
  }
} while( $keep_going ) ;

print DESC "\nClass inventory :\n" ;
print DESC "Name;Package;Comment\n" ;
while ( my ($package, $classes) = each %classes_of_pack ) {
   my @lst = split /\s+/, $classes ;
   for ( @lst ) {
      print DESC "$_;$package\_PACK;;\n" ;
   }
}
close( DESC ) ;




##
# POD Documentation
#
__END__

=head1 NAME

predoc - possible preparation before using "peldoc"

=head1 SYNOPSIS

pel predoc [-help|-man]

pel predoc [options...] F<descfile> F<datfile> appli F<dirs>

=head1 DESCRIPTION

C<pel predoc> will traverse a list of directory hierarchies
from which it will infer a particular
packaging and create a data file, both
preparing a subsequent use of the PELICANS-based application C<peldoc>.

=head1 ARGUMENTS

=over

=item B< F<descfile> >

Name of the description file to be produced.

=item B< F<datfile> >

Name of the C<peldoc> data file to be produced.

=item B<appli>

Name given by C<peldoc> for the application to be documented.

=item B< F<dirs> >

List of the directories from which the packaging of the
classes will be inferred.

=back

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-Wno_unresolved>

Add options to the data file of C<peldoc> so that it will inhibit
warning messages for assertions expressed with functions that are
implemented outside the current application.

=back

=head1 EXAMPLE

=over

=item C<pel predoc doc/description.txt doc/data.pel beauty .>

A description file F<description.txt> and a C<peldoc> data file
F<data.pel> will be produced  in the subdirectory
F<doc> of the current directory, for an application that will be
called "beauty", by recursively scanning all subdirectories
of the current directory.

=back

=cut
