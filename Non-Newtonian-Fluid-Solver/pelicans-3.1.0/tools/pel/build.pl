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
use Defs;
use POSIX ;
use Arch;
my $verbosity = '' ;

##
# Command line Parsing

my $man = '' ;
my $H = '' ;
my $target = '' ;
my $machine = POSIX::uname() ; chomp $machine ;
my $make = "make" ;
my $link_mode = "cc" ;
my $makefile = '' ;

my @with=();
my @without=();


#prepend the PELBUILD to ARGV ...
unshift(@ARGV, split(/ /,$ENV{PELBUILD})) if ($ENV{PELBUILD});

my $arch_exe = "exe" ;
if( Arch::is_posix() == 0 ) {
  $arch_exe .= ".exe" ;
} 
my $result = GetOptions ('exe' => sub { $target = $arch_exe ; }
			 , 'man' => \$man
			 , 'help' => \$H
			 , 'make=s' => \$make
			 , 'archive=s' => \$target
			 , 'object=s' =>  \$target
                         , 'link_mode=s' => \$link_mode
			 , 'verbose' =>  \$verbosity
			 , 'with=s' => \@with
			 , 'without=s' => \@without
                         , 'makefile=s' => \$makefile
			) ;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(-verbose => 1 ) if $H  ;
if( scalar @ARGV != 1 || !$result || !$target ) {
  pod2usage( -verbose => 0 ) ;  
}
@with = split(/,/,join(',',@with));
@without = split(/,/,join(',',@without));

my $lib_path = shift @ARGV ;

if( $makefile eq "" ) {
  $makefile = File::Spec->catfile( $lib_path, "Makefile" ) ;
}
if( ! -f $makefile ) {
  my $mesg = "Unknown file : $makefile\n\n" ;
  $mesg .= "Use \"pel --depend\" to build the makefile \"Makefile\" \n" ;
  $mesg .= "in the directory \"$lib_path\" prior to the\n" ;
  $mesg .= "execution of \"pel --depend\"" ;
  Defs::Error( $mesg ) ;
}

if( $verbosity ) {
  print "Compile $target on $machine \n" ;
  print "\tin $lib_path \n" ;
  print "\twith $makefile\n" ;
  print "\tenabling extra : " . join(',',@with) . "\n" if (@with);
  print "\tdisabling extra : " . join(',',@without) . "\n" if (@without);
}

Defs::Error( "Unable to find Makefile option : $makefile" ) unless -f $makefile ;
##
my @makecmd = ( $make, "-f", $makefile, "TARGET=$target", "LINK_MODE=$link_mode" ) ;
map {push @makecmd, "WITH_".uc($_)."=1"} @with;
map {push @makecmd, "WITH_".uc($_)."=0"} @without;
if( $verbosity ) {
   print "Executing command : \n   " ;
   for ( @makecmd ) { print " $_" ; }
   print "\n"
}
exec(@makecmd);
##
# POD Documentation
#
__END__

=head1 NAME

build - generation of executables, objects and libraries related to PELICANS-based applications

=head1 SYNOPSIS

pel build [-help|-man]

pel build [options...] -exe F<bindir>

pel build [options...] -object F<filename.o> F<bindir>

pel build [options...] -archive F<archname.so> F<bindir>

=head1 DESCRIPTION

The generation of executables, objects and libraries associated to
a particular PELICANS-based application can be simply performed by using
the GNU C<make> utility. The first step consists in building
a suitable makefile with the C<pel depend> utility. The second
step consists in running GNU C<make>  with suitable arguments,
a task whose responsiblity is assigned to C<pel build>.

Essentially, GNU C<make> is invoked with the
following instructions:

=over

=item 1.

read the makefile called F<Makefile>, located in the directory F<bindir>
(unless otherwise specified with the C<-makefile> options);

=item 2.

update the target determined by the mutually exclusive instructions
C<-exe,-object,-archive>;

=item 3.

use a specific set of make options that are determined by
the calling options of C<pel build>.

=back

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-e, -E, -exe>

Ask C<make> to
update the target defined by 
the executable of the PELICANS-based application.

=item B<-o, -O, -object> F<filename.o>

Ask C<make> to update the target
defined by the module object F<filename.o>
(typically, C<filename> is the name of one of
the classes making up
the considered application).

=item B<-a, -A, -archive> F<archname.ext>

Ask C<make> to update the target
defined by the archive F<archname.ext>
(the associated command depends on
the extension F<.ext>).

=item B<-link_mode> link_mode

Set the linker:
    link_mode=cc : use the c++ linker (default)
    link_mode=c  : use the c linker
    link_mode=f  : use the fortran linker

=item B<-make> makename

Use the command of name C<makename> as 
the make command (default: C<make>).
This option is typically used for systems
on which the GNU C<make> utility is
called C<gmake>.

=item B<-makefile> makefile

Set the makefile used for all tasks of compilation,
assembly and linking. Default is: F<bindir>/Makefile.

=item B<-with> EXTlist

Notify that the current application DOES require the
packages of C<EXTlist> so that the archives
associated to those packages
should be added to the list of files to link.
C<EXTlist> is a comma separated list denoting external APIs
that might be used by some components of PELICANS.
This option is meaningful only for the targets
of the options C<-exe> and C<-archive>.
Note that C<pel depend> defines some defaults in the
generated makefile for all possible external libraries.

=item B<-without> EXTlist

Notify that the current application DOES NOT require the
packages of C<EXTlist> so that the archives
associated to those packages
should NOT be added to the list of files to link.
C<EXTlist> is a comma separated list denoting external APIs
that might be used by some components of PELICANS.
This option is meaningful only for the targets
of the options C<-exe> and C<-archive>.
Note that C<pel depend> defines some defaults in the
generated makefile for all possible external libraries.

=back

=head1 ARGUMENT

=over

=item B< F<bindir> >

Directory containing the makefile generated by C<pel depend>,
called F<Makefile>. Any file produced by the execution
of C<pel build> (that is: when updating a target of F<Makefile>)
will be created in that directory (objects, libraries, executable,
dependency files ...). The executable will be called F<exe>.

=back

=head1 EXAMPLES

=over

=item C<pel depend -l pel1 dbg bin/dbg .>

=item C<pel build -exe bin/dbg>

The current application is made of all the header and source files
located in any subdirectory of the working directory. The
generated F<Makefile>, the executable F<exe> and all the files
created during the compiling process will be located in the
subdirectory F<bin/dbg> of the working directory.

=item C<pel build -without petsc,opengl -exe lib>

All the files generated for and during the compilation process
are located in the subdirectory lib of the working directory
(and possibly in the working directory itself).
Linking will be performed without the archives associated
to the PETSc and OpenGL libraries.

=item C<pel build -object myclass.o /usr/smith/appli/bin>

Update the target F<myclass.o> of the makefile F<Makefile>
located in the directory F</usr/smith/appli/bin>, that
is: build object file F<myclass.o> in that directory.

=back

=head1 ENVIRONMENT

It is possible to store arguments and options, overwritable by the command
line arguments, in the environment variable PELBUILD.

=cut
