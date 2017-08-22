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
use Defs;
use Arch;
use File::Basename ;
use File::Find ;
use File::Spec ;
use Getopt::Long;
use Pod::Usage;
use Util;

#-------------------------------------------------------------------------
# Global variables
#-------------------------------------------------------------------------
my @libs = () ;
my @libpaths = () ;
my @includedir = () ;
my @dirlist = () ;
my @dirprecomp = () ;
my @objprecomp = () ;
my $vcproj = '' ;


#-----------------------------------------------------------------------------
# Directory parsing to find source and include directories
#-----------------------------------------------------------------------------
sub mkopt{
  use vars qw/*name *dir *prune/;
  *name   = *File::Find::name;
  *dir    = *File::Find::dir;
  *prune  = *File::Find::prune;

  my $result = {} ;
  my $suffx;
  my %srcdir=();
  my $suffixesh = "h|hh";
  my $suffixes = "cpp|cc|c|F|f|$suffixesh";
  my $suffixeso = ( $vcproj ? "obj" : "o" ) ;
  my $prunedirs = "doc|tests|tools|CVS|lib|grammar";

  foreach my $f (@_) {
    File::Find::find({wanted => sub {
			/^($prunedirs)\z/s && ($File::Find::prune = 1);
			if (/^.*\.($suffixes)\z/s && ! ${$srcdir{$1}}{$dir}) {
			  #      print STDOUT "dir : $dir\nrelpath : $relpath\n" ;
			  ${$srcdir{$1}}{$dir} = $dir ; }
		      }, follow=>1}, $f) if (-d $f);
  }

  foreach my $f (@_) {
    my $dirname = File::Basename::dirname($f);
    if (-f $f && $f =~ /^.*\.($suffixes)\z/s && ! ${$srcdir{$1}}{$dirname}) {
      #      print STDOUT "dir : $dir\nrelpath : $relpath\n" ;
      ${$srcdir{$1}}{$f} = $f ;
    }
  }

  foreach $suffx (sort keys %srcdir) {
    my %tab = %{$srcdir{$suffx}} ;
    if( $suffx =~ /$suffixesh/ ) {
      map { push @includedir, $_ ; } sort keys %tab ;
    } else { 
      # all files included in dirprecomp ?
      my $dir ;
      my %obj = () ;
      foreach $dir (@dirprecomp) {
	print "Searching for precompiled objects in $dir\n" ;
	File::Find::find( { wanted => sub { if(/^(.*.$suffixeso)\z/s) {$obj{$1}=$File::Find::name;}}},$dir ) }
      my @list = keys %obj ;
      my %sorted_tab = () ;
      map { $sorted_tab{$_} = $_ ; } ( sort keys %tab ) ;
      if(scalar keys %obj) {
	foreach $dir ( keys %sorted_tab) {
	  my %prefix = () ;
	  my @tmpprecomp = () ;
	  my $o ;
	  File::Find::find( { wanted => sub {
				if(/^(.*)\.$suffx\z/s) {$prefix{$1}=$1 ;} }}, $dir);
	
	  foreach $o (keys %prefix) {
	    if($obj{"$o.$suffixeso"}){
	      push @tmpprecomp,$obj{"$o.$suffixeso"} ;
	      delete $prefix{"$o"};
	    }
	  }
	  if(scalar(keys %prefix)==0) {
	    print "Using precompiled objects for source files from $dir\n" ;
	    push @objprecomp, @tmpprecomp ;
	    delete $sorted_tab{$dir} ;
	  } elsif(scalar @tmpprecomp) {
	    print "Incomplete list of precompiled objects for $dir\n" ;
	    print "Missing following objects : \n" ;
	    map { print "$_.$suffixeso\n" ; } (keys %prefix) ;
	  }
	}
      }
      my $dirs = [] ;
      map { push @{$dirs}, $_ ; } keys %sorted_tab ;

      $result->{$suffx} = $dirs ;
    }
  }
  $result;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub inline {
  my $result = "# Inlining " . $_[0] . "\n";
  open MAKE, $_[0] ;
  $result .= $_ while (<MAKE>) ;
  close MAKE;
  $result;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub include {
  "include " . $_[0] . "\n";
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub message {
  my $mess = '#' . '-' x 79 . "\n" ;
  map {$mess .= "# " .$_ . "\n";} @_;
  $mess .= '#' . '-' x 79 . "\n";
  $mess;
}


#-------------------------------------------------------------------------
#  MAIN starts here
#-------------------------------------------------------------------------

my $verbosity = '' ;
my $man = '' ;
my $H = '' ;
my $src_pel = '' ;
my $src_std = '' ;
my $profile = '' ;
my $coverage = '' ;
my $root=$Defs::pelicanshome ;
my $import=\&inline;
my $compiler = "gcc" ;  #default
my $output = '' ;
my %defs ;


#prepend the PELDEPEND to ARGV ...
unshift(@ARGV, split(/ /,$ENV{PELDEPEND})) if ($ENV{PELDEPEND});

my $result = GetOptions (man => \$man
			 , help => \$H
			 , verbose =>  \$verbosity
			 , mPELICANS =>  \$src_pel
			 , mSTD =>  \$src_std
			 , profile => \$profile
			 , coverage => \$coverage
			 , include => sub {$import=\&include;}
			 , inline => sub {$import=\&inline;}
			 , 'root=s' => \$root
			 , 'l=s' => \@libs
			 , 'path=s' => \@libpaths
			 , 'i=s' => \@includedir
                         , 'makefile=s' => \$output
                         , 'precomp=s' => \@dirprecomp
			 , 'compiler=s' =>  \$compiler
			 , 'D:s%' =>  \%defs
			 , 'vcproj=s' => \$vcproj ) ;

pod2usage( -exitstatus => 0, -verbose => 2) if $man;
pod2usage( -verbose => 1 ) if $H  ;
if( ( !$vcproj && ( scalar @ARGV < 2 ) ) ||
      !$result ) {
  pod2usage( -verbose => 0 ) ;  
}
my %pelPackages=Defs::package($root) ;

my $option = '' ;
my $outputdir ;

if( !$vcproj ) {
  $option = shift @ARGV ;
  $outputdir = shift @ARGV ;
  $output = "$outputdir/Makefile" if !$output ;
} else {
  $vcproj = File::Spec->canonpath(File::Spec->rel2abs( $vcproj )) ;
  $outputdir = dirname( $vcproj ) ;
  $output = $vcproj ;
}

##
# Source recovering
map { push @dirlist, File::Spec->canonpath(File::Spec->rel2abs($_)) ; } @ARGV ;

my @package_list ;
my $is_posix = Arch::is_posix() ;
if( $src_pel ) { 
  @package_list=( "PEL", "GE", "LA", "PDE", "FE", "RS", "LIB" ) ; 
} elsif( $is_posix == 1 && ! $vcproj ) {
  push @includedir, "$Defs::pelicanshome/include" ;
} else {
  my @packages=( "PEL", "GE", "LA", "PDE", "FE", "RS", "LIB" ) ; 
  map { push @includedir, map({$_.="/include"} @{$pelPackages{$_}}) ; } @packages ;
}
if( $src_std ) { 
  push @package_list, "STD" ; 
}

map { push @dirlist, @{$pelPackages{$_}} ; } @package_list ;

# @dirlist or Defs::Error( "Empty directory list" ) ;
pod2usage( "No source directory provided" ) if ! @dirlist ;

if( $verbosity ) {
  print "Option is $option\n" ;
  print "Compiler is $compiler\n" ;
  print "Binary directory is $outputdir\n" ;
  print "Package read : @package_list \n" ;
  print "Directory read : \n" ;
  map { print "$_\n" ; } @dirlist ;
}

if (! $vcproj) {
  write_makefile() ;
} else {
  write_vcproj() ;
}

exit ;

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub write_makefile {

##
# Compiler option level
my $level;
if ($option eq 'dbg') {
  $level = 2;
} elsif ($option eq 'optpg') {
  # profile mode <=> opt0 + profile options
  $level = 0;
  $option = 'opt0';
  $profile = 1;  
} elsif ($option eq 'optcov') {
  # test coverage mode <=> no (optimisation or debug) + coverage options
  $level = 0;
  $coverage = 1;
} elsif ($option =~ /^opt([012])$/) {
  $level = $1;
} else {
  pod2usage( "Unrecognized option $option" );
}

##
# Open file and write it
open HAND, ">$output" or Defs::Error ( "Unable to open $output" ) ;
my $mkbody = message("Makefile generated by PELICANS depend utility"
	      , "USER DEFINED PART : START") ;

Arch::set_mak_dir(File::Spec->catfile($Defs::pelicanshome, "etc"));
my @best_arch = Arch::query_best_arch($compiler);
my $arch_mkfile = Arch::get_makefile( $best_arch[1], $best_arch[2] );
my $arch_name = Arch::architecture_name( $compiler ) ;
$mkbody .= "PELICANSHOME := " . $Defs::pelicanshome."\n";
$mkbody .= "LIBDIR        = \$(PELICANSHOME)/lib/$arch_name\n";
$mkbody .= "BINDIR        = " . File::Spec->rel2abs( "$outputdir" )."/\n" ;  #\$(LIBDIR)/$option/
$mkbody .= "OPT           = \$(flags$option)\n";
$mkbody .= "CPPFLAGS      = -DLEVEL=$level\n";
map {
  $mkbody .= "CPPFLAGS     += -D$_"; 
  $mkbody .= "=$defs{$_}" if($defs{$_});
  $mkbody .= "\n" } (keys(%defs)) ;


# Add flag to tell that pelican lib is currently building
$mkbody .= "\n" . status('MAKE_PEL', $src_pel );
$mkbody .= "\n" . status('LINK_PEL', $src_pel || $profile);
$mkbody .= "\n" . status('WITH_PROFILE', $profile);
$mkbody .= "\n" . status('WITH_COVERAGE', $coverage);
$mkbody .=<<OVER;

# runtime consistency : MAKE_PEL=1 => LINK_PEL=1
ifeq (\$(MAKE_PEL),1)
  LINK_PEL = 1
endif
OVER

$mkbody .= "\nSRC = \n";
if( @dirlist ) {
  my $list = mkopt( @dirlist ) ;
  my $suffx ;
  foreach $suffx ( keys %{$list} ) {
    map {
      my $src = $_;
      s/$Defs::pelicanshome/\$(PELICANSHOME)/g ;
      $mkbody .= "SRC += \$(wildcard ".$_."/*.".$suffx.")\n" if -d $src;
      $mkbody .= "SRC += $_\n" if -f $src;
    } @{$list->{$suffx}};
  }
}
$mkbody .= "\nPRECOMP_OBJ = \n";
map { $mkbody .= "PRECOMP_OBJ += $_\n"; } @objprecomp ;
$mkbody .= "\nINC = \n";
map { $mkbody .= "INC += $_\n"; } @includedir ;

push( @libpaths, q:$(LIBDIR): );
$mkbody .= "\nLIBPATH= \n";
map { $mkbody .= "LIBPATH += $_\n"; } @libpaths ;

$mkbody .= "\nLDLIBS=\n";
my $dirmatch = File::Spec->catfile('.*','.*');
map {
  if ( Arch::is_posix() == 1 && (-d $_ || $_ =~ /$dirmatch/) ) {
    $mkbody .= "LDLIBS += $_/*.o\n";
  } else {
    $mkbody .= "LDLIBS += -l$_\n";
  }
 } @libs ;

$mkbody .= message("USER DEFINED PART : END");

$mkbody .= message("EXTRA COMPONENT PART : START");
Arch::set_prefix("extra-");
my @best_extra = Arch::query_best_arch($compiler) ;
$mkbody .= &$import(Arch::get_makefile($best_extra[1],$best_extra[2])) if ( defined $best_extra[0] )  ;
$mkbody .= message("EXTRA COMPONENT PART : END");

$mkbody .= message("ARCHITECTURE PART : START");
Defs::Error( "Unable to find Makefile : $arch_mkfile" ) unless -f $arch_mkfile ;
$mkbody .= &$import($arch_mkfile);
$mkbody .= message("ARCHITECTURE PART : END");



my $generic_targets;
my $is_posix = Arch::is_posix();
if( $is_posix == 1 ) {
  $generic_targets =  File::Spec->catfile( $Defs::pelicanshome, "etc", "generic_targets.mak" ) ;
} else {
  $generic_targets =  File::Spec->catfile( $Defs::pelicanshome, "etc", "Windows_targets.mak" ) ;
}
Defs::Error( "Unable to find Makefile : $generic_targets" ) unless -f $generic_targets ;

$mkbody .= message("GENERIC PART : START");
$mkbody .= &$import($generic_targets);
$mkbody .= message("GENERIC PART : END");

print HAND $mkbody;
close HAND;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub write_vcproj {
    ##
    # Open file and modify it
    open VCPROJ, "<$vcproj" or Defs::Error ( "Unable to open $vcproj" ) ;
    my $proj ;
    my $out = 1 ;
    my $list = mkopt( @dirlist ) ;
    while(<VCPROJ>) {
      if( $_ =~ /AdditionalIncludeDirectories/ ) {
        $proj .= "AdditionalIncludeDirectories = \"" ;
        my @il ;
        if( $src_pel || $src_std ) {
          @il = @includedir ;
	} else {
          my $suffx ;
          foreach $suffx (@includedir) {
            my $inc = $suffx ;
            $inc =~ s/\\/\//g ;
            my $pelhome = $Defs::pelicanshome ;
            $pelhome =~ s/\\/\//g ;
            $inc =~ s/$pelhome/\$(PELICANSHOME)/g ;
            push @il, $inc ;
          }
        }
	my $last = pop @il ;
	map { $proj .= Util::win_relative_path($_,$vcproj).";"; } @il ;
	$proj .= Util::win_relative_path($last,$vcproj)."\"\n" ;
      } elsif( $_ =~ "<Files>" ) {
	$out = 0 ;
	$proj .= "<Files>\n" ;
	my @files = @objprecomp ;
        my $suffx ;
        foreach $suffx ( keys %{$list} ) {
          my $dirs = $list->{$suffx} ;
          File::Find::find({wanted => sub {
                              if (/^.*\.($suffx)\z/s ) {
			      	push @files, $File::Find::name ;
                              }
                            }, follow=>0}, @{$dirs}) if(@{$dirs});
        }
        @files = sort @files ;
        map { $proj .= "<File RelativePath=\"".Util::win_relative_path($_,$vcproj)."\"> </File>\n" } ( @files ) ;

      } elsif( $_ =~ "</Files>" ) {
	$out = 1 ;
	$proj .= "</Files>\n" ;
      }
      else {
	$proj .= $_ if ( $out ) ;
      }
    }
    close VCPROJ ;
    open HAND, ">$output" or Defs::Error ( "Unable to open $output" ) ;
    print HAND $proj ;
    close HAND ;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub status {
  my ($str, $cond) = @_;
  $str . ' = ' .($cond ? 1 : 0);
}


##
# POD Documentation
#
__END__

=head1 NAME

depend - generation of makefiles for PELICANS-based applications

=head1 SYNOPSIS

=over

pel depend [-help|-man]

pel depend [options...] arguments...

pel depend [-l lib|dir] opt F<bindir> F<sources> [-precomp directory]

pel depend -vcproj <project.vcproj> F<sources> [-precomp directory]

=back

=head1 DESCRIPTION

The task of compiling a PELICANS-based application is highly
simplified by using the GNU C<make> utility. The aim
of C<pel depend> is to write a suitable makefile
that describes the relationships among
files in the considered application and provides
commands for compiling each file and linking with the appropriate
PELICANS library. Once this makefile exists, C<pel build>
can be used to build the desired executable file.

C<pel depend> has been specially designed to handle
sources located in multiple directories, and to handle
multiple compilers on the same file system. Moreover
the generated makefile includes special targets and commands
to determine or re-determine the dependencies between the files if necessary.

The object oriented methodology of PELICANS strongly relies on
built-in assertions : preconditions, postconditions, invariants and
checks which are implemented using respectively the four
macros : C<PEL_CHECK_PRE>, C<PEL_CHECK_POST>, C<PEL_CHECK_INV>
and C<PEL_CHECK>.
These assertions can be enabled/desabled at translation time and by
using command-line switches. Some options and
arguments of C<pel depend> and C<pel run> are specially devoted to
this task.

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-D> name

In the generated makefile, predefine C<name> as a macro,
with definition 1.
This option may be used any number of times.

=item B<-D> name=value

In the generated makefile, predefine C<name> as a macro,
with definition C<value>.
This option may be used any number of times.

=item B<-compiler> comp

In the generated makefile, use the
the compiler denoted by C<comp> for all
tasks of preprocessing, compilation, assembly
and linking. Default is: C<gcc>.

=item B<-makefile> makefile

Set the name of generated makefile used for all
tasks of compilation, assembly and linking.
Default is: F<bindir>/Makefile.

=item B<-l> F<lib|dir>

Add archive file F<lib> to the list of files to link.
The commands of the generated makefile will invoke
the linker with "-lF<lib>" on the command line.

When the argument is a directory, add all objects files found in the F<dir> directory
to the list of files to link.

This option may be used any number of times.

=item B<-precomp> F<dir>

Examine directory to use precompiled objects if any.

=item B<-path> F<searchdir>

Add directory F<searchdir> to the list of paths that the linker
will search for archive libraries when invoked by the 
commands of the generated makefile. 
This option may be used any number of times.

=item B<-I> F<searchdir>

Add directory F<searchdir> to the list of paths that the compiler will search
for header files when invoked by the commands of the generated
makefile. Note that any subdirectory of <sources> that contains a header
file is automatically added to that list of paths.
This option may be used any number of times.

=item B<-profile>

In the generated makefile, activate the profiling option
when invoking the compiler and linker.
This allows the profiling analysis of the application 
with tools such as gprof.

=item B<-coverage>

In the generated makefile, activate the basic block coverage analysis option
when invoking the compiler and linker.
This allows the structural analysis of the application 
with tools such as tcov.

=item B<-inline | -include>

Decide how external makefiles (eg F<extra-Linux.mak>) taken from the PELICANS
repository are accounted for in the generated makefile.
If C<-inline>, they are copied line by line  (default
behaviour). If C<-include>, the C<include> directive of gmake is used.

=item B<-mSTD>

Add the directories of the PELICANS library "ExamplesOfApplication"
to the list of paths to be searched for sources.

=item B<-mPELICANS>

Add the directories of the PELICANS repository (except
those of the library "ExamplesOfApplication") to the list
of paths to be searched for sources. This option is
reserved to the internal administration of PELICANS.

=item B<-vcproj project.vcproj>

Update Visual C++ project by updating list of include directories
and list of files.

=back

=head1 ARGUMENTS

=over

=item B<opt>

Decide the compilation level that will be set in the
generated makefile. This option influences on the one
hand the optimization used by the compiler when generating
the binary code and on the other hand the assertions that
will be evaluated. Allowed values for C<opt>
are: C<dbg>, C<opt2>, C<opt1>, C<opt0>, C<optpg> or C<optcov>.

=over

=item B<dbg>

The commands of the generated makefile will ask the compiling system
to generate a binary code prepared for debugging.

=item B<opt2>,B<opt1>,B<opt0>

The commands of the generated makefile will ask the compiling system
to use an extensive set of optimization techniques
when generating the binary code.

=item B<dbg>,B<opt2>

In the generated makefile, the commands invoking the compiler
will define the preprocessor name C<LEVEL> as C<2>. Thus,
when running the application,
the preconditions will be evaluated, and the postconditions, invariants
and checks will be possibly evaluated (depending on command-line switches,
see C<pel run>).

=item B<opt1>

In the generated makefile, the commands invoking the compiler
will define the preprocessor name C<LEVEL> as C<1>. Thus
the preconditions will be evaluated when running the application,
but the statement associated to
postconditions, invariants and checks will be removed during
preprocessing stage.

=item B<opt0>

In the generated makefile, the commands invoking the compiler
will define the preprocessor name C<LEVEL> as C<0>. Thus
any statement associated to a precondition, a postcondition, an
invariant or a check will be removed during the
preprocessing stage.

=item B<optpg>

Same as B<opt0> combined with B<-profile>. Sets the most agressive optimisation level and sets the profiling option.

=item B<optcov>

Same as B<-coverage>. Neither optimisation options nor debug information are generated by the compiling system.

=back

=item B< F<bindir> >

Any file produced (objects, libraries, executables, dependency files ...)
will be located in F<bindir>. By default, the generated makefile,
called F<Makefile>, is created in the directory F<bindir>.

=item B< F<sources> >

A list of directories and source files.
C<pel depend> will add the given source files and will add the source files found in the given directories.
Any header file or source file found in
a subdirectory of F<sources> is considered to be part of the application.
Header files are those having a F<.h> or F<.hh> extension whereas source
files are those having a F<.cpp>, F<.cc>, F<.c>, F<.F> or F<.f> extension.

=back

=head1 EXAMPLES

=over

=item C<pel depend -l pel1 dbg lib .>

The current application is made of all the header and source files
located in any subdirectory of the working directory. 
The generated F<Makefile> will be created
in the directory F<lib>. Further compilations will be performed with the
F<dbg> compilation level, and linking will be performed with the library
F<libpel1.so>.

=item C<pel depend -l pel0 -I hea -compiler gcc opt1 bin src>

The current application is made of all the header and source files
located in any subdirectory of F<src> . The generated F<Makefile> will be
created in the directory F<bin>. Further compilations will be performed
by C<gcc> with the C<opt1> compilation level, and linking will be performed
with the library F<libpel0.so>. The current application probably uses header
files that are not in the directory F<src> and that are not PELICANS header
files since the options C<-I hea> is used.

=item C<pel depend opt1 bin src pelsrc /home/users/algo.cc -precomp lib>

Compile the single source file /home/users/algo.cc and the source directories src and pelsrc using precompiled objects in lib.

=back

=head1 ENVIRONMENT

It is possible to store arguments and options, overwritable by the command
line arguments, in the environment variable PELDEPEND.

=cut
