#!/usr/bin/perl

use strict ;
use File::Spec ;
use File::Basename ;
use Arch ;

my $program = File::Spec->rel2abs( $0 ) ;
my @directory_list = ( dirname( $program ), ".." ) ;
my $PELICANSHOME = File::Spec->catdir( @directory_list ) ;

#-----------------------------------------------------------------------------
#    License
#-----------------------------------------------------------------------------
my $print_license = 1 ;
if( scalar @ARGV == 2 && $ARGV[1] =~ "accept" ) 
{
  $print_license = 0 ;
}
if( $print_license == 1 ) 
{
  print "**********************************************************************\n" ;
  print "  PELICANS - Installation under the condition that:\n" ;
  print "      1. the license is read,\n" ;
  print "      2. the terms of the license are accepted.\n" ;
  print "**********************************************************************\n" ;
  print "Do you want to proceed? (y/n)" ;
  chomp( my $answer = <STDIN> ) ;
  exit 1 if $answer !~ /^[Yy]/ ;

  system( "more Licence_CeCILL-C_V1-en.txt" ) == 0 or die "can't print license\n" ;
  print "**********************************************************************\n" ;
  print "  You may want to re-read the CeCILL-C license agreement:\n" ;
  print "      1. in the file  Licence_CeCILL-C_V1-en.txt (english),\n" ;
  print "      2. in the file  Licence_CeCILL-C_V1-fr.txt (french),\n" ;
  print "      3. at the URL  http://www.cecill.info (with additional infos).\n" ;
  print "**********************************************************************\n" ;
  print "Do you accept the terms of the CeCILL-C license? (y/n)" ;
  chomp( $answer = <STDIN> ) ;
  exit 1 if $answer !~ /^[Yy]/ ;
  
  print "**********************************************************************\n" ;
  print "  PELICANS is the copyrighted work (1995-2010) of:\n" ;
  print "           Institut de Radioprotection et de Surete Nucléaire (IRSN)\n" ;
  print "  This software is an application framework, with a set of integrated\n" ;
  print "  reusable components, whose purpose is to simplify the task of\n" ;
  print "  developing softwares of numerical mathematics and\n" ;
  print "  scientific computing.\n" ;
  print "**********************************************************************\n" ;
  print "Do you want to proceed? (y/n)" ;
  chomp( $answer = <STDIN> ) ;
  exit 1 if $answer !~ /^[Yy]/ ;

  print "**********************************************************************\n" ;
  print "  PELICANS is governed by the CeCILL-C license under French law and\n" ;
  print "  abiding by the rules of distribution of free software.\n" ;
  print "  You can use, modify and/or redistribute the software under the terms\n" ;
  print "  of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the\n" ;
  print "  following URL http://www.cecill.info.\n" ;
  print "**********************************************************************\n" ;
  print "Do you want to proceed? (y/n)" ;
  chomp( $answer = <STDIN> ) ;
  exit 1 if $answer !~ /^[Yy]/ ;

  print "**********************************************************************\n" ;
  print "  As a counterpart to the access to the source code and rights to copy,\n" ;
  print "  modify and redistribute granted by the license, users are provided\n" ;
  print "  only with a limited warranty and the software's author, the holder\n" ;
  print "  of the economic rights, and the successive licensors have only\n" ;
  print "  limited liability.\n" ;
  print "**********************************************************************\n" ;
  print "Do you want to proceed? (y/n)" ;
  chomp( $answer = <STDIN> ) ;
  exit 1 if $answer !~ /^[Yy]/ ;

  print "**********************************************************************\n" ;
  print "  In this respect, the user's attention is drawn to the risks\n" ;
  print "  associated with loading, using, modifying and/or developing or\n" ;
  print "  reproducing the software by the user in light of its specific status\n" ;
  print "  of free software, that may mean that it is complicated to manipulate,\n" ;
  print "  and that also therefore means that it is reserved for developers and\n" ;
  print "  experienced professionals having in-depth computer knowledge.\n" ;
  print "  Users are therefore encouraged to load and test the software's\n" ;
  print "  suitability as regards their requirements in conditions enabling\n" ;
  print "  the security of their systems and/or data to be ensured and,\n" ;
  print "  more generally, to use and operate it in the same conditions as\n" ;
  print "  regards security.\n" ;
  print "**********************************************************************\n" ;
  print "Do you want to proceed? (y/n)" ;
  chomp( $answer = <STDIN> ) ;
  exit 1 if $answer !~ /^[Yy]/ ;

  print "**********************************************************************\n" ;
  print "  The fact that you are proceeding further means that you have had\n" ;
  print "  knowledge of the CeCILL-C license and that you accept its terms.\n" ;
  print "**********************************************************************\n" ;
  print "Do you want to proceed? (y/n)" ;
  chomp( $answer = <STDIN> ) ;
  exit 1 if $answer !~ /^[Yy]/ ;
  print "**********************************************************************\n" ;
}

if( Arch::is_posix() == 1 )
{
  my $output_sh = $ARGV[0].".sh" ;
  open SH, "> $output_sh" ;
   
  print SH "PELICANSHOME=$PELICANSHOME\n" ;
  print SH "export PELICANSHOME\n" ;
  print SH ". \$PELICANSHOME/bin/setvar.sh\n" ;
  
  close SH ;
 
  my $output_csh = $ARGV[0].".csh" ;
  open CSH, "> $output_csh" ;
  
  print CSH "setenv PELICANSHOME $PELICANSHOME\n" ;
  print CSH "source \$PELICANSHOME/bin/setvar.csh\n" ;
  
  close CSH ;  
}
else
{
  my $output_bat = $ARGV[0].".bat" ;
  open BAT, "> $output_bat" ;
  
  print BAT "\@set PELICANSHOME=$PELICANSHOME\n" ;
  print BAT "\@call %PELICANSHOME%\\bin\\setvar.bat\n" ;
  
  close BAT ;
}
