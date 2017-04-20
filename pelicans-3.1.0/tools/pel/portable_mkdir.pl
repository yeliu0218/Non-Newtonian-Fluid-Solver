#!/usr/bin/perl -w

use strict ;

my $directory = $ARGV[0] ;
#print "try to create directory $directory\n" ;

if( -d $directory )
{
#  print "directory exist, nothing to do\n" ;
}
else
{
  my $permission = 755 ;
  $directory =~ s=\\=\/=g ;
  my @directory_list = split( "/", $directory ) ;
  my $current_tree ;
  foreach my $tree (@directory_list)
  {
    $current_tree .= $tree."/" ;
	if( !( -d $current_tree ) )
	{
	  mkdir( $current_tree, $permission ) or die "cannot create directory $current_tree\n" ;   
	}
#	print "current tree: $current_tree\n" ;
  }
#   print "directory sucessfully created\n" ;
}
