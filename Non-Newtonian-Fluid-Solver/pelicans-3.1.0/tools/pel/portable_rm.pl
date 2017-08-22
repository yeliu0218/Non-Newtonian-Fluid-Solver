#!/usr/bin/perl -w

use strict ;
use File::Path ;

foreach my $removed (@ARGV) 
{
  $removed =~ s=\\=\/=g ;
  if( -f $removed )
  {
    unlink( $removed ) or die "cannot delete $removed\n" ;
  }
  elsif( -d $removed )
  {
    rmtree( $removed ) or die "cannot deleta $removed\n" ;
  }
}  
