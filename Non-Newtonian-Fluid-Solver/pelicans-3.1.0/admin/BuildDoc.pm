package BuildDoc;

use strict ;

#-------------------------------------------------------------
sub do {
#-------------------------------------------------------------
  my ( $name, $dir, $subject ) = @_ ;
  
  my @result ;
  $ENV{BIBINPUTS}='/home/semar/piar/Rapport/BIBdata' ;
  my $PELICANS=$ENV{PELICANSHOME} ;
  if( ! -d $PELICANS ) { 
     die 'invalid environment variable $PELICANSHOME' ; 
  }
  my $latex2html = "latex2html" ;
  if( -d "$dir/doc" && ! -f "$name.pdf" ) {
    
    $ENV{TEXINPUTS}="::$PELICANS/doc/share:$dir/doc" ;
    open OUT, ">doc.tex" ;

    # print OUT "\\def\\PELICANSHOME{$PELICANS}\n" ;
    open( IN, "$PELICANS/doc/share/top.tex" ) || die "cannot open : $!" ;
    while( <IN> ) {
       print OUT $_ ;
    }
    close(IN) ;

    print OUT "\\begin{document}\n" ;
    open( IN, "$dir/doc/$name.tex" ) || die "cannot open : $dir/doc/$name.tex ($!)" ;
    while( <IN> ) {
       print OUT $_ ;
       if( $_ =~ /AdditionalFile\{(.*)\}/ ) { 
	 push @result, $1 ; } 
    }
    close(IN) ;

    print OUT "\\end{document}\n" ;

    close OUT ;

    system( "$latex2html -init_file $PELICANS/doc/share/config_for_$subject.pl doc.tex" ) && die "\n latex2html compilation failed" ;

    system( "find . -name \"*.ps\" -exec rm -f {} \\;" ) ;
    system( "find . -name \"*.tex\" -exec rm -f {} \\;" ) ;
    system( "find . -name \"*.aux\" -exec rm -f {} \\;" ) ;
    system( "find . -name \"*.log\" -exec rm -f {} \\;" ) ;
    system( "find . -name \"*.out\" -exec rm -f {} \\;" ) ;
    system( "find . -name \"*.pl\" -exec rm -f {} \\;" ) ;
    system( "find . -name \"WARNINGS\" -exec rm -f {} \\;" ) ;
    system( "find . -name \"WARNINGS\" -exec rm -f {} \\;" ) ;
    if( -d "doc" ) {
      system( "rmdir doc" ) and die "Unable to rmdir doc\n" ;
    }
  }
  @result ;
}

warn "BuildDoc is successfully loaded!\n" ;
1 ;
