if( $?PATH ) then
   setenv PATH "$PELICANSHOME/bin:$PATH"
else
   setenv PATH "$PELICANSHOME/bin"
endif
if( $?MANPATH ) then
   setenv MANPATH "$PELICANSHOME/doc/man:$MANPATH"
else
   setenv MANPATH "$PELICANSHOME/doc/man"
endif   
setenv SIGALPATH "$PELICANSHOME/ExamplesOfApplication/visu"
