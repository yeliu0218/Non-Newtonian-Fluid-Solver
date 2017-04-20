if [ -z "$PATH" ]
then
   PATH="$PELICANSHOME/bin"
else
   PATH="$PELICANSHOME/bin:$PATH"
fi
export PATH
if [ -z "$MANPATH" ]
then
   MANPATH="$PELICANSHOME/doc/man"
else
   MANPATH="$PELICANSHOME/doc/man:$MANPATH"
fi
export MANPATH
SIGALPATH="$PELICANSHOME/ExamplesOfApplication/visu"
export SIGALPATH
