all:../src/PEL_Yacc.cc ../src/PEL_Lex.cc 

../src/PEL_Lex.cc:Gram.l
	flex -PPEL -d Gram.l
	mv lex.PEL.c ../src/PEL_Lex.cc

../src/PEL_Yacc.cc:Gram.y
	bison -p PEL -d Gram.y -o y.tab.c
	mv y.tab.c ../src/PEL_Yacc.cc
	mv y.tab.h ../include/PEL_Yacc.hh

clean:
	/bin/rm -f lex.PEL.c y.tab.*  y.output ../src/PEL_Lex.cc ../src/PEL_Yacc.cc ../include/PEL_Yacc.hh
