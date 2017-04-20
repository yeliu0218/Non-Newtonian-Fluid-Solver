all:../src/DOC_Docify_Lex.cc ../src/DOC_Docify_Parser.cc

../src/DOC_Docify_Lex.cc:Docify.fl
	flex -d Docify.fl 
	mv lex.yy.c ../src/DOC_Docify_Lex.cc

../src/DOC_Docify_Parser.cc:Docify.bi
	bison -d Docify.bi -o Docify.bi.tab.c
	mv Docify.bi.tab.c ../src/DOC_Docify_Parser.cc 
	mv Docify.bi.tab.h ../include/DOC_Docify.hh

clean:
	/bin/rm -f *.output ../src/DOC_Docify_Parser.cc ../src/DOC_Docify_Lex.cc ../include/DOC_Docify.hh Docify.bi.tab.*

