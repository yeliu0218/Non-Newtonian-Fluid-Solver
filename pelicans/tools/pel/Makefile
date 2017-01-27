section=1
DOC=documentation_tools/
mandir:=../../doc/man/man$(section)/
cli_dir:=../../doc/Site/pel_CLI/
ps_target:=$(cli_dir)/pel_CLI.ps
manprefix:=$(mandir)pel-
# OBJ_PL:=$(wildcard *.pl)
OBJ_PL:= pel.pl depend.pl build.pl run.pl test.pl time.pl arch.pl predoc.pl newclass.pl cmp.pl cvgce.pl
OBJ_MAN:=$(addprefix $(manprefix), $(patsubst %.pl, %.$(section), $(OBJ_PL) ) )
OBJ_HTML:=$(addprefix $(cli_dir), $(patsubst %.pl,   %.html, $(OBJ_PL) ) )

LIB_TEST = pel2
CCC = gcc

ifdef libdir
PELOPT += --path $(libdir)
endif

all:doc

doc:manpage $(ps_target) html


manpage:$(mandir) $(OBJ_MAN)  

html:$(OBJ_HTML) $(cli_dir) $(cli_dir)index.html

$(cli_dir)%.html:%.pl
	pod2html --infile $< --outfile $@ --podpath=$(<D)
	\rm -f *~~

$(ps_dir):
	mkdir -p $@

$(ps_target):pel_CLI.dvi $(cli_dir) 
	dvips -o $@ pel_CLI

pel_CLI.dvi:pel_CLI.tex
	latex $<
	latex $<

pel_CLI.tex:$(OBJ_PL)
	perl -I ${DOC} ${DOC}pod2latex -full -modify -out $@ $^

$(manprefix)%.$(section):%.pl
	pod2man $< > $@

$(mandir):
	mkdir -p $@

$(cli_dir):
	mkdir -p $@


$(cli_dir)index.html:$(OBJ_HTML)
	perl ${DOC}build_html_CLI_index $@ $^

clean:
	\rm -fr $(mandir) $(ps_target) $(cli_dir) pel_CLI.*
