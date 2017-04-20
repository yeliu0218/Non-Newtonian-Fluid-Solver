%option noyywrap
%{
 // File : Gram.l
 // LCH 21.07.1999
 // Object : Lexical analysis description
 //  The valid items are :
 //     PEL__INTEGER , PEL__REAL ,  PEL__STRING
#include <iostream>
using std::istream ;
#include <cstdio>
#include <cstring>

#include <PEL_Lexical.hh>
#define YYSTYPE PEL_LexicalPtr

#include <PEL_Lexical.hh>
#include <PEL_Map.hh>
#include <PEL_Module.hh>
#include <PEL_assertions.hh>
#include <PEL_Root.hh>
#include <PEL_String.hh>
#include <PEL_Double.hh>
#include <PEL_Int.hh>
#include <PEL_Yacc.hh>

#include <stack>
using std::stack ;
#ifdef _WIN32

/* Windows doesn't define the isatty() function.
   As we always read from a file, we can easily
   provide the following simple implementation: */

int isatty (int fd)
{
    return 0;
}
#endif

extern YYSTYPE PELlval ;

/* Forward declarations */
int Identifiers(char * chain) ;
extern void PEL__Buffer(std::string const& chain) ;
extern void PELerror( const char * s) ;
extern PEL_Module* YY_top_module ;

/* Change lex input method to read from C++ streams */
std::istream*my_stream = 0 ;

#define YY_INPUT(buf,result,max_size)\
PEL_ASSERT( my_stream!=0 ) ;\
if ( *my_stream )  { my_stream->read(buf,max_size);result=my_stream->gcount();} else {result=0;}

%}
%%
<<EOF>>                 { return(PEL__EOF) ;}
\/\/.*                  { // Comment description 
                          PEL__Buffer(&yytext[0]) ; 
                        }

\"[^\"]*\"              { // String description 
                          PEL__Buffer(&yytext[0]) ;
                                  // We delete extra ""
                          yytext[ strlen(yytext)-1 ]='\0' ;
                          std::string str((char*) &yytext[1]) ;
                          size_t idx ;
                          while( (idx=str.find("'"))<=str.length() ) str.replace(idx,1,"\"") ;
                          PELlval = PEL_Lexical::create( PEL_String::create( 0, str ) );
                          return(PEL__STRING) ;
                        }

[0-9]+\.[0-9]* |
[0-9]+\.[0-9]*[eE][+-]?[0-9]+      { 
                          PEL__Buffer(&yytext[0]) ;
                          PELlval = PEL_Lexical::create( PEL_Double::create( 0,
									    atof(&yytext[0]) ) );
                          return( PEL__REAL ) ;
                        }

[0-9]+                  { // Integer description
                          PEL__Buffer(&yytext[0]) ;
                          PELlval = PEL_Lexical::create( PEL_Int::create( 0, 
                                                                         atoi(&yytext[0]) ) ) ;
                          return(PEL__INTEGER) ;
                        }

[a-zA-Z_][a-zA-Z_0-9\#]* |
"<<" | 
"||" | 
"&&" | 
"<=" | 
">=" | 
"!="                     { // Identifier description 
                          PEL__Buffer(&yytext[0]) ;
                          return( Identifiers(&yytext[0]) );
                        }

[\t\n \r]               { // Ignore tabular, carriage return white space, lf
                          PEL__Buffer(&yytext[0]) ; 
                        }

[\<\>\[\]{}#+(),?\:=\-\*/$!]   { // Ponctuation
                          PEL__Buffer(&yytext[0]) ; 
                          return( yytext[0] ) ;
                        }

.                       { // Ignore extra-character after them been signified
                          PEL__Buffer(&yytext[0]) ; 
                          std::cerr << "unknown character : " << yytext[0] << std::endl ;
                          PELerror("parse error") ;
                        }
%%

  //
  // Lexical equivalences of key-words 
  //
  
#define PEL__NB_IDENTIF (PEL__LAST-PEL__ZERO-1)

  static char PEL__EQUIVALENCES[PEL__NB_IDENTIF][20] = {
  	"MODULE", "END" , "true" , "false", "include",
         "<<", "||", "&&", "if", "<=", ">=", "!=" } ;

  int Identifiers(char * chain) 
  {
     int ret = PEL__IDENTIF ;
     for( int i=0 ; i<PEL__NB_IDENTIF ; i++ )
     {
       if(strcmp(PEL__EQUIVALENCES[i],chain)==0)
       {
          ret = PEL__ZERO+1+i ;
	  PELlval = 0 ;
	  break ;
       }
     }
     if( ret == PEL__IDENTIF )
     {
        PEL_String* chainStr = PEL_String::create( 0, chain ) ;
	PELlval = PEL_Lexical::create( chainStr ) ;
     }
     return( ret ) ;
  }
  	
static stack<istream*> file_stack ;
static stack<YY_BUFFER_STATE> buffer_stack ;

void switch_to_buffer( istream * file )
{ 
   buffer_stack.push( YY_CURRENT_BUFFER ) ;
   file_stack.push( my_stream ) ;
   my_stream = file ;
   yy_switch_to_buffer( yy_create_buffer( 0, YY_BUF_SIZE ) ) ; 
}

void un_switch_to_buffer( void )
{ 	
   YY_BUFFER_STATE old_buffer = buffer_stack.top() ; buffer_stack.pop() ;
   yy_delete_buffer( YY_CURRENT_BUFFER ) ;
   yy_switch_to_buffer( old_buffer ) ;
   my_stream = file_stack.top() ; file_stack.pop() ;
}


void PEL_re_init_lexer( void )
{ 
  while(!file_stack.empty()) file_stack.pop() ;
  while(!buffer_stack.empty()) buffer_stack.pop() ;
  my_stream = 0 ;
}


