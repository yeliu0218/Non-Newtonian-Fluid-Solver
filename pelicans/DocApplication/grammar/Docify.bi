%{
#include <list>
#include <string>
#include <stdio.h>
#include <DOC_Argument.hh>
#include <DOC_Text.hh>
#include <DOC_Category.hh>
#include <DOC_Class.hh>
#include <DOC_Enum.hh>
#include <DOC_FriendDeclaration.hh>
#include <DOC_Function.hh>
#include <DOC_Sequence.hh>
#include <DOC_Method.hh>
#include <DOC_Struct.hh>
#include <DOC_Tools.hh>
#include <DOC_Type.hh>
#include <DOC_Typedef.hh>
#include <PEL_assertions.hh>
#include <PEL_Root.hh>

extern int yylex( void ) ;
extern void yyerror( const char * s )  ;
extern void switch_mode( void ) ;

DOC_Attribute::Protection  protection_courante = DOC_Attribute::Private ;
DOC_Category const* category_courante = 0 ;
DOC_Text * comment_classe = 0 ;
DOC_Text * comment_classe_old = 0 ;
DOC_Text * comment_courant = 0 ;
DOC_Method * methode_courante = 0 ;
%}

/* Description des symbols terminaux, generes par flexx */
%token YYDOC_COMMENT
%token YYDOC_INCLUDE
%token YYDOC_CLASS
%token YYDOC_IDENTIF
%token YYDOC_CATEGORY 
%token YYDOC_PUBLIC
%token YYDOC_PROTECTED
%token YYDOC_PRIVATE
%token YYDOC_CONST
%token YYDOC_VOID
%token YYDOC_VIRTUAL
%token YYDOC_STATIC
%token YYDOC_OPERATOR
%token YYDOC_TYPE_BASE
%token YYDOC_PEL_ASSERT
%token YYDOC_PEL_CHECK_PRE
%token YYDOC_PEL_CHECK_POST
%token YYDOC_PEL_LABEL
%token YYDOC_USING
%token YYDOC_FRIEND
%token YYDOC_MUTABLE
%token YYDOC_CHIFFRE
%token YYDOC_ENUM
%token YYDOC_STRUCT
%token YYDOC_CHAINE
%token YYDOC_FLECHE
%token YYDOC_FORALLEXISTS
%token YYDOC_TYPEDEF
%token YYDOC_CAST
%token YYDOC_POINTPOINT
%token YYDOC_NEW
%token YYDOC_EXTERN
/* Operateurs booleen */
%token YYDOC_EGALEGAL
%token YYDOC_DIFFEGAL
%token YYDOC_INFEGAL
%token YYDOC_SUPEGAL
%token YYDOC_OROR
%token YYDOC_ANDAND
%left YYDOC_POINTPOINT
%left '+' '-'
%left '*' '/' '%'
%left '|' YYDOC_OROR '&' YYDOC_ANDAND
%left '<' '>' YYDOC_INFEGAL YYDOC_SUPEGAL
%left YYDOC_EGALEGAL YYDOC_DIFFEGAL
%nonassoc UMINUS
%nonassoc COMP NOT
%nonassoc TYPE
%left YYDOC_FLECHE '.'
%start file
%%
   /**********************************************************************
    * Top components
    */
file:     liste_elements  { comment_classe = comment_courant = 0 ; }

liste_elements : 
	| liste_elements element   { protection_courante = DOC_Attribute::Private ; }

element : comment         
        | directive_include
        | class_predeclaration
        | class_declaration       { comment_classe_old =  0 ; }
	| method_implementation
	| using
        | extern_def
        | error 
          {
             yyerror( "unknown element" ) ;
          }

directive_include : YYDOC_INCLUDE { comment_classe = 0 ;
				  comment_courant = 0 ;
                                  DOC_Text const* txt =  $1->to_text() ;
				  if( DOC_Tools::read( DOC_Tools::file_include( txt->text() ) ) ) 
					YYACCEPT ; }

class_predeclaration : a__classe YYDOC_IDENTIF ';'

class_declaration : a__classe YYDOC_IDENTIF '{' liste_directives_classe '}' ';' {
	   string const& name = $2->to_text()->text() ;
	   DOC_Class::create(   name, 
			        comment_classe, 
			     	DOC_Tools::file(), 
				$2->to_text()->line_number(),	      
				$4->to_sequence() ) ;
	   comment_classe = 0 ;
	   $$ = 0 ; }
	| a__classe YYDOC_IDENTIF ':' YYDOC_PUBLIC YYDOC_IDENTIF '{' liste_directives_classe '}' ';' {
	   string const& name = $2->to_text()->text() ;
	   string const& name_mother = $5->to_text()->text() ;
	   DOC_Class * mother = DOC_Class::search( name_mother ) ;
	   if( mother==0 ) mother = 
	     DOC_Class::create( name_mother, 
				0, 
			     	"", 
				0,	      
				0,
				0 ) ;
	   DOC_Class::create( name, 
				comment_classe, 
			     	DOC_Tools::file(), 
				$2->to_text()->line_number(),	      
				$7->to_sequence(),
				mother ) ;
	   
	   comment_classe = 0 ;
	   $$ = 0 ; }

a__classe : YYDOC_CLASS           { 
	if( comment_classe==0 ) // Sinon, comment classes friend
	{
		comment_classe = comment_courant ; 
  		comment_courant=0 ; 
	}       
        else if( comment_classe_old==0 ) 
               comment_classe_old = comment_classe ;

	$$ = $1 ; }

liste_directives_classe : { $$ = DOC_Sequence::create(0) ; }
	| liste_directives_classe directives_classe   { 
	   if( $2 !=0 ) 
	   {
	      $1->to_sequence()->list()->append( $2->to_class_item() ) ;
	      $2->set_owner( $1 ) ;
              $$ = $1 ; 
           } }

	
directives_classe : comment          { $$ = 0 ; }
	| YYDOC_CATEGORY                  { category_courante = DOC_Category::create( 
                                             $1->to_text()->text() ) ; 
                                            $$ = 0 ; }
	| niveau_protection ':'          { $$ = 0 ; }
	| prototype ';'
	| attribute
	| class_friend 			 
	| enum                          
	| struct_dec  ';'                            
	| class_declaration              { $$ = 0 ; comment_classe = comment_classe_old ; }           
	| type_def
	| error                          { $$ = 0 ; }

struct_dec : YYDOC_STRUCT identificateur '{' liste_attribute '}' {
                string const& id = $2->to_text()->text() ;
		$$ = DOC_Struct::create( id, 
				  protection_courante,
				  category_courante,
				  comment_courant,
				  $4->to_sequence()->list() ) ;

	        comment_courant = 0 }


liste_attribute : { $$ = DOC_Sequence::create(PEL_Root::object()) ; }
	| liste_attribute attribute
                  { $1->to_sequence()->list()->append( $2 ) ; 
		    $2->set_owner($1) ;
                    $$=$1 ; }
	| liste_attribute prototype ';'
                  { $1->to_sequence()->list()->append( $2 ) ; 
                    $2->set_owner($1) ;
		    $$=$1 ; }

class_friend : YYDOC_FRIEND  { string const& str = $1->to_text()->text() ; 
			       DOC_FriendDeclaration * fr =  DOC_FriendDeclaration::create( 
					str, 
					protection_courante,
					category_courante,
					comment_courant ) ;
				$$ = fr ;
				comment_courant = 0 }

comment : YYDOC_COMMENT      { if( comment_courant!=0 )
					   {
						comment_courant->append( 
						   DOC_Tools::text_comment( $1->to_text()->text() ) ) ;
					   }
					   else
					   {
						comment_courant=DOC_Text::create( PEL_Root::object(),
						   DOC_Tools::text_comment( $1->to_text()->text() ) ) ;
					   }
					   $$ = 0 ; }
type_def: YYDOC_TYPEDEF type identificateur ';'  { 
		DOC_Type * typ = $2->to_type() ; 
                string const& str = $3->to_text()->text() ;
		$$ = DOC_Typedef::create( str, 
  			          typ,  
				  protection_courante,
				  category_courante,
				  comment_courant ) ;
	        comment_courant = 0 ; }

enum: YYDOC_ENUM identificateur '{' liste_identificateurs '}' ';' 
      {  string const& str = $2->to_text()->text() ;
	 $$ = DOC_Enum::create( str, 
			protection_courante,
			category_courante,
			comment_courant,
                        $4->to_sequence()->list() ) ;
         comment_courant = 0 ; }
	| YYDOC_ENUM '{' liste_identificateurs '}' ';' 
      {  string const& str = "" ;
	 $$ = DOC_Enum::create( str, 
			protection_courante,
			category_courante,
			comment_courant,
                        $3->to_sequence()->list() ) ;
         comment_courant = 0 ; }



liste_identificateurs: identificateur initialiseur_enum 
        { DOC_Sequence * lst= DOC_Sequence::create( PEL_Root::object() ) ;
          lst->list()->append( $1 ) ; 
	  $$ = lst ; }
	| liste_identificateurs ',' identificateur initialiseur_enum
        { DOC_Sequence * lst=$1->to_sequence() ;
          lst->list()->append( $3 ) ; 
	  $$ = lst ; }

initialiseur_enum:
	| '=' expression

niveau_protection : YYDOC_PUBLIC    { category_courante = 0 ; 
				    protection_courante = DOC_Attribute::Public ; 
				    $$ = 0 ; }
	| YYDOC_PRIVATE             { category_courante = 0 ; 
				    protection_courante = DOC_Attribute::Private ; 
				    $$ = 0 ; }
	| YYDOC_PROTECTED           { category_courante = 0 ; 
				    protection_courante = DOC_Attribute::Protected ; 
				    $$ = 0 ; }

prototype :  identificateur '(' liste_arguments ')' liste_modificateur { 
	   DOC_Method * meth =  DOC_Method::create( $1->to_text()->text(), 
					    DOC_Type::create( PEL_Root::object(), "" ), 
					    protection_courante,
					    category_courante,
					    comment_courant,
					    $3->to_sequence(),
					    $5->to_sequence(),
					    $1->to_text()->line_number(),
					    DOC_Tools::file() ) ; 
	   comment_courant = 0 ;
	   $$ = meth ; }
      	| type identificateur '(' liste_arguments ')' liste_modificateur { 
	   DOC_Method * meth = DOC_Method::create( $2->to_text()->text(), 
					 $1->to_type(), 
					 protection_courante,
					 category_courante,
 					 comment_courant,
					 $4->to_sequence(),
					 $6->to_sequence(),
					 $2->to_text()->line_number(),
					 DOC_Tools::file() ) ; 
	   comment_courant = 0 ;
	   $$ = meth ; }	
	| YYDOC_VIRTUAL prototype { DOC_Method * meth = $2->to_class_item()->method() ; 
					  PEL_ASSERT( meth!=0 ) ;
                                          meth->set_virtual() ; 
                                          $$=meth ; }
	| YYDOC_STATIC prototype  { DOC_Method * meth = $2->to_class_item()->method() ; 
					  PEL_ASSERT( meth!=0 ) ;
                                          meth->set_static() ; 
                                          $$=meth ; }

attribute : type YYDOC_IDENTIF defaut_attribute ';' { DOC_Type * typ = $1->to_type() ; 
                                  string const& str = $2->to_text()->text() ; 
				  DOC_Attribute * att =  DOC_Attribute::create( str, 
						     typ,  
						     protection_courante,
						     category_courante,
						     comment_courant ) ;
				  if( $3!=0 ) att->initialize( $3 ) ;
				  $$ = att ;
				  comment_courant = 0 }
	| YYDOC_STATIC attribute  { DOC_Attribute * att = $2->to_class_item()->attribute() ;
				 PEL_ASSERT( att!=0 ) ;
                                 att->set_static() ;
				 $$ = att ; }
	| YYDOC_MUTABLE attribute { DOC_Attribute * att = $2->to_class_item()->attribute() ;
				 PEL_ASSERT( att!=0 ) ;
                                 // att->mutable() ;
				 $$ = att ; }
        | YYDOC_STRUCT attribute { DOC_Attribute * att = $2->to_class_item()->attribute() ;
				 PEL_ASSERT( att!=0 ) ;
				 $$ = att ; }
defaut_attribute: { $$ = 0 ; }
	| '=' expression { $$ = $2 ; }

liste_modificateur :           { $$ = DOC_Sequence::create( PEL_Root::object() ) ; }
	| liste_modificateur YYDOC_CONST { 
   				 DOC_Sequence * lst = $1->to_sequence() ;
				 lst->list()->append( DOC_Text::create( PEL_Root::object(), "const" ) ) ;
				 $$ = lst ; }
	| liste_modificateur '=' YYDOC_CHIFFRE { 
   				 DOC_Sequence * lst = $1->to_sequence() ;
				 lst->list()->append( DOC_Text::create( PEL_Root::object(), "abstract" ) ) ;
				 $$ = lst ; }

identificateur: YYDOC_IDENTIF
	| YYDOC_OPERATOR
	| YYDOC_IDENTIF YYDOC_POINTPOINT identificateur  { string const& str1 = $1->to_text()->text() ;
					     string const& str2 = $3->to_text()->text() ; 
					     $$ = DOC_Text::create( PEL_Root::object(),
								  str1+"::"+str2,
							          $1->to_text()->line_number() ) ; } 
/*	| identificateur '[' ']' { string const& str1 = $1->to_text()->text() ;
                                   $$ = DOC_Text::create( PEL_Root::object(),
				        str1+"[]",$1->to_text()->line_number()  ) ; } */
 
type :   type_simple
	| type '&'             { DOC_Type * typ = $1->to_type() ; typ->set_reference() ; $$ = typ ; }
	| type '*'             { DOC_Type * typ = $1->to_type() ; typ->set_pointer() ; $$ = typ ; }
	| type YYDOC_CONST       { DOC_Type * typ = $1->to_type() ; typ->set_constant() ; $$ = typ ; }
        | identificateur '<' liste_type '>' { $$ = DOC_Type::create( PEL_Root::object(),
						 $1->to_text()->text()+"<"+
                                                 $3->to_text()->text()+">" )  }
	| type YYDOC_POINTPOINT identificateur { $$ = DOC_Type::create( PEL_Root::object(),
						 $1->to_type()->full_type_name()+"::"+
                                                 $3->to_text()->text() )  }

liste_type : type              { $$ = DOC_Text::create( PEL_Root::object(), $1->to_type()->full_type_name() ) ; }
	| liste_type ',' type  { $$ = DOC_Text::create( PEL_Root::object(), 
				      $1->to_text()->text()+","+$3->to_type()->full_type_name() ) ; }

type_simple : YYDOC_VOID              { $$ = DOC_Type::create( PEL_Root::object(), "void" ) ; }
	| YYDOC_TYPE_BASE        { string const& str = $1->to_text()->text() ; 
				$$ = DOC_Type::create( PEL_Root::object(), str ) ; }
	| identificateur       { string const& str = $1->to_text()->text() ; 
				$$ = DOC_Type::create( PEL_Root::object(), str ) ; }

liste_arguments : YYDOC_VOID     { $$ = DOC_Sequence::create( PEL_Root::object() ) ; }
	| argument val_defaut  { DOC_Sequence* lst = DOC_Sequence::create( PEL_Root::object() ) ; 
				 DOC_Argument * arg = $1->to_argument() ;
				 if( $2!=0 ) arg->initialize( $2->text() ) ;
                                 lst->list()->append( arg ) ; 
                                 $$ = lst ; }
	| liste_arguments ',' argument val_defaut  { DOC_Sequence* lst = $1->to_sequence() ; 
				                    DOC_Argument * arg = $3->to_argument() ;
				                    if( $4!=0 ) 
							arg->initialize( $4->text() ) ;
                                                    lst->list()->append( arg ) ; 
                                         	    $$ = lst ; }

val_defaut:                    { $$ = 0 ; }
	| '=' expression       { $$ = $2 ; }

argument : type  identificateur   { string const& str = $2->to_text()->text() ; 
                                 DOC_Type * typ = $1->to_type() ; 
                                 $$ = DOC_Argument::create( PEL_Root::object(), typ, str ) ; } 

method_implementation :  prototype    {  if( $1!=0 && $1->owner()!=0 )
					    { methode_courante = $1->to_class_item()->method() ;}
					 else
					    { methode_courante=0 ; if( $1) $1->destroy() ; } }
	| instruction ';'             { switch_mode() ; } 
        | var_decl

var_decl: type identificateur val_defaut ';'   
	| YYDOC_STATIC var_decl
	| YYDOC_CONST var_decl


expression: YYDOC_CHIFFRE            
	| '~' expression          %prec COMP { $$ = DOC_Function::create( PEL_Root::object(), "~", $2 ) ; }
	| '*' expression          %prec UMINUS { $$ = DOC_Function::create( PEL_Root::object(), "*", $2 ) ; }
	| '&' expression          %prec UMINUS { $$ = DOC_Function::create( PEL_Root::object(), "&", $2 ) ; }
	| '-' expression          %prec UMINUS { $$ = DOC_Function::create( PEL_Root::object(), "-", $2 ) ; }
	| '(' type_simple ')' expression %prec TYPE { $$ = $4 ; }
	| '(' expression ')'      { $$ = DOC_Function::create( PEL_Root::object(), "", $2 ) ; }
	| YYDOC_CHAINE
	| '+' '+' expression	
		{ $$ = DOC_Function::create( PEL_Root::object(), "++", $3 ) ; } 
	| expression '+' '+'	
		{ $$ = DOC_Function::create( PEL_Root::object(), "++", $1 ) ; }
	| expression '+' expression	
		{ $$ = DOC_Function::create( PEL_Root::object(), "+", $1, $3 ) ; }
	| expression '*' expression	
		{ $$ = DOC_Function::create( PEL_Root::object(), "*", $1, $3 ) ; }
	| expression '/' expression	
		{ $$ = DOC_Function::create( PEL_Root::object(), "/", $1, $3 ) ; }
	| expression '%' expression
		{ $$ = DOC_Function::create( PEL_Root::object(), "%", $1, $3 ) ; } 
	| expression '-' expression	
		{ $$ = DOC_Function::create( PEL_Root::object(), "-", $1, $3 ) ; }
	| expression YYDOC_EGALEGAL expression
		{ $$ = DOC_Function::create( PEL_Root::object(), "==", $1, $3 ) ; }
	| expression YYDOC_DIFFEGAL expression
		{ $$ = DOC_Function::create( PEL_Root::object(), "!=", $1, $3 ) ; }
	| expression YYDOC_OROR expression
		{ $$ = DOC_Function::create( PEL_Root::object(), "||", $1, $3 ) ; }
	| expression YYDOC_ANDAND expression	
		{ $$ = DOC_Function::create( PEL_Root::object(), "&&", $1, $3 ) ; }
	| expression '&' expression
		{ $$ = DOC_Function::create( PEL_Root::object(), "&", $1, $3 ) ; }
	| expression '|' expression
		{ $$ = DOC_Function::create( PEL_Root::object(), "|", $1, $3 ) ; } 
	| expression '<' expression
		{ $$ = DOC_Function::create( PEL_Root::object(), "<", $1, $3 ) ; } 
	| expression '>' expression
		{ $$ = DOC_Function::create( PEL_Root::object(), ">", $1, $3 ) ; } 
	| expression YYDOC_INFEGAL expression
		{ $$ = DOC_Function::create( PEL_Root::object(), "<=", $1, $3 ) ; } 
	| expression YYDOC_SUPEGAL expression
		{ $$ = DOC_Function::create( PEL_Root::object(), ">=", $1, $3 ) ; } 
	| '!' expression    %prec NOT      { $$ = DOC_Function::create( PEL_Root::object(), "!", $2 ) ; } 
	| forall_exp
	| identifiant
	| YYDOC_NEW type '(' liste_args ')' { $$ = 0 ; } 

cast : YYDOC_CAST '<' type '>' '(' expression ')' { $$ = DOC_Function::create( PEL_Root::object(), $1->to_text()->text(),
								     $3,
								     $6, false ) ; }

identifiant : identificateur	  { $$ = DOC_Function::create( PEL_Root::object(), $1->to_text()->text() ) ; }
	| identifiant YYDOC_FLECHE identifiant
                                  { $$ = DOC_Function::create( PEL_Root::object(), "->", $1, $3 ) ; }
	| identifiant '.' identifiant
                                  { $$ = DOC_Function::create( PEL_Root::object(), ".", $1, $3 ) ; }
	| identifiant '(' liste_args ')' 
	       { $$ = DOC_Function::create( PEL_Root::object(), $1->to_function(), $3->to_sequence()->list() ) ; }
	| cast

instruction : YYDOC_PEL_CHECK_PRE '(' expression ')'
				 { if( methode_courante ) 
                                       methode_courante->add_condition( $3->to_function(), DOC_Method::pre ) ; }
	| YYDOC_PEL_CHECK_POST  '(' expression ')'
			         { if( methode_courante ) 
                                       methode_courante->add_condition( $3->to_function(), DOC_Method::post ) ; }
	| YYDOC_PEL_ASSERT  '(' expression ')'        
				 { if( methode_courante )
				   {
                                       methode_courante->add_condition( $3->to_function(), DOC_Method::assert ) ; 
				    }}
	| YYDOC_PEL_LABEL '(' YYDOC_CHAINE ')'     
				 { if( methode_courante )
				   {
                                       methode_courante->set_label( $3->to_text()->text() ) ; 
				    }}


using : YYDOC_USING 
extern_def : YYDOC_EXTERN '{' liste_prototype '}' { $$=0 ; }

liste_prototype : { $$ =0 ;}
	| liste_prototype prototype ';'  { if($2) $2->destroy() ; $$=0 ; }
	| liste_prototype directive_include { $$ =0;}

liste_args :  { $$ = DOC_Sequence::create( PEL_Root::object() ) ; }
	| expression                  { DOC_Sequence * lst = DOC_Sequence::create( PEL_Root::object() ) ; 
					lst->list()->append( $1 ) ; 
					$$ = lst ; }
	| liste_args ',' expression   { DOC_Sequence * lst = $1->to_sequence() ; 
					lst->list()->append( $3 ) ; 
					$$ = lst ; }
	
forall_exp : YYDOC_FORALLEXISTS '(' for_desc ',' expression ')' { $$ = DOC_Function::create( PEL_Root::object(), $1->to_text()->text(),
							      $3, $5, false ) ; }

for_desc : '(' for_init ';' expression ';' expression_for ')' { $$ = DOC_Function::create( PEL_Root::object(), "",
							      $2, $4, $6 ) ; }

for_init : 		 { $$ = DOC_Sequence::create( PEL_Root::object() ) ; }
	| for_init_item 
	| for_init ',' for_init_item {
			DOC_Sequence * lst = $1->to_sequence() ;
			lst->list()->append( $3 ) ;
			$$ = lst ; }

for_init_item :  for_init_identif 
	       | for_init_identif '=' expression {
			DOC_Sequence * lst = DOC_Sequence::create( PEL_Root::object() ) ;
			DOC_Function * fonc = DOC_Function::create( PEL_Root::object(), "=", $1, $3 ) ;
			lst->list()->append( fonc ) ;
			$$ = lst ; }

for_init_identif: identifiant | argument

expression_for: { $$ = DOC_Sequence::create( PEL_Root::object() ) ; }
	| expression {  DOC_Sequence * lst = DOC_Sequence::create( PEL_Root::object() ) ;
			lst->list()->append( $1 ) ;
			$$ = lst ; }
	| expression_for ',' expression {
			DOC_Sequence * lst = $1->to_sequence() ;
			lst->list()->append( $3 ) ;
			$$ = lst ; }
