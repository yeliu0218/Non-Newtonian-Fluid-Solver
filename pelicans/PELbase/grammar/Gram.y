%{
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
using std::istream ;
   
#include <PEL_assertions.hh>
#include <PEL_Context.hh>
#include <PEL_Double.hh>
#include <PEL_Map.hh>
#include <PEL_Lexical.hh>
#include <PEL_Expression.hh>
#include <PEL_Module.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_List.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_Root.hh>
#include <PEL_Data.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_IntVector.hh>
#include <PEL_BoolArray2D.hh>
#include <PEL_BoolVector.hh>
#include <PEL_StringArray2D.hh>
#include <PEL_StringVector.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>
#include <PEL_Bool.hh>
#include <PEL.hh>

#include <stringVector.hh>

/*-----------------------------------------------.
| Public functions.                              |
`-----------------------------------------------*/
bool PEL_readFile( PEL_Module * top,
                   istream* input_stream,
                   std::string const& name,
                   bool debug = false ) ;
std::string PEL_current_parsed_module_path_name( void ) ;
int PEL_current_parsed_line( void ) ;
void PEL_re_init_lexer( void ) ;


#define YYSTYPE PEL_LexicalPtr
int yylex( void ) ;
extern int PEL_flex_debug ;

#include <stack>
using std::stack ;

static stack<bool> main_stack ;
static stack<int> nb_line_stack ;
static stack<std::string> path_stack ;
static stack<std::string> name_stack ;
static stack<istream*> file_stack ;

bool endFile( void ) ;
std::string comment( void ) ;
void PELerror( const char * ) ;
bool readFile( PEL_Module * top,
               istream* input_stream,
               std::string const& name,
               bool main_read ) ;
void switch_to_buffer( istream * file ) ;
void un_switch_to_buffer( void ) ;
void substitute_for_assignment( std::string const& key,
                                PEL_Data* data ) ;

static int PEL__NbLines = 1 ;
static std::string buff ;
static std::string relativePath ;
static std::string currentFile = "" ;
static bool parsing = false ;

static PEL_Module * YY_top_module = 0 ;
static PEL_Module * dummy_module = 0 ;
static PEL_List * modules_LILO = 0 ;

static PEL_Context const* CTX_INI = 0 ;

#define CREATE_OP(op,name,lst)\
             stringVector const& exps = \
                PEL_Expression::registered_expressions() ; \
             if( !exps.has( name ) ) \
             { std::string mess = "Unknown operator \"" ; \
               mess += name ; \
               mess += "\"\n" ; \
               mess += "valid ones are :\n" ; \
               stringVector e = exps ; \
               e.sort() ; \
               for( size_t i=0 ; i<e.size() ; ++i ) \
                  mess += "   - \""+e(i)+"\"\n" ; \
               PELerror( mess.c_str()  ) ; } \
             if( !PEL_Expression::valid_arguments_of( name, lst ) ) \
             { std::string mess =  "Valid syntax for \"" ; \
               mess += name ; \
               mess += "\" operator is : \n  " ; \
               mess += PEL_Expression::usage_of( name ) ; \
               PELerror( mess.c_str() ) ; } \
             PEL_Expression * op = PEL_Expression::create( 0, name, lst, comment() ) ; \
             op->unset_external_brackets() ;

#define CREATE_SINGLE_OP(name,arg1,resu) \
             PEL_List * lst = PEL_List::create( 0 ) ; \
             lst->append( arg1->to_data() ) ; \
             arg1->change_owner(lst,arg1->to_data() ) ; \
             CREATE_OP( op, name, lst ) ; \
             lst->set_owner( op ) ; \
             resu=PEL_Lexical::create( op )

#define CREATE_BIN_OP(name,arg1,arg2,resu) \
             PEL_List * lst = PEL_List::create( 0 ) ; \
             lst->append( arg1->to_data() ) ; \
             lst->append( arg2->to_data() ) ; \
             arg1->change_owner(lst,arg1->to_data() ) ; \
             arg2->change_owner(lst,arg2->to_data() ) ; \
             CREATE_OP( op, name, lst ) ; \
             lst->set_owner( op ) ; \
             resu=PEL_Lexical::create( op )

%}
/* Grammar description */
%token PEL__IDENTIF PEL__STRING PEL__REAL PEL__INTEGER PEL__EOF
/* Grammar key-words */
%token PEL__ZERO
%token PEL__MODULE 
%token PEL__END 
%token PEL__TRUE PEL__FALSE
%token PEL_INCLUDE
%token PEL__CONCAT
%token PEL__OR
%token PEL__AND
%token PEL__IF
%token PEL__LE
%token PEL__GE
%token PEL__NEQ
%token PEL__LAST
/* Special grammar key-words */
%start data_file
%left PEL__CONCAT
%left PEL__OR PEL__AND
%left '=' PEL__NEQ
%left '<' '>' PEL__LE PEL__GE
%left '+' '-'
%left '*' '/'
%nonassoc UMINUS
%nonassoc UNOT
%%
   /**********************************************************************
    * Top components
    */
      
data_file : /* empty */
           | data_file item_data_file 

item_data_file : module_def
               | directive

directive: '#' PEL_INCLUDE path
               {
                  if( $3->to_data()->data_type() != PEL_Data::String )
                       PELerror( "Include must refer to string expression" ) ;
                  if( YY_top_module!=dummy_module )
                  {
                    if( ! $3->to_data()->value_can_be_evaluated(
                                               YY_top_module->context() ) )
                    {
                       std::string msg = "Include path can not be evaluated\n"  ;
                       msg += "Undefined variable(s):\n" ;
                       stringVector const& undef =
                          $3->to_data()->undefined_variables(
                                                  YY_top_module->context() ) ;
                       for( size_t i=0 ; i<undef.size() ; ++i )
                       {
                          msg += "   - \""+undef(i)+"\"\n" ;
                       }
                       PELerror( msg.c_str() ) ;
                    }
                    std::string file_data =
                       $3->to_data()->to_string( YY_top_module->context() ) ;
                    readFile( YY_top_module, 0, file_data, false ) ;
                  }
               }
          | PEL__EOF { if( endFile() ) YYACCEPT ; }

path : '(' something ')' { $$ = $2 ; }
     | PEL__STRING
     
   /**********************************************************************
    * Utilities specifications
    */
    
   /*
    * Free list specifications
    */
    
free_list : /* Empty list */ 
        | free_list subsitute
        | free_list assignment
        | free_list variable_def
        | free_list module_def
        | free_list directive

assignment : PEL__IDENTIF '=' something
             {
                if( YY_top_module!=dummy_module )
                {
                   if( $1->to_data()->data_type()!=PEL_Data::String )
                      PELerror( "Key must be a string" ) ;
                   std::string const& str = $1->to_data()->to_string() ;
                   if( YY_top_module->has_entry( str ) )
                   {
                     std::string mess = "Module " ;
                     mess+=YY_top_module->name() + " already contains "+
                        str + " assignment (use ==)" ;
                     PELerror( mess.c_str() ) ;
                   }
                   PEL_Data* data = $3->to_data() ;
                   $3->change_owner(YY_top_module, data) ;
                   if( YY_top_module->has_module( str ) )
                   {
                      std::string mess = "Can't create entry and module with the same name "+str+" in module " + YY_top_module->name() ;
                      PELerror( mess.c_str() ) ;
                   }
                   YY_top_module->add_entry( str, data ) ;
                }
                $$ = 0 ;
             }

variable_def : variable '=' something
             {
                if( YY_top_module!=dummy_module )
                {
                   PEL_Variable* var = dynamic_cast<PEL_Variable*>(
                      $1->to_data() ) ;
                   $1->change_owner(YY_top_module,var) ;
                   PEL_ASSERT( var!=0 ) ;
                   if( var->data_type()!=$3->to_data()->data_type() )
                   {
                      std::string mess = "Bad type for expression assigned to " ;
                      mess+=var->name() ;
                      PELerror( mess.c_str() ) ;
                   }
                   if( !YY_top_module->context()->has_variable(var) )
                   {
                      PEL_Data* data = $3->to_data() ;
                      $3->change_owner(0,data) ;
                      YY_top_module->add_variable( var, data ) ;
                   }
                   else if( !CTX_INI->has_variable(var) )
                   {
                      std::string mess = "Module " ;
                      mess+=YY_top_module->name() + " already contains "+
                            var->name() + " variable (use ==)" ;
                      PELerror( mess.c_str() ) ;
                   }
// Module initial context is not modified...
//                   else
//                   {
//                      PEL_Data* data = $3->to_data() ;
//                      $3->change_owner(0,data) ;
//                      YY_top_module->modify_variable( var, data ) ;
//                   }
                }
                $$ = 0 ;
             }
variable_def : variable '=''=' something
             {
                if( YY_top_module!=dummy_module )
                {
                   PEL_Variable* var = dynamic_cast<PEL_Variable*>(
                      $1->to_data() ) ;
                   PEL_ASSERT( var!=0 ) ;
                   if( var->data_type()!=$4->to_data()->data_type() )
                   {
                      std::string mess = "Bad type for expression assigned to " ;
                      mess+=var->name() ;
                      PELerror( mess.c_str() ) ;
                   }
                   if( !YY_top_module->context()->has_variable(var) )
                   {
                      std::string mess = var->name() + " doesn't already exist " ;
                      PELerror( mess.c_str() ) ;
                   }
                   PEL_Data* data = $4->to_data() ;
                   $4->change_owner(0,data) ;
                   YY_top_module->modify_variable( var, data ) ;
                }
                $$ = 0 ;
             }

subsitute : PEL__IDENTIF '=''=' something
             {
                if( YY_top_module!=dummy_module )
                {
                   if( $1->to_data()->data_type()!=PEL_Data::String )
                      PELerror( "Key must be a string" ) ;
                   PEL_String * str = static_cast<PEL_String*>( $1->to_data() ) ;
                   if( !YY_top_module->has_entry( str->to_string() ) )
                   {
                      std::string mess = str->to_string() + " doesn't already exist " ;
                      PELerror( mess.c_str() ) ;
                   }
                   substitute_for_assignment( str->to_string(),
                                              $4->to_data() ) ;
                }
                $$ = 0 ;
             }

something : SimpleType 
     | vector
     | array
     | function
     | variable
     | test
     | binary_operator
     | unary_operator
     |'(' something ')'
          {
             $$ = $2 ;
             if( $$->is_data() )
             {
                PEL_Expression* op = dynamic_cast<PEL_Expression*>( $$->to_data() ) ;
                if( op != 0 ) op->set_external_brackets() ;
             }
          }
   

function: PEL__IDENTIF '(' liste_args ')' 
          {
             std::string const& op_name = $1->to_data()->to_string() ;
             if( op_name=="this_file_dir" )
             {
                $$=PEL_Lexical::create(
                   PEL_String::create( 0, relativePath ) ) ;
             }
             else
             {
                CREATE_OP( op, op_name, $3->to_list() ) ; 
                $$=PEL_Lexical::create( op ) ;
                $3->change_owner(op,$3->to_list()) ;
             }
          }

variable: '$' PEL__IDENTIF 
          {
             std::string const& str = $2->to_data()->to_string() ;
             if( str.length()<2 ||
                 ( str[0]!='I' && str[0]!='D' && str[0]!='B' && str[0]!='S' )
                 || ( str[1]!='S' && str[1]!='V' && str[1]!='A' ) )
             {
                std::string msg = "\""+str+"\" is not a valid variable name\n" ;
                msg += "A valid name is \"XY_name\"\n" ;
                msg += "   where \"X\" is the scalar type of the variable :\n" ;
                msg += "       - \"I\" : integer\n" ;
                msg += "       - \"D\" : double\n" ;
                msg += "       - \"B\" : boolean\n" ;
                msg += "       - \"S\" : string\n" ;
                msg += "   and \"Y\" defined its dimension :\n" ;
                msg += "       - \"S\" : simple (only one element)\n" ;
                msg += "       - \"V\" : vector\n" ;
                msg += "       - \"A\" : array2D\n" ;
                msg += "Examples : \"DV_coordinates\", \"SS_name\", \"IA_connectivity\"\n" ;
                PELerror( msg.c_str() ) ;
             }
             PEL_Variable const* var = PEL_Variable::object(
                $2->to_data()->to_string() ) ;
             $$=PEL_Lexical::create( var->create_clone(0) ) ;
          }

test: '(' switches_list something ')' 
          {
             PEL_List * lst = $2->to_list() ;
             $3->change_owner(lst,$3->to_data() ) ;
             lst->append( $3->to_data() ) ;
             CREATE_OP( op, "(?:)",lst ) ;
             op->set_external_brackets() ;
             $2->change_owner(op,lst) ;
             $$=PEL_Lexical::create( op ) ;
          }

switches_list :  something '?' something ':'
              {
                 PEL_LABEL( "Gram.y :: switches_list1" ) ;
                 PEL_List * lst = PEL_List::create( 0 ) ;
                 $1->change_owner(lst,$1->to_data() ) ;
                 lst->append( $1->to_data() ) ;
                 $3->change_owner(lst,$3->to_data() ) ;
                 lst->append( $3->to_data() ) ;
                 $$ = PEL_Lexical::create( lst ) ;
              }
  | switches_list something '?' something ':'
              {
                 PEL_LABEL( "Gram.y :: switches_list2" ) ;
                 PEL_List * lst = $1->to_list() ;
                 $2->change_owner(lst,$2->to_data() ) ;
                 lst->append( $2->to_data() ) ;
                 $4->change_owner(lst,$4->to_data() ) ;
                 lst->append( $4->to_data() ) ;
                 $$ = $1 ;
              }

liste_args:   {  $$=PEL_Lexical::create( PEL_List::create( 0 ) ) ; }
  | something
              {
                 PEL_List * lst = PEL_List::create( 0 ) ;
                 $1->change_owner(lst,$1->to_data() ) ;
                 lst->append( $1->to_data() ) ;
                 $$ = PEL_Lexical::create( lst ) ;
              }
  | liste_args ',' something
              {
                 PEL_List * lst = $1->to_list() ;
                 $3->change_owner(lst,$3->to_data() ) ;
                 lst->append( $3->to_data() ) ;
                 $$ = $1 ;
              }

unary_operator: '-' something %prec UMINUS{
                    if( $2->to_data()->is_raw_data() )
                    {
                       PEL_Data::Type dt = $2->to_data()->data_type() ;
                       
                       if( dt == PEL_Data::Double )
                       {
                          $$=PEL_Lexical::create(
                             PEL_Double::create( 0, - $2->to_data()->to_double() ) ) ;
                       }
                       else if( dt == PEL_Data::Int )
                       {
                          $$=PEL_Lexical::create(
                             PEL_Int::create( 0, - $2->to_data()->to_int() ) ) ;
                       }
                       else
                       {
                          PELerror( "Undefined unary minus operator" ) ;
                       }
                    }
                    else
                    {
                       CREATE_SINGLE_OP("unary_minus",$2,$$) ; 
                    } }
              | '!' something %prec UNOT { CREATE_SINGLE_OP("!",$2,$$) ; }

binary_operator: something '<' something { CREATE_BIN_OP("<",$1,$3,$$) ; }
               | something '>' something { CREATE_BIN_OP(">",$1,$3,$$) ; }
               | something '+' something { CREATE_BIN_OP("+",$1,$3,$$) ; }
               | something '-' something { CREATE_BIN_OP("-",$1,$3,$$) ; }
               | something '*' something { CREATE_BIN_OP("*",$1,$3,$$) ; }
               | something '=' something { CREATE_BIN_OP("=",$1,$3,$$) ; }
               | something '/' something { CREATE_BIN_OP("/",$1,$3,$$) ; }
               | something PEL__NEQ something { CREATE_BIN_OP("!=",$1,$3,$$) ; }
               | something PEL__CONCAT something { CREATE_BIN_OP("<<",$1,$3,$$) ; }
               | something PEL__OR  something { CREATE_BIN_OP("||",$1,$3,$$) ; }
               | something PEL__AND something { CREATE_BIN_OP("&&",$1,$3,$$) ; }
               | something PEL__LE something { CREATE_BIN_OP("<=",$1,$3,$$) ; }
               | something PEL__GE something { CREATE_BIN_OP(">=",$1,$3,$$) ; }

module_def : module
         |  error

module : module_deb free_list module_fin
         {
            PEL_ASSERT( $3==0 || YY_top_module!=$3->to_module() ) ;
            $$ = $3 ;
         }

if_module : { $$ = PEL_Lexical::create( PEL_Bool::create( 0, true ) ) ; }
          |  PEL__IF '(' something ')'
         {
            if( $3->to_data()->data_type() != PEL_Data::Bool )
               PELerror( "Conditional module if must refer to boolean expression" ) ;
            bool cond = false ;
            if( YY_top_module!=dummy_module )
            {
               if( ! $3->to_data()->value_can_be_evaluated(
                                             YY_top_module->context() ) )
               {
                  std::string msg = "Conditional module if can not be evaluated\n" ;
                  msg += "Undefined variable(s):\n" ;
                  stringVector const& undef =
                     $3->to_data()->undefined_variables(
                                                  YY_top_module->context() ) ;
                  for( size_t i=0 ; i<undef.size() ; ++i )
                  {
                     msg += "   - \""+undef(i)+"\"\n" ;
                  }
                  PELerror( msg.c_str() ) ;
               }
               cond = $3->to_data()->to_bool( YY_top_module->context() ) ;
            }
            $$ = PEL_Lexical::create( PEL_Bool::create( 0, cond ) ) ;
         }

module_deb : if_module PEL__MODULE PEL__IDENTIF
       {
          modules_LILO->append( YY_top_module ) ;
          if( $1->to_data()->to_bool() && YY_top_module!=dummy_module )
          {
             PEL_Module* old = YY_top_module ;
             std::string a_name = $3->to_data()->to_string() ;
             std::string path = "/" + a_name ;
             if( old->has_module( path ) )
             {
                YY_top_module = old->module( path ) ;
             }
             else
             {
                YY_top_module = PEL_Module::create( old, a_name ) ;
                if( old->has_entry( a_name ) )
                {
                   std::string mess = "Can't create entry and module with the same name "+a_name+" in module " + old->name() ;
                   PELerror( mess.c_str() ) ;
                }
                old->add_module( YY_top_module ) ;
             }
          }
          else
          {
             YY_top_module = dummy_module ;
          }
       }


module_fin : PEL__END PEL__MODULE  PEL__IDENTIF
       {
          if( YY_top_module!=dummy_module )
          {
             if( $3->to_data()->to_string()!=YY_top_module->name() )
             {
                std::string mess = "Module " + YY_top_module->name() +
                   " can't be closed with " + $3->to_data()->to_string() ;
                PELerror( mess.c_str() ) ;
             }
             $$=PEL_Lexical::create( YY_top_module ) ;
          }
          else
          {
             $$=0 ;
          }
          size_t last = modules_LILO->count()-1 ;
          YY_top_module = static_cast<PEL_Module*>( modules_LILO->at( last ) ) ;
          modules_LILO->remove_at( last ) ;
          if( modules_LILO->has( dummy_module ) )
          {
             YY_top_module = dummy_module ;
          }
       }  
    
   /*
    * Vector specifications
    */
vector : '<' simple_item_list '>'
         {
            PEL_List * lst = $2->to_list() ;
            PEL_Data* res=0 ;
            if( lst->count() > 0  )
            {
               PEL_Data const* item =
                  dynamic_cast<PEL_Data const*>(lst->at( 0 )) ;
               if( item!=0  )
               {
                  switch( item->data_type() )
                  {
                     case PEL_Data::Double : 
                        res = PEL_DoubleVector::create( 0, lst ) ;
                        break ;
                     case PEL_Data::Int : 
                        res = PEL_IntVector::create( 0, lst ) ;
                        break ;
                     case PEL_Data::Bool : 
                        res = PEL_BoolVector::create( 0, lst ) ;
                        break ;
                     case PEL_Data::String : 
                        res = PEL_StringVector::create( 0, lst ) ;
                        break ;
                     default :
                        break ;
                  }
               }
            }
            if( res==0 )
            {
		PELerror( "invalid list of values enclosed in < .. >" ) ;
            }
            $$ = PEL_Lexical::create( res ) ;
         }

simple_item_list : { $$ = PEL_Lexical::create( PEL_List::create( 0 ) ) ; }
	| simple_item_list item_vector
          {
             $1->to_list()->append( $2->to_data() ) ;
             $$=$1;
          }

array : '[' simple_vector_list ']'
{
            PEL_List * lst = $2->to_list() ;
            PEL_Data* res=0 ;
            if( lst->count() > 0  )
            {
               PEL_Data const* item =
                  dynamic_cast<PEL_Data const*>(lst->at( 0 )) ;
               if( item!=0  )
               {
                  switch( item->data_type() )
                  {
                     case PEL_Data::DoubleVector : 
                        res = PEL_DoubleArray2D::create( 0, lst ) ;
                        break ;
                     case PEL_Data::IntVector : 
                        res = PEL_IntArray2D::create( 0, lst ) ;
                        break ;
                     case PEL_Data::BoolVector : 
                        res = PEL_BoolArray2D::create( 0, lst ) ;
                        break ;
                     case PEL_Data::StringVector : 
                        res = PEL_StringArray2D::create( 0, lst ) ;
                        break ;
                     default :
                        break ;
                  }
               }
            }
            if( res==0 )
            {
		PELerror( "invalid list of vectors enclosed in [ .. ]" ) ;
            }
            $$ = PEL_Lexical::create( res ) ;
}

simple_vector_list : vector
          { $$ = PEL_Lexical::create( PEL_List::create( 0 ) ) ;
            $$->to_list()->append($1->to_data() ) ;
          }
	| simple_vector_list ',' vector
          {
             $1->to_list()->append( $3->to_data() ) ;
             $$=$1;
          }

item_vector: unary_operator | SimpleType
/*
 * Scalar types
 */

SimpleType: PEL__REAL
          | PEL__INTEGER
          | PEL__STRING
          | PEL__TRUE      { $$ = PEL_Lexical::create( PEL_Bool::create( 0, true ) ) ; }
	  | PEL__FALSE     { $$ = PEL_Lexical::create( PEL_Bool::create( 0, false ) ) ; }      

%%
// Read a file        
// read recursivly data file throught yyparse function
//----------------------------------------------------------------------
bool PEL_readFile( PEL_Module * top,
                   istream* input_stream,
                   std::string const& name,
                   bool debug )
//----------------------------------------------------------------------
{
   PEL_LABEL( "Gram.y::PEL_readFile" ) ;
   PEL_ASSERT( EQUIVALENT( input_stream!=0, name.empty() ) ) ;
   PEL_ASSERT( top!=0 ) ;
   PEL_ASSERT( IMPLIES( input_stream!=0, input_stream->good() ) ) ;

   if( parsing ) // Direct recursive called is forbidden
   {
      PELerror( "Try to read a new file, while a data file is already being parsed." ) ;
   }
   PEL_flex_debug = ( debug ? 1 : 0 ) ;
   CTX_INI = top->context() ;
   bool result = readFile( top, input_stream, name, true ) ;
   CTX_INI = 0 ;

   return( result ) ;
}

//----------------------------------------------------------------------
bool readFile( PEL_Module * top,
               istream* input_stream,
               std::string const& name,
               bool main_read )
//----------------------------------------------------------------------
{
   PEL_LABEL( "Gram.y::readFile" ) ;
   PEL_ASSERT( EQUIVALENT( input_stream!=0, name.empty() ) ) ;
   PEL_ASSERT( top!=0 ) ;
   PEL_ASSERT( IMPLIES( input_stream!=0, input_stream->good() ) ) ;
   
   YY_top_module = top ;
   if( modules_LILO==0 ) modules_LILO = PEL_ListIdentity::create( 0 ) ;
   if( dummy_module==0 ) dummy_module = PEL_Module::create( 0, "dummy_module" ) ;
   if( main_read )
   {
      parsing = true ;
      relativePath = "." ;
   }
   
   main_stack.push( main_read ) ;
   path_stack.push( relativePath ) ;
   nb_line_stack.push( PEL__NbLines ) ;
   name_stack.push( currentFile ) ;
   PEL__NbLines = 1 ;
   
   istream* newFile = 0 ;
   std::ifstream* createdStream = 0 ;
   if( input_stream!=0 )
   {
      newFile = input_stream ;
      currentFile = "stream reading" ;
   }
   else if( name=="-stdin" )
   {
      newFile = &std::cin ;
      currentFile = "standard input stream" ;
   }
   else
   {
      char separator = PEL_System::path_name_separator() ;
      
      bool absolute_dir = name.find_first_of( separator )==0 ||
         // Special case for windows-like path name
         ( name.length() > 2 && name[1]==':' && name[2]=='\\' ) ;
      size_t id = name.find_last_of( separator ) ;
      if( absolute_dir )
      {
         currentFile = name ;
         relativePath = name.substr( 0, id ) ;
      }
      else
      {
         currentFile = relativePath + separator + name ;
         if( id < name.length() )
         {
            relativePath = relativePath + separator + name.substr( 0, id ) ;
         }
      }
      newFile = createdStream = new std::ifstream( currentFile.c_str() ) ;
      if( newFile->fail() )
         PEL_Error::object()->raise_plain( "Unable to open file "+currentFile ) ;
      
   }
   
   file_stack.push(newFile) ;
   switch_to_buffer( newFile ) ;
   
   if( main_read )
   {
      PELparse() ;
      parsing = false ;
   }

   if( createdStream!=0 && main_read ) // Other streams are destroyed in endFile
   {
      createdStream->close() ;
      delete createdStream ; createdStream=0 ;
   }
   
   return true ;
   
}

//----------------------------------------------------------------------
bool endFile( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "Gram.y::endFile" ) ;
   PEL_ASSERT( !main_stack.empty() ) ;
   PEL_ASSERT( !path_stack.empty() ) ;
   PEL_ASSERT( !nb_line_stack.empty() ) ;
   PEL_ASSERT( !name_stack.empty() ) ;
   PEL_ASSERT( !file_stack.empty() ) ;
   
   un_switch_to_buffer() ;
   bool main_read = main_stack.top() ; main_stack.pop() ;
   relativePath = path_stack.top() ; path_stack.pop() ;
   PEL__NbLines = nb_line_stack.top() ; nb_line_stack.pop() ;
   currentFile = name_stack.top() ; name_stack.pop() ;
   istream* file = file_stack.top() ; file_stack.pop() ;
   
   bool res = false ;
   
   if( !main_read )
   { 
      delete file ; file = 0 ;
      res = false ;
   }
   else
   {
      PEL_ASSERT( main_stack.empty() ) ;
      PEL_ASSERT( path_stack.empty() ) ;
      PEL_ASSERT( nb_line_stack.empty() ) ;
      PEL_ASSERT( name_stack.empty() ) ;
      PEL_ASSERT( file_stack.empty() ) ;
      if( modules_LILO->count()!=0 )
      {
         std::string mess = "When reading PELICANS data structure, following modules are not correclty closed\n" ;
         for( size_t i=0 ; i<modules_LILO->count() ; i++ )
            mess += "\"" + static_cast<PEL_Module*>( modules_LILO->at( i ) )->name() + "\"\n" ;
         PELerror( mess.c_str() ) ;
      }
      
      PEL_Lexical::remove_all_lexical() ;
      modules_LILO->destroy() ; modules_LILO=0 ;
      dummy_module->destroy() ; dummy_module = 0 ;
      res = true ;
   }
   
   return res ;
   
}

//----------------------------------------------------------------------
void PEL_re_init_parser( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_re_init_parser" ) ;
   
   while(!main_stack.empty()) main_stack.pop() ;
   while(!path_stack.empty()) path_stack.pop() ;
   while(!nb_line_stack.empty()) nb_line_stack.pop() ;
   while(!name_stack.empty()) name_stack.pop() ;
   while(!file_stack.empty()) {
      //if( file_stack.top()!=0 ) delete file_stack.top() ;
      file_stack.pop() ; 
   }
   
   PEL_re_init_lexer() ;
   PEL_Lexical::remove_all_lexical() ;
   if(modules_LILO!=0 ) {
      modules_LILO->destroy() ;
      modules_LILO=0 ;
   }
   if( dummy_module!=0 ) 
   {      
      dummy_module->destroy() ;
      dummy_module = 0 ;
   }
   parsing = false ;
   
}

//
// Buffering used to help in error case 
//
//----------------------------------------------------------------------
void PEL__Buffer( std::string const& chain ) 
//----------------------------------------------------------------------
{
   size_t cr=chain.find('\n') ;
   size_t cr_last = 0 ;
   for( ;
        cr<chain.length() ;
        cr=chain.find('\n', cr+1 ) )
   {
      PEL__NbLines++ ;
      cr_last = cr+1 ;
      buff="" ;
   }
   buff += chain.substr( cr_last, chain.length()-cr_last ) ;
}

//----------------------------------------------------------------------
void PELerror( const char * s)
//----------------------------------------------------------------------
{
   PEL_re_init_parser() ;
   
   PEL_Error::object()->raise_read_syntax_error(
      currentFile,
      PEL__NbLines,
      buff,
      s ) ;
}

//----------------------------------------------------------------------
std::string comment( void )
//----------------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl
       << "Last line number " << PEL__NbLines
       << " read in file " << currentFile << " was:" << std::endl
       << ">> " ;
   if( buff.length()<40 )
   {
      msg << buff << " <<" << std::endl ;
   }
   else
   {
      msg << buff.substr( 0, 40 ) << " ... (TO BE CONTINUED) <<" << std::endl ;
   }
   
   return msg.str() ;
}

//----------------------------------------------------------------------
std::string PEL_current_parsed_module_path_name( void )
//----------------------------------------------------------------------
{
   std::string result = "" ;
   if( YY_top_module!=dummy_module && YY_top_module!=0 )
   {      
      result = YY_top_module->absolute_path_name() ;
   }
   return result ;
}

//----------------------------------------------------------------------
int PEL_current_parsed_line( void )
//----------------------------------------------------------------------
{
   return PEL__NbLines ;
}

//----------------------------------------------------------------------
void substitute_for_assignment( std::string const& key,
                                PEL_Data* data )
//----------------------------------------------------------------------
{
   PEL_CHECK( YY_top_module!=dummy_module ) ;
   
   if( modules_LILO->index_limit()<1 )
   {
      std::string mess = "No path to apply substitution of " ;
      mess += key ;
      PELerror( mess.c_str() ) ;
   }
   
   PEL_Module * mod = static_cast<PEL_Module *>( modules_LILO->at( 0 ) ) ;
   PEL_Iterator * it = modules_LILO->create_iterator( 0 ) ;
   it->go_next() ;
   
   for( ; it->is_valid() ; it->go_next() )
   {
      PEL_Module * child = static_cast<PEL_Module *>( it->item() ) ;
      if( !mod->has_module( child->name() ) )
      {
         std::string mess = "No path to apply substitution of " ;
         mess += key + " for module " + child->name() ;
         mod->print( PEL::out(), 0 ) ;
         PELerror( mess.c_str() ) ;
      }
      mod = mod->module( child->name() ) ;
   }
   if( !mod->has_module( YY_top_module->name() ) )
   {
      std::string mess = "No path to apply substitution of " ;
      mess += key + " for module " + YY_top_module->name() ;
      mod->print( PEL::out(), 0 ) ;
      
      PELerror( mess.c_str() ) ;
   }
   mod = mod->module( YY_top_module->name() ) ;
   std::string root_key = "/" +key ;
   if( !mod->has_entry( root_key ) )
   {
      std::string mess = key + " doesn't already exist " ;
      PELerror( mess.c_str() ) ;
   }
   if( data->owner()!=0 ) data = data->create_clone(0) ;
   
   data->set_owner(mod) ;
   mod->replace_data_of_entry( root_key, data ) ;
   
   it->destroy() ;
      
}
