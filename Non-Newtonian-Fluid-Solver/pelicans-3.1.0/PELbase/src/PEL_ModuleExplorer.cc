/*
 *  Copyright 1995-2010 by IRSN
 *
 *  This software is an application framework, with a set of integrated  
 *  reusable components, whose purpose is to simplify the task of developing 
 *  softwares of numerical mathematics and scientific computing.
 * 
 *  This software is governed by the CeCILL-C license under French law and 
 *  abiding by the rules of distribution of free software. You can use, modify 
 *  and/or redistribute the software under the terms of the CeCILL-C license  
 *  as circulated by CEA, CNRS and INRIA at the following URL 
 *  "http://www.cecill.info". 
 *
 *  As a counterpart to the access to the source code and rights to copy,  
 *  modify and redistribute granted by the license, users are provided only 
 *  with a limited warranty and the software's author, the holder of the  
 *  economic rights, and the successive licensors have only limited liability. 
 *
 *  In this respect, the user's attention is drawn to the risks associated  
 *  with loading, using, modifying and/or developing or reproducing the  
 *  software by the user in light of its specific status of free software,
 *  that may mean that it is complicated to manipulate, and that also  
 *  therefore means that it is reserved for developers and experienced 
 *  professionals having in-depth computer knowledge. Users are therefore 
 *  encouraged to load and test the software's suitability as regards their 
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and, more generally, to use and operate it in the same 
 *  conditions as regards security. 
 *
 *  The fact that you are presently reading this means that you have had 
 *  knowledge of the CeCILL-C license and that you accept its terms.
 */

#include <PEL_ModuleExplorer.hh>

#include <PEL_assertions.hh>
#include <PEL_Context.hh>
#include <PEL_ContextPair.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_List.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_ModulePattern.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>
#include <boolVector.hh>
#include <boolArray2D.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>
#include <stringArray2D.hh>

#include <fstream>
#include <sstream>
using std::string ;

//--------------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_ModuleExplorer:: create( PEL_Object* a_owner,
                             PEL_Module const* mm,
                             PatternStatus status )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: create" ) ;
   PEL_CHECK_PRE( mm != 0 ) ;
   PEL_CHECK_PRE( IMPLIES( status==verify || status==build,
                           PEL_ModulePattern::has_pattern() ) ) ;
   PEL_CHECK_PRE( IMPLIES( status==build,
                           PEL_ModulePattern::build_pattern() ) ) ;

   PEL_ModulePattern* pat = 0 ;
   if( status!=ignore )
   {
      pat = PEL_ModulePattern::create( 0, mm, mm->name() ) ;
   }
   
   PEL_ModuleExplorer* result =
                       new PEL_ModuleExplorer( a_owner, mm, status, pat ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->pattern_status() == status ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_ModuleExplorer:: PEL_ModuleExplorer( PEL_Object* a_owner,
                                         PEL_Module const* mm,
                                         PatternStatus status,
                                         PEL_ModulePattern* pat )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner ),
     mod( mm ),
     submodule_iterator( mm->create_module_iterator( this) ),
     keyword_iterator( mm->create_entry_iterator( this) ),
     tmp_context( PEL_ContextPair::create( this,0,0) ),
     MP( pat ),
     keyword_iterator_started(false),
     my_status( status )
{
   if( pat!=0 ) pat->set_owner( this ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//--------------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_ModuleExplorer:: create_subexplorer( PEL_Object* a_owner,
                                         std::string const& path_and_name ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: create_subexplorer( PEL_Object*, std::string const& )" ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   if( !mod->has_module( path_and_name ) ) 
   {
      PEL_Error::object()->raise_missing_module( this, path_and_name ) ;
   }
   PEL_Module const* mm = mod->module( path_and_name ) ;
   PEL_ModulePattern* pat = 0 ;
   if( pattern()!=0 )
   {
      pat = pattern()->create_subpattern( 0, mod, path_and_name ) ;
   }
   
   PEL_ModuleExplorer* result =
                      new PEL_ModuleExplorer( a_owner, mm, my_status, pat ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->pattern_status() == pattern_status() ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK( result->mod == mm ) ;
   PEL_CHECK_POST( result->name() == PEL_Module::basename( path_and_name ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_ModuleExplorer:: create_clone( PEL_Object* a_owner ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: create_clone" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_ModulePattern* pat = 0 ;
   if( pattern()!=0 )
   {
      pat = pattern()->create_clone( 0 ) ;
   }
   PEL_ModuleExplorer* result =
                     new PEL_ModuleExplorer( a_owner, mod, my_status, pat ) ;
   
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   PEL_CHECK_POST( result->pattern_status() == pattern_status() ) ;
   PEL_CHECK_INV( invariant() ) ;
   return result ;
}

//---------------------------------------------------------------------------
PEL_ModuleExplorer:: ~PEL_ModuleExplorer( void )
//---------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   mod = 0 ;
   MP = 0 ;
}

//--------------------------------------------------------------------------
PEL_ModuleExplorer::PatternStatus
PEL_ModuleExplorer:: pattern_status( void ) const
//--------------------------------------------------------------------------
{
   if( my_status != ignore && !PEL_ModulePattern::has_pattern() )
   {
      PEL_Error::object()->raise_internal(
         "*** PEL_ModuleExplorer error:\n"
         "    pattern file has already been closed" ) ;
   }
   return( my_status ) ;
}

//---------------------------------------------------------------------------
PEL_ModulePattern*
PEL_ModuleExplorer:: pattern( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: pattern" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_ModulePattern* result = MP ;
   
   PEL_CHECK_POST( IMPLIES( build_pattern(), result!=0 ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
bool
PEL_ModuleExplorer:: build_pattern( void ) const
//--------------------------------------------------------------------------
{
   bool result = my_status == build &&
                 !PEL_Assertion::is_checking() ;
   if( my_status == build && !PEL_ModulePattern::build_pattern() )
   {
      PEL_Error::object()->raise_internal(
         "*** PEL_ModuleExplorer error:\n"
         "    pattern file has already been closed" ) ;
   }
   
   return( result ) ;
}

//--------------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_ModuleExplorer:: validity( PEL_Object* a_owner,
                               bool recurse ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: validity" ) ;
   PEL_CHECK_PRE( pattern_status()==verify ) ;
   PEL_CHECK_PRE( PEL_ModulePattern:: has_pattern() ) ;

   PEL_ModuleExplorer* result  = 0 ;
   PEL_ModulePattern const* pat =
      PEL_ModulePattern::create_pattern_description( 0, mod->name() ) ;
   if( pat == 0 )
   {
      PEL_Error::object()->raise_plain(
         "*** pattern verification failed\n"
         "    unable to find description of:\n"
         "       module: "+mod->name()+"\n"
         "    the control file may not exist or may be corrupted" ) ;
   }
   else
   {
      PEL_Module* m = pat->validity( 0, mod, recurse ) ;
      result  = PEL_ModuleExplorer::create( a_owner, m ) ;
      m->set_owner( result ) ;
      pat->destroy() ; pat = 0 ;
   }

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;   
   return( result ) ;
}

//---------------------------------------------------------------------------
bool 
PEL_ModuleExplorer:: has_module( std::string const& a_path_and_name ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: has_module" ) ;
   PEL_CHECK_PRE( !a_path_and_name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = mod->has_module( a_path_and_name ) ;

   if( build_pattern() && result )
   {
      std::string dir_name = PEL_Module::dirname( a_path_and_name ) ;
      std::string base_name = PEL_Module::basename( a_path_and_name ) ;
      PEL_ASSERT( !base_name.empty() ) ;
      if( dir_name.empty() )
      {
         pattern()->add_pattern( a_path_and_name,
                                 mod->module( a_path_and_name ),
                                 PEL_ModulePattern::optional ) ;
      }
      else
      {
         PEL_ModuleExplorer const* exp = create_subexplorer( 0, dir_name ) ;
	 PEL_ASSERT( exp->has_module( base_name ) ) ;
	 exp->destroy() ; exp = 0 ;
      }
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ModuleExplorer:: is_empty( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: is_empty" ) ;
   return( mod->is_empty() ) ;
}

//---------------------------------------------------------------------------
PEL_Context const*
PEL_ModuleExplorer:: context( std::string const& a_path_and_name,
                              PEL_Context const* ct ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: context" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Module const* m = mod ;
   if( !a_path_and_name.empty() )
   {
      m = mod->module( a_path_and_name ) ;
   }
   PEL_Context const* result = m->context() ;
   if( ct!=0 )
   {
      if( build_pattern() )
      {
         stringVector vars( ct->nb_variables() ) ;
         for( size_t i=0 ; i<ct->nb_variables() ; ++i )
         {
            PEL_Variable const* var = ct->variable(i) ;
            vars(i) = var->name() ;
         }
         vars.sort() ;
         for( size_t i=0 ; i<vars.size() ; ++i )
         {
            std::string const& str = vars(i) ;
            PEL_Variable const* var = PEL_Variable::object( str ) ;
            pattern()->add_provided_variable( str, m, var->data_type() ) ;
         }
      }
      tmp_context->re_initialize( m->context(), ct ) ;
      result = tmp_context ;
   }
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
bool 
PEL_ModuleExplorer:: has_entry( std::string const& a_path_and_keyword ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: has_entry" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = mod->has_entry( a_path_and_keyword ) ;
   if( result && build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::Undefined ) ;      
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: start_module_iterator( void ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: start_module_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   submodule_iterator->start() ;
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ModuleExplorer:: is_valid_module( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: is_valid_module" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return submodule_iterator->is_valid() ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: go_next_module( void ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: go_next_module" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( is_valid_module() ) ;
   submodule_iterator->go_next() ;
}

//---------------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_ModuleExplorer:: create_subexplorer( PEL_Object* a_owner ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: create_subexplorer( PEL_Object* )" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( is_valid_module() ) ;

   PEL_ModulePattern* pat = 0 ;
   if( pattern()!=0 )
   {
      pat = pattern()->generic_pattern( 0, submodule_iterator->item() ) ;
   }
   
   PEL_ModuleExplorer* result =
      new PEL_ModuleExplorer( a_owner, submodule_iterator->item(), my_status, pat ) ;
   
   PEL_CHECK_POST( has_module( result->name() ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   return result ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: start_entry_iterator( void ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: start_entry_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   keyword_iterator->start() ;
   keyword_iterator_started = true ;
}

//---------------------------------------------------------------------------
bool
PEL_ModuleExplorer:: is_valid_entry( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: is_valid_entry" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return keyword_iterator->is_valid() ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: go_next_entry( void ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: go_next_entry" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( is_valid_entry() ) ;
   keyword_iterator->go_next() ;
}

//---------------------------------------------------------------------------
std::string const&
PEL_ModuleExplorer:: keyword( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: keyword" ) ;
   PEL_CHECK_PRE( is_valid_entry() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   string const& result = keyword_iterator->item()->keyword() ;
   
   PEL_CHECK_POST( has_entry( result ) ) ;
   return result ;
}

//---------------------------------------------------------------------------
PEL_DataWithContext*
PEL_ModuleExplorer:: data( PEL_Object* a_owner,
                           PEL_Context const* ct ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: data" ) ;
   PEL_CHECK_PRE( is_valid_entry() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Data const* a_data = keyword_iterator->item()->data() ;
   PEL_DataWithContext* result =
      PEL_DataWithContext::create( a_owner,
                                   a_data,
                                   context( "", ct ) ) ;
   if( build_pattern() )
   {
      pattern()->add_generic_keyword( keyword_iterator->item()->keyword(),  mod, a_data->data_type() ) ;
   }
   
   PEL_CHECK_POST( result !=0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_DataWithContext*
PEL_ModuleExplorer:: abstract_data( PEL_Object* a_owner,
                                    std::string const& a_path_and_keyword,
                                    PEL_Context const* ct ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: abstract_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   PEL_Data const* a_data = mod->data_of_entry( a_path_and_keyword ) ;
   
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, a_data->data_type() ) ;
   }
   PEL_DataWithContext* result =
      PEL_DataWithContext::create(
         a_owner,
         a_data,
         context( PEL_Module::dirname( a_path_and_keyword ), ct ) ) ;
   
   PEL_CHECK_POST( result !=0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ModuleExplorer:: bool_data( std::string const& a_path_and_keyword,
                                PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: bool_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   PEL_Data const* val = mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::Bool ) ;
   }
   if( val->data_type()!=PEL_Data::Bool )
   {
      PEL_Error::object()->raise_bad_data_type(
                              this, a_path_and_keyword, PEL_Data::Bool  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   return( val->to_bool( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
int
PEL_ModuleExplorer:: int_data( std::string const& a_path_and_keyword,
                               PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: int_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   PEL_Data const* val = mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::Int ) ;
   }
   if( val->data_type()!=PEL_Data::Int )
   {
      PEL_Error::object()->raise_bad_data_type(
                                this, a_path_and_keyword, PEL_Data::Int ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_int( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
double
PEL_ModuleExplorer:: double_data( std::string const& a_path_and_keyword,
                                  PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: double_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   PEL_Data const* val = mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::Double ) ;
   }
   if( val->data_type()!=PEL_Data::Double )
   {
      PEL_Error::object()->raise_bad_data_type(
                               this, a_path_and_keyword, PEL_Data::Double ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
      
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_double( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
std::string const&
PEL_ModuleExplorer:: string_data( std::string const& a_path_and_keyword,
                                  PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: string_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::String ) ;
   }
   if( val->data_type()!=PEL_Data::String )
   {
      PEL_Error::object()->raise_bad_data_type(
                                this, a_path_and_keyword, PEL_Data::String  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_string( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
intVector const&
PEL_ModuleExplorer:: intVector_data( std::string const& a_path_and_keyword,
                                     PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: intVector_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::IntVector ) ;
   }
   if( val->data_type()!=PEL_Data::IntVector )
   {
      PEL_Error::object()->raise_bad_data_type(
                            this, a_path_and_keyword, PEL_Data::IntVector  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
              this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_int_vector( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
doubleVector const&
PEL_ModuleExplorer:: doubleVector_data( std::string const& a_path_and_keyword,
                                        PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: doubleVector_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::DoubleVector ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=PEL_Data::DoubleVector )
   {
      PEL_Error::object()->raise_bad_data_type(
                         this, a_path_and_keyword, PEL_Data::DoubleVector  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
              this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_double_vector( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
doubleArray2D const&
PEL_ModuleExplorer:: doubleArray2D_data( std::string const& a_path_and_keyword,
                                         PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: doubleArray2D_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::DoubleArray2D ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=PEL_Data::DoubleArray2D )
   {
      PEL_Error::object()->raise_bad_data_type(
                        this, a_path_and_keyword, PEL_Data::DoubleArray2D  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_double_array2D( a_ctx ) ) ;
}


//---------------------------------------------------------------------------
boolArray2D const&
PEL_ModuleExplorer:: boolArray2D_data( std::string const& a_path_and_keyword,
                                         PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: boolArray2D_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::BoolArray2D ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=PEL_Data::BoolArray2D )
   {
      PEL_Error::object()->raise_bad_data_type(
                        this, a_path_and_keyword, PEL_Data::BoolArray2D  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_bool_array2D( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
stringArray2D const&
PEL_ModuleExplorer:: stringArray2D_data( std::string const& a_path_and_keyword,
                                         PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: stringArray2D_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::StringArray2D ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=PEL_Data::StringArray2D )
   {
      PEL_Error::object()->raise_bad_data_type(
                        this, a_path_and_keyword, PEL_Data::StringArray2D  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_string_array2D( a_ctx ) ) ;
}


//---------------------------------------------------------------------------
intArray2D const&
PEL_ModuleExplorer:: intArray2D_data( std::string const& a_path_and_keyword,
                                      PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: intArray2D_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::IntArray2D ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=PEL_Data::IntArray2D )
   {
      PEL_Error::object()->raise_bad_data_type(
                           this, a_path_and_keyword, PEL_Data::IntArray2D  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_int_array2D( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
doubleArray3D const&
PEL_ModuleExplorer:: doubleArray3D_data( std::string const& a_path_and_keyword,
                                         PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: doubleArray3D_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::DoubleArray3D ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=PEL_Data::DoubleArray3D )
   {
      PEL_Error::object()->raise_bad_data_type(
                      this, a_path_and_keyword, PEL_Data::DoubleArray3D  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_double_array3D( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
intArray3D const&
PEL_ModuleExplorer:: intArray3D_data( std::string const& a_path_and_keyword,
                                      PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: intArray3D_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::IntArray3D ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=PEL_Data::IntArray3D )
   {
      PEL_Error::object()->raise_bad_data_type(
                          this, a_path_and_keyword, PEL_Data::IntArray3D  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_int_array3D( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
boolVector const&
PEL_ModuleExplorer:: boolVector_data( std::string const& a_path_and_keyword,
                                      PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: boolVector_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::BoolVector ) ;
   }
   if( val->data_type()!=PEL_Data::BoolVector )
   {
      PEL_Error::object()->raise_bad_data_type(
                           this, a_path_and_keyword, PEL_Data::BoolVector  ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
              this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_bool_vector( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
stringVector const&
PEL_ModuleExplorer:: stringVector_data( std::string const& a_path_and_keyword,
                                        PEL_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: stringVector_data" ) ;
   PEL_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   PEL_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, PEL_Data::StringVector ) ;
   }
   if( val->data_type()!=PEL_Data::StringVector )
   {
      PEL_Error::object()->raise_bad_data_type(
                          this, a_path_and_keyword, PEL_Data::StringVector ) ;
   }
   PEL_Context const* a_ctx =
                  context( PEL_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      PEL_Error::object()->raise_not_evaluable(
               this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   return( val->to_string_vector( a_ctx ) ) ;
}

//---------------------------------------------------------------------------
PEL_Module*
PEL_ModuleExplorer:: create_clone_of_attached_module(
   PEL_Object* a_owner ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: create_clone_of_attached_module" ) ;

   PEL_Module* result = mod->create_clone( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   
   return result ;
}

//---------------------------------------------------------------------------
std::string const&
PEL_ModuleExplorer:: name( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::string const& result = mod->name() ;
   PEL_CHECK_POST( !result.empty() ) ;
   PEL_CHECK_POST( PEL_Module::dirname( result ).empty() ) ;
   PEL_CHECK_POST( result==PEL_Module::basename( absolute_path_name() ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
std::string const&
PEL_ModuleExplorer:: absolute_path_name( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: absolute_path_name" ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::string const& result = mod->absolute_path_name() ;
   PEL_CHECK_POST( !result.empty() ) ;
   PEL_CHECK_POST( PEL_Module::basename( result )==name() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_Object const*
PEL_ModuleExplorer:: owner_of_attached_module( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: owner_of_attached_module" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( mod->owner() ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   mod->print( os, indent_width ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: write( std::string const& file,
                            std::string const& format ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: write" ) ;
   PEL_CHECK_INV( invariant() ) ;
   mod->write( file, format ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: test_data( std::string const& a_keyword,
                                std::string const& expression,
                                PEL_Context const* ct ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: test_data" ) ;
   PEL_CHECK_PRE( !a_keyword.empty() ) ;
   PEL_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   PEL_CHECK_PRE( !expression.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   if( !mod->has_entry( a_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_keyword ) ;
   }

   PEL_Data const* eval = mod->create_evaluation( 0, expression, ct ) ;
   if( eval==0 || eval->data_type()!=PEL_Data::Bool )
   {
      PEL_Error::object()->raise_data_error(
         this, a_keyword,
         "   condition ("+expression+")\n"
         "   can't be evaluated or is not a boolean expression" ) ;
   }

   if( !eval->to_bool() )
   {
      std::ostringstream msg ;
      msg << "   condition ("+expression+")\n" ;
      msg << "   is not verified:\n" ;
      msg << "      (" ;
      eval->print( msg, 0 ) ;
      msg << ")" ;
      PEL_Error::object()->raise_data_error(
                                     this , a_keyword, msg.str() ) ;
   }
   eval->destroy() ;

   if( build_pattern() )
   {
      pattern()->attach_verify_data( a_keyword, expression ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: test_data_in( std::string const& a_keyword,
                                   stringVector const& choices ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: test_data_in" ) ;
   PEL_CHECK_PRE( !a_keyword.empty() ) ;
   PEL_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   PEL_CHECK_PRE( choices.size() > 0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_keyword ) ;
   }

   PEL_Data const* eval = mod->data_of_entry( a_keyword ) ;

   if( eval->data_type() != PEL_Data::String &&
       eval->data_type() != PEL_Data::Bool &&
       eval->data_type() != PEL_Data::Int  &&
       eval->data_type() != PEL_Data::StringVector &&
       eval->data_type() != PEL_Data::BoolVector &&
       eval->data_type() != PEL_Data::IntVector )
   {
      // Double are forbidden because the equality between two doubles
      // is not clear...
      PEL_Error::object()->raise_data_error(
         this, a_keyword,
         "*** PEL_ModuleExplorer:: test_data_in error\n"
         "    test only available for:\n"
         "        - string values\n"
         "        - boolean values\n"
         "        - integer values\n"
         "        - string vector values\n"
         "        - boolean vector values\n"
         "        - integer vector values" ) ;
   }

   std::string str = eval->value_as_string( mod->context() ) ;
   PEL::remove_enclosing_characters(str,'\"') ;
   
   if( !choices.has( str ) )
   {
      std::string mess ;
      for( size_t i=0 ; i<choices.size() ; i++ )
      {
         mess += "   - \""+choices(i)+"\"\n" ;
      }
      PEL_Error::object()->raise_bad_data_value(
         this, a_keyword, mess ) ;
   }

   if( build_pattern() && a_keyword != "type" && a_keyword != "concrete_name" )
   {
      pattern()->attach_list_of_valid_choices( a_keyword, choices ) ;
   }
   
   
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: test_data_as( std::string const& a_keyword,
                                   std::string const& regexp,
                                   std::string const& where ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: test_data_as" ) ;
   PEL_CHECK_PRE( !regexp.empty() ) ;
   PEL_CHECK_PRE( !a_keyword.empty() ) ;
   PEL_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_List* list = mod->create_data_selection( 0, regexp, 0, where ) ;
         
   stringVector choices(0) ;
   PEL_ListIterator * it = list->create_iterator( list ) ;
   for( it->start() ; it->is_valid() ; it->go_next() ) 
   {
      PEL_Data const* adata = static_cast<PEL_Data*>( it->item() ) ;
      if( adata->data_type()!=PEL_Data::String )
      {
         PEL_Error::object()->raise_bad_data_type(
            this, regexp, PEL_Data::String ) ;
      }
      choices.append( adata->to_string() ) ;
   }
   
   list->destroy() ;
   if( !mod->has_entry( a_keyword ) ) 
   {
      PEL_Error::object()->raise_missing_keyword( this, a_keyword ) ;
   }

   PEL_Data const* eval = mod->data_of_entry( a_keyword ) ;
   if( eval->data_type()!=PEL_Data::String )
   {
      PEL_Error::object()->raise_bad_data_type(
         this, a_keyword, PEL_Data::String ) ;
   }

   std::string const& str = eval->to_string( mod->context() ) ;
   
   if( !choices.has( str ) )
   {
      std::string mess ;
      for( size_t i=0 ; i<choices.size() ; i++ )
      {
         mess += "   - \""+choices(i)+"\"\n" ;
      }
      PEL_Error::object()->raise_bad_data_value(
         this, a_keyword, mess ) ;
   }

   if( build_pattern() )
   {
      pattern()->attach_list_of_dynamic_choices( a_keyword, regexp, where ) ;
   }
   
   
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: set_default( std::string const& a_keyword,
                                  std::string const& expression ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: set_default" ) ;
   PEL_CHECK_PRE( !a_keyword.empty() ) ;
   PEL_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   PEL_CHECK_PRE( !expression.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( build_pattern() )
   {
      if( !mod->has_entry( a_keyword ) ) 
      {
         PEL_Error::object()->raise_missing_keyword( this, a_keyword ) ;
      }
      pattern()->attach_default_data( a_keyword, expression ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: set_help( std::string const& a_keyword,
                               std::string const& expression ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: set_help" ) ;
   PEL_CHECK_PRE( !a_keyword.empty() ) ;
   PEL_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   PEL_CHECK_PRE( !expression.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( build_pattern() )
   {
      if( !mod->has_entry( a_keyword ) ) 
      {
         PEL_Error::object()->raise_missing_keyword( this, a_keyword ) ;
      }
      pattern()->attach_help_data( a_keyword, expression ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}


//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: test_file( std::string const& filename,
                                std::string const& mode ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: test_file" ) ;
   PEL_CHECK_PRE( !filename.empty() ) ;
   PEL_CHECK_PRE( filename.find("/") > filename.length() ) ;
   PEL_CHECK_PRE( mode=="read" || mode=="write" ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const& fn = string_data( filename ) ;
   
   if( ( mode=="read" && !( PEL_System::can_read( fn ) ) ) ||
       ( mode=="write" && !( PEL_System::can_write( fn ) ) ) )
   {
      PEL_Error::object()->raise_bad_file( this, filename, mode ) ;
   }  
   
   if( build_pattern() )
   {
      if( !mod->has_entry( filename ) ) 
      {
         PEL_Error::object()->raise_missing_keyword( this, filename ) ;
      }
      pattern()->attach_file_extension( filename, mode ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModuleExplorer:: declare_data( std::string const& a_path_and_keyword,
                                   PEL_Data::Type type ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModuleExplorer:: declare_data" ) ;
   PEL_CHECK( !a_path_and_keyword.empty() ) ;
   PEL_CHECK( my_status == build ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_ModulePattern::Access acc =  PEL_ModulePattern::mandatory ;
   if( type==PEL_Data::Undefined )
   {
      acc=PEL_ModulePattern::optional ;
   }
   
   if( keyword_iterator_started &&
       is_valid_entry() &&
       keyword()==PEL_Module::basename( a_path_and_keyword ) )
   {
      acc=PEL_ModulePattern::generic ;
   }
   
   pattern()->add_entry( a_path_and_keyword, mod, type, acc ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ModuleExplorer:: invariant( void ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( mod!=0 ) ;
   PEL_ASSERT( IMPLIES( my_status==build, MP!=0  ) ) ;
   PEL_ASSERT( IMPLIES( MP!=0, my_status!=ignore ) ) ;
   return( true ) ;
}
