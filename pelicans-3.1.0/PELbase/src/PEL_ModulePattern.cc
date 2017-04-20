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

#include <PEL_ModulePattern.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_BoolVector.hh>
#include <PEL_Communicator.hh>
#include <PEL_Context.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exceptions.hh>
#include <PEL_Exec.hh>
#include <PEL_Int.hh>
#include <PEL_IntVector.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_List.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_Root.hh>
#include <PEL_String.hh>
#include <PEL_StringVector.hh>
#include <PEL_System.hh>
#include <boolVector.hh>
#include <boolArray2D.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <intArray2D.hh>
#include <intVector.hh>
#include <stringArray2D.hh>
#include <stringVector.hh>
#include <fstream>
#include <sstream>

#define LL "**********************************************" << std::endl
using std::string ;

struct PEL_ModulePattern_ERROR
{
   static void n0( std::string const& pattern,
                   std::string const& error ) ;
   static void n1( std::string const& keyword,
                   std::string const& old,
                   std::string const& newv ) ;
} ;

//--------------------------------------------------------------------------
std::string PEL_ModulePattern:: pattern_filename = "" ;
bool PEL_ModulePattern:: HAS = false ;
bool PEL_ModulePattern:: BUILD = false ;
PEL_Module* PEL_ModulePattern:: MP_File = 0 ;
PEL_Module* PEL_ModulePattern:: MP_Base = 0 ;
PEL_Module* PEL_ModulePattern:: MP_Pattern = 0 ;
std::string const indirection_to_keyword = "_to" ;
std::string const type_keyword = "_type" ;
std::string const access_keyword = "_access" ;
std::string const name_keyword = "_name" ;
std::string const module_type_keyword = "S" ;
std::string const variable_type_keyword = "V" ;
std::string const conditional_type_keyword = "C" ;
std::string const condition_keyword = "_if" ;
std::string const in_keyword = "_in" ;
std::string const vector_in_keyword = "_vector_in" ;
std::string const select_in_keyword = "_select_in" ;
std::string const where_keyword = "_where" ;
std::string const test_keyword = "_test" ;
std::string const default_keyword = "_default" ;
std::string const help_keyword = "_help" ;
std::string const file_keyword = "_file" ;
std::string const unique_keyword = "_unique" ;
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
bool
PEL_ModulePattern:: has_pattern( void )
//--------------------------------------------------------------------------
{
   return( HAS ) ;
}

//--------------------------------------------------------------------------
bool
PEL_ModulePattern:: build_pattern( void )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: build_pattern" ) ;
   
   bool result = BUILD ;
   
   PEL_CHECK_POST( IMPLIES( result, has_pattern() ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
void
PEL_ModulePattern:: build_pattern_base( std::string const& filename )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: build_pattern_base" ) ;
   PEL_CHECK_PRE( !has_pattern() ) ;
   PEL_CHECK_PRE( !build_pattern() ) ;
   
   BUILD = true ;
   open_pattern_base( filename ) ;
   
   PEL_CHECK_POST( has_pattern() ) ;
   PEL_CHECK_POST( build_pattern() ) ;
}
//--------------------------------------------------------------------------
void
PEL_ModulePattern:: open_pattern_base( std::string const& filename )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: open_pattern_base" ) ;
   PEL_CHECK_PRE( !has_pattern() ) ;

   static std::string const revision = "2.0" ;
   
   HAS = true ;
   pattern_filename = filename ;
   std::ifstream in( filename.c_str() ) ;
   if( in )
   {
      in.close() ;
      MP_File = PEL_Module::create( PEL_Root::object(), "MAIN", filename ) ;
      if( ! MP_File->has_module( "Base" ) )
      {
         PEL_ModulePattern_ERROR:: n0(
            pattern_filename,
            "    In module: /MAIN\n"
            "    missing module: Base" ) ;
      }
      MP_Base = MP_File->module( "Base" ) ;
      check_evaluable( true, MP_Base, "revision", PEL_Data::String ) ;
      if( MP_Base->data_of_entry( "revision" )->to_string() != revision )
      {
         PEL_ModulePattern_ERROR:: n0(
            pattern_filename,
            "    bad revision encountered\n"
            "     file revision: "
                    +MP_Base->data_of_entry( "revision" )->to_string()+"\n"
            "     expected revision: "+revision ) ;
      }
      if( MP_Base->has_module( "Pattern" ) )
      {
         MP_Pattern = MP_Base->module( "Pattern" ) ;
      }
      else
      {
         MP_Pattern = PEL_Module::create( MP_Base, "Pattern" ) ;
         MP_Base->add_module( MP_Pattern ) ;
      }
   }
   else
   {
      MP_File = 0 ;
      MP_Base = PEL_Module::create( PEL_Root::object(),
                                    "Base" ) ;
      MP_Pattern = PEL_Module::create( MP_Base,
                                       "Pattern" ) ;
      MP_Base->add_module( MP_Pattern ) ;
      MP_Base->add_entry( "revision",
                          PEL_String::create( MP_Base, revision ) ) ;
   }
   
   PEL_CHECK_POST( has_pattern() ) ;
}
 
//--------------------------------------------------------------------------
void
PEL_ModulePattern:: close_pattern_base( void )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: close_pattern_base" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   
   HAS = false ;
   BUILD = false ;
   pattern_filename = "" ;
   if( MP_File != 0 )
   {
      PEL_Root::object()->destroy_possession( MP_File ) ;
      MP_File = 0 ;
      MP_Base = 0 ;
      MP_Pattern = 0 ;
   }
   else
   {
      PEL_Root::object()->destroy_possession( MP_Base ) ;
      MP_Base = 0 ;
      MP_Pattern = 0 ;
   }
   
   PEL_CHECK_POST( !has_pattern() ) ;
   PEL_CHECK_POST( !build_pattern() ) ;
}

//--------------------------------------------------------------------------
void
PEL_ModulePattern:: save_pattern( void )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: save_pattern" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;

   if( build_pattern() )
   {
      complete_indirection_choices( MP_Base ) ;
   }
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   for( size_t i=0 ; i<com->nb_ranks() ; i++ )
   {
      if( i==com->rank() )
      {
         if( i>0 ) 
         {
            PEL_Module const* prec = PEL_Module::create(0,"MAIN",pattern_filename) ;
            merge(MP_Base,prec->module("Base")) ;
            prec->destroy() ;
         }
         
         std::ofstream out( pattern_filename.c_str() ) ;
         out.close() ;
         MP_Base->write( pattern_filename, "text" ) ;
      }
      com->barrier() ;
   }
   
   BUILD = false ;
   HAS = false ;

   PEL_CHECK_POST( !has_pattern() ) ;
   PEL_CHECK_POST( !build_pattern() ) ;
}

//--------------------------------------------------------------------------
PEL_ModulePattern*
PEL_ModulePattern:: create_pattern_description(
                       PEL_Object* a_owner, std::string const& class_name ) 
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: create_pattern_description" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !class_name.empty() ) ;
   
   PEL_ModulePattern * result = 0 ;
   if( MP_Base->has_module( class_name ) )
   {
      result = new PEL_ModulePattern( a_owner, MP_Base->module( class_name ) ) ;
   }

   PEL_CHECK_POST( IMPLIES( result!=0, result->owner()==a_owner ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
PEL_ModulePattern*
PEL_ModulePattern:: create( PEL_Object* a_owner,
                            PEL_Module const* way,
                            std::string const& a_name )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: create" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   
   PEL_ModulePattern* root_pattern = new PEL_ModulePattern( 0, MP_Base ) ;
   if( build_pattern() )
   {
      root_pattern->add_pattern( a_name, way, mandatory ) ;
   }
   PEL_Module* a_father ;
   PEL_Module* pat = root_pattern->find_subpattern( way, a_name, a_father ) ;
   PEL_ModulePattern * result = 0 ;
   if( pat!=0 ) result = new PEL_ModulePattern( a_owner, pat, a_father ) ;
   if( build_pattern() ) PEL_ASSERT( result!=0 ) ;
   root_pattern->destroy() ; root_pattern = 0 ;

   PEL_CHECK_POST( IMPLIES( build_pattern(), result != 0 ) ) ;
   PEL_CHECK_POST( IMPLIES( result != 0, result->owner() == a_owner ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
PEL_ModulePattern*
PEL_ModulePattern:: create_clone( PEL_Object* a_owner ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: create_clone" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModulePattern * result = new PEL_ModulePattern( a_owner, mod ) ;
   
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_ModulePattern:: PEL_ModulePattern( PEL_Object* a_owner,
                                       PEL_Module* mm,
                                       PEL_Module const* a_father )
//---------------------------------------------------------------------------
    : PEL_Object( a_owner )
    , mod( mm )
    , father( a_father )
{
   PEL_LABEL( "PEL_ModulePattern:: PEL_ModulePattern" ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//--------------------------------------------------------------------------
PEL_ModulePattern*
PEL_ModulePattern:: create_subpattern( PEL_Object* a_owner,
                                       PEL_Module const* way,
                                       std::string const& path_and_name ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: create_subpattern" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   PEL_CHECK_PRE( ( way!=0 && way->has_module( path_and_name ) ) ||
                  ( allowed_module( path_and_name, way )!=0 ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModulePattern* result = 0 ;
   if( way!=0 )
   {
      string first, mn ;
      split( path_and_name, first, mn ) ;
      if( first.empty() ) first = mn ;
      way = way->module( first ) ;
      if( build_pattern() )
      {
         add_pattern( path_and_name, way, mandatory ) ;
      }
      PEL_Module* a_father ;
      PEL_Module* m = find_subpattern( way, path_and_name, a_father ) ;
      if( m!=0 ) result = 
         new PEL_ModulePattern( a_owner, m, a_father ) ;
   }
   else
   {
      PEL_ASSERT( child( mod, path_and_name ) != 0 ) ;
      result = new PEL_ModulePattern( a_owner, child( mod, path_and_name ), mod ) ;
   }
   
   
   PEL_CHECK_POST( IMPLIES( build_pattern(), (result != 0) ) ) ;
   PEL_CHECK_POST( IMPLIES( result != 0, result->owner() == a_owner ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_ModulePattern:: ~PEL_ModulePattern( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
bool 
PEL_ModulePattern:: is_mutable( PEL_Module const* module )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: is_mutable" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( module != 0 ) ;

   bool result = false ;
   PEL_Module const* start = module ;
   std::string const& indirection = indirection_to( module ) ;
   if( !indirection.empty() ) 
   {
      if( MP_Pattern->has_module( indirection ) )
      {
         start = MP_Pattern->module( indirection ) ;
      }
      else
      {
         start = 0 ;
      }
   }
   if( start!= 0 )
   {
      PEL_ModuleIterator * pel_it = start->create_module_iterator( 0 ) ;
      for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
      {
         PEL_Module* sub_mod = pel_it->item() ;
         if( is_condition_description( sub_mod ) )
         {
            result = true ;
            break ;
         }
      }
      pel_it->destroy() ;
   }
   
   return( result ) ;
}

//---------------------------------------------------------------------------
std::string const&
PEL_ModulePattern:: indirection_to( PEL_Module const* module )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: indirection_to" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( module != 0 ) ;
   
   return( property( module, indirection_to_keyword ) ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ModulePattern:: is_module_description( PEL_Module const* module ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: is_module_description" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( module != 0 ) ;
   
   return( property( module, type_keyword ) == module_type_keyword ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ModulePattern:: is_entry_description( PEL_Module const* module ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: is_entry_description" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( module != 0 ) ;
   
   std::string const& prop = property( module, type_keyword ) ;
   
   bool result = prop!=variable_type_keyword &&
                 prop!=module_type_keyword &&
                 prop!=conditional_type_keyword ;
   
   return( result ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ModulePattern:: is_variable_description( PEL_Module const* module ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: is_variable_description" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( module != 0 ) ;
   
   return( property( module, type_keyword )==variable_type_keyword ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ModulePattern:: is_condition_description( PEL_Module const* module ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: is_condition_description" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( module != 0 ) ;
   
   return property( module, type_keyword )==conditional_type_keyword ;
   
}


//---------------------------------------------------------------------------
std::string const&
PEL_ModulePattern:: variable_access( PEL_Module const* module,
                                     PEL_Module const* way ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: variable_access" ) ;
   static std::string result ;
   result = "" ;
   std::string const& prop = property( module, access_keyword ) ;
   // Is property an expression
   if( !access_name().has( prop ) )
   {
      if( way!=0 )
      {
         PEL_Data const* eval = way->create_evaluation( 0, prop, 0 ) ;
         if( eval!=0 )
         {
            if( eval->data_type()==PEL_Data::String )
            {
               result = eval->to_string() ;
            }
            eval->destroy() ; eval = 0 ;
         }
      }
   }
   else
   {
      result = prop ;
   }
   if( !access_name().has( result ) )
   {
      std::ostringstream os ;
      os << "Bad entry for access property of " << module->name() << " : " << result 
         << std::endl << " Valid ones are : " << access_name() << std::endl ;
      
      PEL_Error::object()->raise_internal( os.str() ) ;
   }
   
   return result ;
}


//---------------------------------------------------------------------------
stringVector const& 
PEL_ModulePattern:: mandatory_modules( PEL_Module const* way,
                                       bool first ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: mandatory_modules" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   static stringVector result(0) ;
   if( first ) result.re_initialize( 0 ) ;
   if( is_mutable(mod) && way!=0 )
   {
      PEL_List* list = list_of_selected_conditional_modules( 0, mod, way ) ;
      PEL_Iterator * it = list->create_iterator(list) ;
      for(it->start();it->is_valid();it->go_next())
      {
         PEL_Module* indirected = const_cast<PEL_Module*>( static_cast<PEL_Module const*>( it->item() )) ;
         PEL_ModulePattern const* pat = new PEL_ModulePattern(it,indirected) ;
         pat->mandatory_modules(  way, false ) ;
      }
      list->destroy() ;
   }   
   PEL_ModuleIterator * pel_it = mod->create_module_iterator( 0 ) ;
   for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
   {
      PEL_Module const* sub_mod = pel_it->item() ;
      if( is_module_description( sub_mod ) )
      {
         if( variable_access( sub_mod, way ) == access_name()( mandatory ) )
            result.extend( sub_mod->name() ) ;
      }
   }
   pel_it->destroy() ;
   result.sort() ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_Module*
PEL_ModulePattern:: allowed_module( std::string const& a_name,
                                    PEL_Module const* way ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: allowed_module" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Module* result = 0 ;
   PEL_Module* a_child = child( mod, a_name ) ;
   
   if( a_child!=0 && is_module_description( a_child ) )
   {
      result = a_child ;
   }
   else if( ( result = generic_module(way) ) != 0 )
   {
   }
   else if( is_mutable(mod) && way!=0 )
   {
      PEL_List* list = list_of_selected_conditional_modules( 0, mod, way ) ;   
      
      PEL_Iterator * it = list->create_iterator(list) ;
      for(it->start();result==0 && it->is_valid();it->go_next())
      {
         PEL_Module* indirected = const_cast<PEL_Module*>( static_cast<PEL_Module const*>( it->item() )) ;
         PEL_ModulePattern const* pat = new PEL_ModulePattern(it,indirected) ;
         result = pat->allowed_module( a_name, way ) ;         
      }
      list->destroy() ;
   }   
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_Module*
PEL_ModulePattern:: generic_module( PEL_Module const* way ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: generic_module" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Module* result = 0 ;
   PEL_ModuleIterator * pel_it = mod->create_module_iterator( 0 ) ;
   for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
   {
      PEL_Module* sub_mod = pel_it->item() ;
      if( is_module_description( sub_mod ) )
      {
         if( variable_access( sub_mod, way ) == access_name()( generic ) )
         {
            result =sub_mod ;
         }
      }
   }
   pel_it->destroy() ;
   
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_Module*
PEL_ModulePattern:: generic_entry( PEL_Module const* way ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: generic_entry" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Module * result = 0 ;
   PEL_ModuleIterator * pel_it = mod->create_module_iterator( 0 ) ;
   for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
   {
      PEL_Module * sub_mod = pel_it->item() ;
      if( is_entry_description( sub_mod ) )
      {
         if( variable_access( sub_mod, way ) == access_name()( generic ) )
         {
            result =sub_mod ;
         }
      }
   }
   pel_it->destroy() ;
            
   return( result ) ;
}

//---------------------------------------------------------------------------
stringVector const& 
PEL_ModulePattern:: mandatory_entries( PEL_Module const* way, bool first ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: mandatory_entries" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   static stringVector result(0) ;
   if( first ) result.re_initialize( 0 ) ;
   if( is_mutable(mod) && way!=0 )
   {
      PEL_List* list = list_of_selected_conditional_modules( 0, mod, way ) ;
      PEL_Iterator * it = list->create_iterator(list) ;
      for(it->start();it->is_valid();it->go_next())
      {
         PEL_Module* indirected = const_cast<PEL_Module*>( static_cast<PEL_Module const*>( it->item() )) ;
         PEL_ModulePattern const* pat = new PEL_ModulePattern(it,indirected) ;
         pat->mandatory_entries(  way, false ) ;
      }
      list->destroy() ;
   }   
   PEL_ModuleIterator * pel_it = mod->create_module_iterator( 0 ) ;
   for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
   {
      PEL_Module const* sub_mod = pel_it->item() ;
      if( is_entry_description( sub_mod ) )
      {
         if( variable_access( sub_mod, way ) == access_name()( mandatory ) )
         {
            std::string const& str = property( sub_mod, name_keyword ) ;
            if( !str.empty() )
            {
               result.extend( property( sub_mod, name_keyword ) ) ;
            }
            else
            {
               PEL::out() << "No property " << name_keyword << " defined for " <<
                  sub_mod->absolute_path_name() << " control module" << std::endl ;
            }            
         }
      }
   }
   pel_it->destroy() ;
   result.sort() ;
   
   return( result ) ;
}

//---------------------------------------------------------------------------
stringVector const& 
PEL_ModulePattern:: provided_variables( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: provided_variables" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   static stringVector result(0) ;
   result.re_initialize( 0 ) ;
   PEL_ModuleIterator * pel_it = mod->create_module_iterator( 0 ) ;
   for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
   {
      PEL_Module const* sub_mod = pel_it->item() ;
      if( is_variable_description( sub_mod ) )
      {
         result.extend( property( sub_mod, name_keyword ) ) ;
      }
   }
   pel_it->destroy() ;
   result.sort() ;
   
   return( result ) ;
}

//---------------------------------------------------------------------------
std::string const& 
PEL_ModulePattern:: type_of_entry( std::string const& a_name ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: type_of_entry" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( PEL_Module::dirname(a_name).empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   static string result ;
   PEL_Module* a_child = child( mod, a_name ) ;
   if( a_child!=0 && is_entry_description( a_child ) )
   {
      result = property( a_child, type_keyword ) ;
   }
   else
   {
      result = "no type defined" ;
   }
   
   return( result ) ;
}

//---------------------------------------------------------------------------
std::string const& 
PEL_ModulePattern:: name( void ) const
//---------------------------------------------------------------------------
{
   return( mod->name() ) ;
}

//---------------------------------------------------------------------------
PEL_Module const*
PEL_ModulePattern:: allowed_entry( std::string const& a_name,
                                   PEL_Module const* way ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: allowed_entry" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( PEL_Module::dirname(a_name).empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Module const* result = 0 ;
   PEL_Module const* a_child = child( mod, a_name ) ;
   
   if( a_child!=0 && is_entry_description(a_child) )
   {
      result = a_child ;
   }
   else if( ( result = generic_entry(way) ) != 0 )
   {
   }
   else if( is_mutable(mod) && way!=0 )
   {
      PEL_List* list = list_of_selected_conditional_modules( 0, mod, way ) ;
      PEL_Iterator * it = list->create_iterator(list) ;
      for(it->start();result==0 && it->is_valid();it->go_next())
      {
         PEL_Module* indirected = const_cast<PEL_Module*>( static_cast<PEL_Module const*>( it->item() )) ;
         PEL_ModulePattern const* pat = new PEL_ModulePattern(it,indirected) ;
         result = pat->allowed_entry( a_name, way ) ;
      }
      list->destroy() ;
   }   
   
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: check_allowed_value_for_entry( PEL_Module* result,
                                                   PEL_Module const* model,
                                                   std::string const& a_name,
                                                   PEL_Module const* way ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: check_allowed_value_for_entry" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( result != 0 ) ;
   PEL_CHECK( model != 0 ) ;
   PEL_CHECK( !a_name.empty() ) ;
   PEL_CHECK( PEL_Module::dirname(a_name).empty() ) ;
   PEL_CHECK( way != 0 ) ;
   PEL_CHECK( is_entry_description( model ) ) ;
   PEL_CHECK( way->has_entry( a_name ) ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   
   
   std::string const& type_of_model = property( model, type_keyword ) ;
   PEL_Data const* value = way->data_of_entry( a_name ) ;
   std::string const& type_of_data = PEL_Data::type_name( value->data_type() ) ;

   if( type_of_model.empty() )
   {
      PEL_Error::notify( result,
                         "Entry "+a_name+" is not allowed",
                         way->absolute_path_name() ) ;
   }
   else if( type_of_model != type_of_data )
   {
      PEL_Error::notify( result,
                         "Entry "+a_name+" is of type " + type_of_data +
                         " but should be of type "+ type_of_model,
                         way->absolute_path_name() ) ;
   }
   else
   {
      if( !property( model, test_keyword ).empty() )
      {
         std::string const& expression = property( model, test_keyword ) ;
      
         PEL_Data const* eval = way->create_evaluation( 0, expression, 0 ) ;
         if( eval==0 )
         {
            PEL_Error::notify( result,
                               "Unable to evaluate ("+expression+") for entry "+a_name,
                               way->absolute_path_name() ) ;        
         }
         else
         {
            if( eval->data_type()!=PEL_Data::Bool )
            {
               PEL_Error::notify( result,
                                  "Test (" + expression + ") for entry "+a_name+
                                  " is not a valid boolean expression",
                                  way->absolute_path_name() ) ;        
         
            }
            else if( !eval->to_bool() )
            {
               PEL_Error::notify( result,
                                  "Test for entry "+a_name+
                                  " : (" + expression + ") failed",
                                  way->absolute_path_name() ) ;         
            }
            eval->destroy() ;
         }
         
      }
      if( model->has_entry( in_keyword ) )
      {
         check_in_spec( result,
                        model,
                        a_name,
                        way ) ;
      }
      if( model->has_entry( vector_in_keyword ) )
      {
         check_vector_in_spec( result,
                               model,
                               a_name,
                               way ) ;
      }
      if( model->has_entry( unique_keyword ) )
      {
         check_unique_spec( result,
                            model,
                            a_name,
                            way ) ;
      }
      if( model->has_entry( select_in_keyword ) )
      {
         check_evaluable( true, model, select_in_keyword, PEL_Data::String ) ;
         std::string const& regexp =
            model->data_of_entry( select_in_keyword )->to_string(
                                                          model->context() ) ;
         std::string where = "true" ;
         if( model->has_entry( where_keyword ) )
         {
            check_evaluable( true, model, where_keyword, PEL_Data::String ) ;
            where =  model->data_of_entry( where_keyword )->to_string(
                                                          model->context() ) ;
         }
         PEL_List* list = way->create_data_selection( 0, regexp, 0, where ) ;
         
         stringVector in(0) ;
         PEL_ListIterator * it = list->create_iterator( list ) ;
         for( it->start() ; it->is_valid() ; it->go_next() )
         {
            std::string item = static_cast<PEL_Data*>( it->item() )->value_as_string(0) ;
            PEL::remove_enclosing_characters(item,'\"') ;
            in.append( item ) ;
         }
         list->destroy() ;
         check_evaluable( false, way, a_name, PEL_Data::Undefined ) ;
         std::string val = value->value_as_string( way->context() ) ;
         PEL::remove_enclosing_characters(val,'\"') ;
         if( !in.has( val ) )
         {
            PEL_Error::notify( result,
                               "\"" + val + "\" is not a valid value for "
                               + a_name,
                               way->absolute_path_name(),
                               in ) ;
         }
      }
      if( model->has_entry( file_keyword ) )
      {
         check_evaluable( true, model, file_keyword, PEL_Data::String ) ;
         string const& mode =
            model->data_of_entry( file_keyword )->to_string(
                                                         model->context() ) ;
         check_evaluable( false, way, a_name, PEL_Data::String ) ;
         std::string const& filename = value->to_string( way->context() ) ;
         if( ( mode=="read" && ! PEL_System::can_read( filename ) ) ||
             ( mode=="write" && !PEL_System::can_write( filename ) ) )
         {
            PEL_Error::notify( result,
                               "Entry " + a_name + " = " + filename + " : file can't be accessed with mode : " + mode,
                               way->absolute_path_name() ) ;
         }
      }
      
   }   
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: check_in_spec( PEL_Module* result,
                                   PEL_Module const* model,
                                   std::string const& a_name,
                                   PEL_Module const* way )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: check_in_spec" ) ;
   
   PEL_Data const* value = way->data_of_entry( a_name ) ;
   PEL_Data::Type dtype=value->data_type() ;
   bool ok = true ;
   if( dtype==PEL_Data::Bool || dtype==PEL_Data::BoolVector )
   {
      check_evaluable( true, model, in_keyword, PEL_Data::BoolVector ) ;
      boolVector const& in =
         model->data_of_entry( in_keyword )->to_bool_vector(
            model->context() ) ;

      check_evaluable( false, way, a_name, dtype ) ;
      
      if( dtype==PEL_Data::Bool )
      {
         bool dval = value->to_bool( way->context() ) ;
         ok = in.has( dval ) ;
      }
      else
      {
         boolVector const& dval = value->to_bool_vector( way->context() ) ;
         for( size_t i=0 ; ok && i<dval.size() ; i++ )
            ok = in.has( dval(i) ) ;
      }
   }
   else if( dtype==PEL_Data::Int || dtype==PEL_Data::IntVector )
   {
      check_evaluable( true, model, in_keyword, PEL_Data::IntVector ) ;
      intVector const& in =
         model->data_of_entry( in_keyword )->to_int_vector(
            model->context() ) ;

      check_evaluable( false, way, a_name, dtype ) ;
      
      if( dtype==PEL_Data::Int )
      {
         int dval = value->to_int( way->context() ) ;
         ok = in.has( dval ) ;
      }
      else
      {
         intVector const& dval = value->to_int_vector( way->context() ) ;
         for( size_t i=0 ; ok && i<dval.size() ; i++ )
            ok = in.has( dval(i) ) ;
      }
   }
   else if( dtype==PEL_Data::Double || dtype==PEL_Data::DoubleVector )
   {
      check_evaluable( true, model, in_keyword, PEL_Data::DoubleVector ) ;
      doubleVector const& in =
         model->data_of_entry( in_keyword )->to_double_vector(
            model->context() ) ;

      check_evaluable( false, way, a_name, dtype ) ;
      
      if( dtype==PEL_Data::Double )
      {
         double dval = value->to_double( way->context() ) ;
         ok = in.has( dval ) ;
      }
      else
      {
         doubleVector const& dval = value->to_double_vector( way->context() ) ;
         for( size_t i=0 ; ok && i<dval.size() ; i++ )
            ok = in.has( dval(i) ) ;
      }
   }
   else if( dtype==PEL_Data::String || dtype==PEL_Data::StringVector )
   {
      check_evaluable( true, model, in_keyword, PEL_Data::StringVector ) ;
      stringVector const& in =
         model->data_of_entry( in_keyword )->to_string_vector(
            model->context() ) ;

      check_evaluable( false, way, a_name, dtype ) ;
      
      if( dtype==PEL_Data::String )
      {
         string dval = value->to_string( way->context() ) ;
         ok = in.has( dval ) ;
      }
      else
      {
         stringVector const& dval = value->to_string_vector( way->context() ) ;
         for( size_t i=0 ; ok && i<dval.size() ; i++ )
            ok = in.has( dval(i) ) ;
      }
   }
   else  
   {
      PEL_Error::object()->raise_plain(
         PEL_Data::type_name( dtype ) +
         " : unsupported data type for _in specification" ) ;
   }
   
   if( !ok )
   {
      std::string val = value->value_as_string( way->context() ) ;
      string tab = model->data_of_entry( in_keyword )->value_as_string( model->context() ) ;
      stringVector inval( tab.substr(2, tab.length()-4), ' ' ) ;
      PEL_Error::notify( result,
                         val + " is not a valid value for "+ a_name,
                         way->absolute_path_name(),
                         inval ) ;
   }
   
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: check_unique_spec( PEL_Module* result,
                                       PEL_Module const* model,
                                       std::string const& a_name,
                                       PEL_Module const* way )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: check_unique_spec" ) ;

   PEL_Data const* value = way->data_of_entry( a_name ) ;
   std::string val = value->value_as_string( way->context() ) ;
   stringVector tab(val, ' ' ) ;
   bool ok = true ;

   for( size_t i=0 ; ok && i<tab.size() ; i++ )
   {
      for( size_t j=i+1 ; ok && j<tab.size() ; j++ )
      {
         ok = ok && tab(i) != tab(j) ;
      }
   }
 
   if( !ok )
   {
      PEL_Error::notify( result,
                         val + " is not a valid value for "+ a_name
                         + ": values must be uniques",
                         way->absolute_path_name() ) ;
   } 
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: check_vector_in_spec( PEL_Module* result,
                                          PEL_Module const* model,
                                          std::string const& a_name,
                                          PEL_Module const* way )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: check_vector_in_spec" ) ;
   
   PEL_Data const* value = way->data_of_entry( a_name ) ;
   PEL_Data::Type dtype=value->data_type() ;
   bool ok = true ;
   if( dtype==PEL_Data::BoolVector )
   {
      check_evaluable( true, model, vector_in_keyword, PEL_Data::BoolArray2D ) ;
      boolArray2D const& in =
         model->data_of_entry( vector_in_keyword )->to_bool_array2D(
            model->context() ) ;

      check_evaluable( false, way, a_name, dtype ) ;
      
      boolVector const& val = value->to_bool_vector( way->context() ) ;
      ok = in.index_bound(1)==val.size() ;
      if( ok )
      {
         for( size_t i=0 ; i<in.index_bound(0) ; i++ )
         {
            ok = true ;
            for( size_t j=0 ; ok && j<val.size() ; j++ )
               ok = in(i,j) == val(j) ;
            if( ok ) break ;
         }
      }     
   }
   else if( dtype==PEL_Data::IntVector )
   {
      check_evaluable( true, model, vector_in_keyword, PEL_Data::IntArray2D ) ;
      intArray2D const& in =
         model->data_of_entry( vector_in_keyword )->to_int_array2D(
            model->context() ) ;

      check_evaluable( false, way, a_name, dtype ) ;
      
      intVector const& val = value->to_int_vector( way->context() ) ;
      ok = in.index_bound(1)==val.size() ;
      if( ok )
      {
         for( size_t i=0 ; i<in.index_bound(0) ; i++ )
         {
            ok = true ;
            for( size_t j=0 ; ok && j<val.size() ; j++ )
               ok = in(i,j) == val(j) ;
            if( ok ) break ;
         }
      }     
   }
   else if( dtype==PEL_Data::DoubleVector )
   {
      check_evaluable( true, model, vector_in_keyword, PEL_Data::DoubleArray2D ) ;
      doubleArray2D const& in =
         model->data_of_entry( vector_in_keyword )->to_double_array2D(
            model->context() ) ;

      check_evaluable( false, way, a_name, dtype ) ;
      
      doubleVector const& val = value->to_double_vector( way->context() ) ;
      ok = in.index_bound(1)==val.size() ;
      if( ok )
      {
         for( size_t i=0 ; i<in.index_bound(0) ; i++ )
         {
            ok = true ;
            for( size_t j=0 ; ok && j<val.size() ; j++ )
               ok = in(i,j) == val(j) ;
            if( ok ) break ;
         }
      }     
   }
   else if( dtype==PEL_Data::StringVector )
   {
      check_evaluable( true, model, vector_in_keyword, PEL_Data::StringArray2D ) ;
      stringArray2D const& in =
         model->data_of_entry( vector_in_keyword )->to_string_array2D(
            model->context() ) ;

      check_evaluable( false, way, a_name, dtype ) ;
      
      stringVector const& val = value->to_string_vector( way->context() ) ;
      ok = in.index_bound(1)==val.size() ;
      if( ok )
      {
         for( size_t i=0 ; i<in.index_bound(0) ; i++ )
         {
            ok = true ;
            for( size_t j=0 ; ok && j<val.size() ; j++ )
               ok = in(i,j) == val(j) ;
            if( ok ) break ;
         }
      }     
   }
   else  
   {
      PEL_Error::object()->raise_plain(
         PEL_Data::type_name( dtype ) +
         " : unsupported data type for _vector_in specification" ) ;
   }
   
   if( !ok )
   {
      std::string val = value->value_as_string( way->context() ) ;
      std::string tab = model->data_of_entry( vector_in_keyword )->value_as_string( model->context() ) ;
      PEL::replace(tab,"\"","\'") ;
      
      stringVector inval( tab.substr(2, tab.length()-4), ',' ) ;
      PEL_Error::notify( result,
                         val + " is not a valid value for "+ a_name,
                         way->absolute_path_name(),
                         inval ) ;
   }
   
}

//---------------------------------------------------------------------------
PEL_Module*  
PEL_ModulePattern:: validity( PEL_Object* a_owner,
                              PEL_Module const* checked,
                              bool recurse,
                              PEL_Module* result ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: validity" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( checked!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   bool failed = false ;
   
   if( result==0 ) result = PEL_Module::create( a_owner,
                                                "validity_result" ) ;
   if( !failed )
   {
      
      PEL_ModuleIterator * pel_it = checked->create_module_iterator( 0 ) ;
      for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
      {
         PEL_Module const* sub_mod = pel_it->item() ;
         std::string const& sub_name = sub_mod->name() ;
         PEL_Module const* allowed = allowed_module( sub_name, checked ) ;
         
         if( allowed==0 )
         {
            PEL_Error::notify( result,
                               "Module "+sub_name+" is not allowed",
                               checked->absolute_path_name() ) ;        
         }
         else if( recurse )
         {
            PEL_ModulePattern * pat =
               new PEL_ModulePattern( 0, const_cast<PEL_Module*>(allowed) ) ;
            pat->validity( a_owner, sub_mod, true, result ) ;
            pat->destroy() ; pat=0 ;
         }
      }
      pel_it->destroy() ; pel_it=0 ;
      stringVector const& modules = mandatory_modules(checked) ;

      for( size_t i = 0 ; i<modules.size() ; i++ )
      {
         if( !checked->has_module( modules(i) ) )
         {
            PEL_Error::notify( result,
                               "Module "+
                               modules(i)+" is mandatory",
                               checked->absolute_path_name() ) ;        
         }
      }
      
      PEL_KeywordDataIterator * key_it = checked->create_entry_iterator( 0 ) ;
      for( key_it->start() ; key_it->is_valid() ; key_it->go_next() )
      {
         std::string const& key = key_it->item()->keyword() ;
         PEL_Module const* entry_description = allowed_entry( key, checked ) ;
         
         if( entry_description==0 )
         {
            PEL_Error::notify( result,
                          "Entry "+key+" is not allowed",
                          checked->absolute_path_name() ) ;        
         }
         else 
         {
            check_allowed_value_for_entry( result,
                                           entry_description, key, checked ) ;
         }
      }
      key_it->destroy() ; key_it=0 ;
      if( generic_entry(checked)==0 )
      {
         stringVector const& entries = mandatory_entries(checked) ;
         
         for( size_t i = 0 ; i<entries.size() ; i++ )
         {
            std::string const& key = entries(i) ;
            PEL_ASSERT( !key.empty() ) ;
            if( !checked->has_entry( key ) )
            {
               PEL_Error::notify( result,
                             "Entry "+key+" is mandatory",
                             checked->absolute_path_name() ) ;        
            }
         }
      }
   }
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_ModulePattern*
PEL_ModulePattern:: generic_pattern( PEL_Object* a_owner,
                                     PEL_Module const* way ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: generic_pattern" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( build_pattern() ) add_pattern( way->name(), way, generic ) ;
   PEL_Module* a_father ;
   PEL_Module* m = find_subpattern(way,way->name(),a_father) ;
   PEL_ModulePattern* result = 0 ;
   if( m!=0 ) result = new PEL_ModulePattern( a_owner, m, a_father ) ;
   if( build_pattern() ) PEL_ASSERT( result!=0 ) ;
   
   return( result ) ;
}

//---------------------------------------------------------------------------
void 
PEL_ModulePattern:: split( std::string const& a_name,
                           std::string& first_dir,
                           std::string& last )
//---------------------------------------------------------------------------
{
   first_dir = "" ;
   last = a_name ;
   while( last[0]=='/' ) last = last.substr( 1, last.length()-1 ) ;
   size_t i = last.find_first_of( '/' ) ;
   if( i< last.length() )
   {
      first_dir = last.substr( 0, i ) ;
      last = last.substr( i+1, last.length()-i-1 ) ;
   }
}

//---------------------------------------------------------------------------
std::string 
PEL_ModulePattern:: class_part( std::string const& class_name )
//---------------------------------------------------------------------------
{
   std::string result = class_name ;
   size_t i = result.find_first_of( "#" ) ;
   
   if( i < result.length() ) result = result.substr( 0, i ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void 
PEL_ModulePattern:: add_pattern( std::string const& path_and_name,
                                 PEL_Module const* way,
                                 Access access ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: add_pattern" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( build_pattern() ) ;
   PEL_CHECK_PRE( IMPLIES( access==mandatory, way!=0 ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   string first, mn ;
   split( path_and_name, first, mn ) ;
   if( !first.empty() )
   {
      PEL_ASSERT( way!=0 ) ;
      add_pattern( first, way, unspecified ) ;
      PEL_Module* a_father ;
      PEL_Module* m = find_subpattern(way,way->name(),a_father) ;
      PEL_ModulePattern* pat = new PEL_ModulePattern( 0, m, a_father) ;
      PEL_ASSERT( pat!=0 ) ;
      
      string f, l ;
      split( mn, f, l ) ;
      if( f.empty() ) f = l ;
      PEL_Module const* next_way = 0 ;
      if( way->has_module(f) ) next_way = way->module( f ) ;
      
      pat->add_pattern( mn, next_way, access ) ;
      pat->destroy() ;
   }
   else
   {
      add_pattern_simple( mn, way, access ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}
   
//---------------------------------------------------------------------------
bool 
PEL_ModulePattern:: is_module_matching( std::string const& module_name,
                                        std::string const& pattern_name ) 
//---------------------------------------------------------------------------
{
   return( module_name == pattern_name ) ;
}

//---------------------------------------------------------------------------
void 
PEL_ModulePattern:: check_name_validity( std::string const& a_name )
//---------------------------------------------------------------------------
{
   bool result = a_name.find( "-" ) > a_name.length() &&
                 a_name.find( " " ) > a_name.length() ;
   if( !result )
   {
      PEL_Error::object()->raise_plain(
         a_name +
         " is not a valid entry name, it contains white space or special caracter " ) ;
   }
   
}

//---------------------------------------------------------------------------
PEL_Module*
PEL_ModulePattern:: module_of_pattern( std::string const& class_name,
                                       std::string const& instance_name,
                                       bool build_if_not_exist ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern::module_of_pattern" ) ;
   PEL_CHECK( has_pattern() ) ;
   
   PEL_Module* result = 0 ;
   check_name_validity( instance_name ) ;
   
   if( build_if_not_exist )
   {
      if( !MP_Pattern->has_module( class_name ) )
      {
         MP_Pattern->add_module( PEL_Module::create( MP_Pattern, class_name ) ) ;
      }
      PEL_Module* classmod = MP_Pattern->module( class_name ) ;
      if( !classmod->has_module( instance_name ) )
      {
         classmod->add_module( PEL_Module::create( classmod, instance_name ) ) ;
      }
   }
   if( MP_Pattern->has_module( class_name ) &&
       MP_Pattern->module( class_name )->has_module( instance_name ) )
   {
      result = MP_Pattern->module( class_name )->module( instance_name ) ;
   }
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: set_property( PEL_Module* mod,
                                  std::string const& name,
                                  std::string const& property ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: set_property" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( !name.empty() ) ;
   PEL_CHECK( PEL_Module::dirname(name).empty() ) ;
   PEL_CHECK( !property.empty() ) ;
   
   if( !mod->has_entry( name ) )
   {
      mod->add_entry( name, PEL_String::create( mod, property ) ) ;
   }
   else
   {
      check_evaluable( true, mod, name, PEL_Data::String ) ;
      if( !( mod->data_of_entry( name )->to_string() == property ) )
      {
         PEL_Error::object()->raise_internal(
            "Unconsistent module description for module " +
            mod->name() + "\n for which status of " + name
            + " was declared to be " + mod->data_of_entry( name )->to_string()
            + " and is now declared to be " + property ) ;
      }
   }
}

//---------------------------------------------------------------------------
std::string const& 
PEL_ModulePattern:: property( PEL_Module const* mod, std::string const& name ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: property" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( mod!=0 ) ;
   PEL_CHECK( !name.empty() ) ;
   PEL_CHECK( PEL_Module::dirname(name).empty() ) ;
   
   static std::string const result_null = "" ;

   if( mod->has_entry( name ) ) check_evaluable( true, mod, name, PEL_Data::String ) ;
   std::string const& result =
      mod->has_entry( name ) ?
         mod->data_of_entry( name )->to_string( mod->context() ) :
         result_null ;

   PEL_CHECK_POST( EQUIVALENT( mod->has_entry( name ), !result.empty() ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void 
PEL_ModulePattern:: add_pattern_simple( std::string const& a_name,
                                        PEL_Module const* way,
                                        Access access ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: add_pattern_simple" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( IMPLIES( access==mandatory, way!=0 ) ) ;
   PEL_CHECK( !a_name.empty() ) ;
   PEL_CHECK( PEL_Module::dirname(a_name).empty() ) ;
   PEL_CHECK( IMPLIES( way!=0 ,
                       is_module_matching( class_part(way->name() ),
                                           class_part(a_name) ) ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   std::string mod_name = a_name ;
   if(access==generic) // We searching for a "good" name
   {
      mod_name = inferred_generic_name( a_name ) ;      
   }
   PEL_Module * n = allowed_module( mod_name, 0 ) ;
   
   if( n==0 )
//   if( !mod->has_module( a_name ) )
   {
      // generic_module_keyword Can't be added with other modules
      std::string bad_mixed = "" ;
      
      if( access==generic  )
      {
         PEL_ModuleIterator * pel_it = mod->create_module_iterator( 0 ) ;
         for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
            if( is_module_description( pel_it->item() ) &&
                variable_access( pel_it->item() ,way)!=access_name()(generic) )
            {
               bad_mixed += " " + pel_it->item()->name() ;
            }
         pel_it->destroy() ;
      }
      else
      {
         if( generic_module(way)!=0 )
            bad_mixed = a_name ;
      }
      if( !bad_mixed.empty() )
      {
         PEL_Error::object()->raise_plain(
            "When building pattern " + mod->absolute_path_name() +
            "\n from module " + way->absolute_path_name() +
            "\n we can't mix generic module construction \n" +
            " with non-generic one : " + bad_mixed ) ;
      }
      
      n = PEL_Module::create( mod, mod_name )  ;
      mod->add_module( n ) ;
   }

   set_property( n, type_keyword, module_type_keyword ) ;
   if( access!=generic || !n->has_entry( name_keyword ) )
      set_property( n, name_keyword, mod_name ) ;
   
   if( way!=0 )
   {
      std::string on ="" ;
      std::string to ="" ;
      std::string val ="" ;
      
      if( way->has_entry( "concrete_name" ) || way->has_entry( "type" ) )
      {
         if(  way->has_entry( "concrete_name" ) )
         {
            to = class_part( way->name() ) ;
            on = "concrete_name" ;
         }
         else
         {
            on = "type" ;
         }
         check_evaluable( false, way, on, PEL_Data::String ) ;
         val = way->data_of_entry( on )->to_string( way->context() ) ;
      }
      if( !on.empty() )
      {
         PEL_ASSERT( !val.empty() ) ;
         check_name_validity( val ) ;

         if( !n->has_module(on) )
         {
            PEL_ModulePattern* pat = new PEL_ModulePattern( 0, n ) ;
            pat->add_entry_simple( on, PEL_Data::String, mandatory ) ;
            pat->destroy() ;
         }
         
         PEL_Module* on_mod =  n->module(on) ;
         extend_entries( on_mod, in_keyword, val ) ;
         
         if( !to.empty() )
         {
            set_property( n, indirection_to_keyword, to ) ;
         }
         
         PEL_ModulePattern* pat = new PEL_ModulePattern(0,n) ;
         pat->add_entry_simple( on, PEL_Data::String, mandatory ) ;
         pat->destroy() ;
         
         PEL_Module* indirected = 0 ;
         if( to.empty() )
         {
            if( !n->has_module( val ) )
               n->add_module( PEL_Module::create( n, val ) ) ;
            indirected = n->module(val) ;
         }
         else
         {
            indirected = module_of_pattern( to,
                                            val,
                                            true ) ;
         }
         set_property( indirected, type_keyword, conditional_type_keyword ) ;
         set_property( indirected, condition_keyword, "#("+on+")=\""+ val +"\"") ;

      }
   }
   
      
   std::string acc_str = access_name()( access ) ;
      
   if( access!=unspecified )
   {
      std::string const& prop = property( n, access_keyword ) ;
   
      if( prop.empty() )
      {
         set_property( n, access_keyword, access_name()( access ) ) ;
      }
      else if( access == optional && prop==access_name()(mandatory) )
      {
         PEL_Error::object()->raise_plain(
            "When building pattern, module " + a_name + " is accessed \n" +
            " as both a mandatory and an optional module" ) ;
      }
   }
   
   PEL_CHECK_INV( invariant() ) ;
}



//---------------------------------------------------------------------------
std::string
PEL_ModulePattern:: inferred_generic_name( std::string const& a_name ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: inferred_generic_name" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( !a_name.empty() ) ;
   
   std::string mod_name = class_part( a_name ) ;
   if( class_part( a_name ) == a_name )
   {
      std::string const& f = mod->name() ;
      if( f.find( "list_of_" ) == 0 && f.length()>8 ) 
      {
         mod_name = f.substr( 8 ) ;
            
      }
      else if( f.rfind( "s" )==f.length()-1 )
      {
         mod_name = f ;
      }

      if( mod_name.length()>1 &&
          mod_name.rfind( "s" )==mod_name.length()-1 )
         mod_name = mod_name.substr( 0, mod_name.length()-1 ) ;
   }
   return( mod_name ) ;
}



//---------------------------------------------------------------------------
PEL_Module* 
PEL_ModulePattern:: find_subpattern( PEL_Module const* way,
                                     std::string const& a_name,
                                     PEL_Module *& a_father ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: find_subpattern" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( !a_name.empty() ) ;
   PEL_CHECK( way!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* result = mod ;
   string first, mn ;
   split( a_name, first, mn ) ;
   if( !first.empty() )
   {
      PEL_ASSERT( way!=0 ) ;
      PEL_Module* step = find_subpattern( way, first, a_father) ;
      PEL_ModulePattern* pat = new PEL_ModulePattern( 0, step ) ;
      string f, l ;
      split( mn, f, l ) ;
      if( f.empty() ) f = l ;
      result = pat->find_subpattern( way->module( f ), mn, a_father ) ;
      pat->destroy() ;
   }
   else
   {
      result = find_subpattern_simple( way, mn, a_father ) ;
   }
   
   
   PEL_CHECK_INV( invariant() ) ;
   return result ;
}

//---------------------------------------------------------------------------
PEL_Module* 
PEL_ModulePattern:: find_subpattern_simple( PEL_Module const* way,
                                            std::string const& a_name,
                                            PEL_Module *& a_father ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: find_subpattern_simple" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( way!=0 ) ;
   PEL_CHECK( !a_name.empty() ) ;
   PEL_CHECK( PEL_Module::dirname(a_name).empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Module* result = 0 ;
   PEL_Module* a_child = child( mod, a_name ) ;
   if( is_module_matching( way->name(), a_name ) &&
       a_child!=0 &&
       is_module_description(a_child) )
   {
      result = a_child ;
   }
   if( result==0 ) result = generic_module(way) ;
   a_father = result ;
   
   if( result!=0 && is_mutable( result ) )
   {
      PEL_List * list = list_of_selected_conditional_modules( 0, result, way ) ;
      if( list->count() == 0 ) 
      {
         std::string mess = "Reading module " + way->absolute_path_name() ;
         mess += "\n   comparing to pattern " + result->absolute_path_name() ;
         mess += "\n Indirection is not allowed" ;
         PEL_Error::object()->raise_plain( mess ) ;
      }
      else if( list->count() != 1 ) 
      {      
         std::string mess = "Reading module " + way->absolute_path_name() ;
         mess += "\n   comparing to pattern " + result->absolute_path_name() ;
         mess += "\nMore than only one indirection is not allowed in pattern: " ;
         PEL_Iterator * it = list->create_iterator(0) ;
         for(it->start();it->is_valid();it->go_next())
            mess += " " + static_cast<PEL_Module const*>( it->item() )->name() ;
         it->destroy() ;
         PEL_Error::object()->raise_plain( mess ) ;
      }
      result =  static_cast<PEL_Module *>( list->at(0) ) ;
      list->destroy() ;
   }
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}



//---------------------------------------------------------------------------
PEL_List* 
PEL_ModulePattern:: list_of_selected_conditional_modules(
                           PEL_Object * a_owner,
                           PEL_Module const* module, PEL_Module const* way ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: list_of_selected_conditional_modules" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( module != 0 ) ;
   PEL_CHECK( is_mutable( module ) ) ;
   
   PEL_List* result = PEL_List::create( a_owner ) ;
   PEL_List* list = list_of_conditional_modules(0,module) ;
   
   PEL_Iterator * pel_it = list->create_iterator( list ) ;
   for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
   {
      PEL_Module* cond_mod = static_cast<PEL_Module*>(pel_it->item()) ;
      PEL_CHECK( is_condition_description( cond_mod ) ) ;
      std::string expression = property( cond_mod, condition_keyword ) ;

      // Recognize old style assertions
      size_t idx ;
      if( ( idx = expression.find( "==" ) ) < expression.length() ) 
      {
         std::ostringstream os ;
         
         os << "In controler, old style conditional expression "<<expression ;
         std::istringstream is(expression) ;
         std::string key, eq, res ;
         is >> key >> eq >> res ;
         expression = "#(" + key + ") = \"" + res + "\"" ;
         os << " should be modified to " << expression << std::endl ;
         PEL_Error::object()->raise_plain( os.str() ) ;
         
      }
      // PEL::out()<<"Evaluation de : "<<expression<<std::endl ;
      PEL_Data const* eval = 0 ;
      try 
      {
         eval = way->create_evaluation( 0, expression, 0 ) ;
      }
      catch( PEL_Exceptions::Error e )
      {
         PEL::out() << "Unable to interpret condition " << expression <<
            " in controler at " << cond_mod->absolute_path_name() << " with context : " << std::endl ;
         way->print( PEL::out(), 0 ) ;
         
         eval = 0 ;        
         throw e ;
      }
      if( eval!=0 )
      {
         if( eval->data_type()!=PEL_Data::Bool )
         {
            PEL_Error::object()->raise_plain( 
               "Test for conditional module "+cond_mod->name()+ " : "+expression+
               "\n is not a valid boolean expression when checking "+
               way->absolute_path_name() ) ;        
            
         }
         else if( eval->to_bool() )
         {
            result->append( cond_mod ) ;
         }
         eval->destroy() ;
      }
   }
   list->destroy() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}


//---------------------------------------------------------------------------
PEL_List* 
PEL_ModulePattern:: list_of_conditional_modules( PEL_Object * a_owner,
                                                 PEL_Module const* module ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: list_of_conditional_modules" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( module != 0 ) ;
   PEL_CHECK( is_mutable( module ) ) ;

   PEL_List*  result = PEL_List::create( a_owner ) ;
   PEL_Module const* start = module ;
   std::string const& indirection = indirection_to( module ) ;
   if( !indirection.empty() )
   {
      if( MP_Pattern->has_module( indirection ) )
      {
         start = MP_Pattern->module( indirection ) ;
      }
      else
      {
         start = 0 ;
      }
   }
   if( start!= 0 )
   {
      PEL_ModuleIterator * pel_it = start->create_module_iterator( 0 ) ;
      for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
      {
         PEL_Module* sub_mod = pel_it->item() ;
         if( is_condition_description( sub_mod ) )
         {
            result->append( sub_mod ) ;
         }
      }
      pel_it->destroy() ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: add_entry( std::string const& a_keyword,
                               PEL_Module const* way,
                               PEL_Data::Type type,
                               Access acc )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: add_entry" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !a_keyword.empty() ) ;
   
   PEL_CHECK_PRE( IMPLIES( !PEL_Module::dirname(a_keyword).empty(),
                           way!=0 ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   if( a_keyword!="type" && a_keyword!="concrete_name" )
   {
      string dir = PEL_Module::dirname(a_keyword) ;
      if( !dir.empty() )
      {
         PEL_ASSERT( way!=0 ) ;
         if( !way->has_module( dir ) )
         {
            PEL_Error::object()->raise_plain(
               "Can't access to keyword " + a_keyword + " from module " + way->name() +
               " when path to this keyword doesn't exist " ) ;
         }
         std::string first, last ;
         split( dir, first, last ) ;
         if( first.empty() ) first = last ;
      
         add_pattern( dir, way->module( first ), unspecified ) ;
         PEL_ModulePattern* pat = create_subpattern( 0, way, dir ) ;
         pat->add_entry_simple( PEL_Module::basename(a_keyword), type, acc ) ;
         pat->destroy() ;
      }
      else
      {
         add_entry_simple( a_keyword, type, acc ) ;
      }
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: add_entry_simple( std::string const& a_keyword,
                                      PEL_Data::Type type,
                                      Access access )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: add_entry_simple" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( !a_keyword.empty() ) ;
   PEL_CHECK( PEL_Module::dirname(a_keyword).empty() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* mod_entry = child( mod, a_keyword ) ;
   if( mod_entry==0 )
   {
      mod_entry = PEL_Module::create( mod, a_keyword ) ;      
      mod->add_module( mod_entry ) ;
   }
   
   set_property( mod_entry, name_keyword, a_keyword ) ;
   
   if( type!=PEL_Data::Undefined )
   {
      set_property( mod_entry, type_keyword, PEL_Data::type_name( type ) ) ;
   }
   else
   {
      PEL_ASSERT( access==optional )  ;
   }
   
   std::string const& prop = property( mod_entry, access_keyword ) ;
   
   if( prop.empty() )
   {
      set_property( mod_entry, access_keyword, access_name()(access) ) ;
   }
   else if( access == optional && prop!=access_name()(optional) )
   {
      PEL_Error::object()->raise_plain(
         "When building pattern, data " + a_keyword + " is accessed \n" +
         " as both a "+prop+" and an optional entry" ) ;
   }  
   PEL_CHECK_INV( invariant() ) ;  
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: add_generic_keyword( std::string const& a_keyword,
                                         PEL_Module const* way,
                                         PEL_Data::Type type )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: add_generic_keyword" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !a_keyword.empty() ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   
   add_entry( "item", way, type, generic ) ; 
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: add_provided_variable( std::string const& a_keyword,
                                           PEL_Module const* way,
                                           PEL_Data::Type type )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: add_provided_variable" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !a_keyword.empty() ) ;
   PEL_CHECK_PRE( PEL_Module::dirname(a_keyword).empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Module* mod_entry = child( mod, a_keyword ) ;
   if( child( mod, a_keyword )==0 )
   {
      mod_entry = PEL_Module::create( mod, a_keyword ) ;
      mod->add_module( mod_entry ) ;
   }
   
   set_property( mod_entry, name_keyword, a_keyword ) ;
   set_property( mod_entry, type_keyword, variable_type_keyword ) ;
}
 
//--------------------------------------------------------------------------
void
PEL_ModulePattern:: expand_then_simplify( void )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: expand_then_simplify" ) ;
   PEL_CHECK( has_pattern() ) ;
   
   PEL_CHECK( MP_Base!=0 ) ;
   PEL_CHECK( MP_Pattern !=0 ) ;
   PEL_CHECK( MP_Base->has_module( "Pattern" ) ) ;
   
   MP_Base->remove_module( "Pattern" ) ;
   
   PEL_ModulePattern * pat = new PEL_ModulePattern( 0, MP_Base ) ;
   pat->expand_then_simplify_one() ;
   pat->destroy() ; pat = 0 ;

   BUILD = true ;

   PEL_CHECK_POST( build_pattern() ) ;
}

//--------------------------------------------------------------------------
void
PEL_ModulePattern:: attach_help_data( std::string const& keyword,
                                      std::string const& help )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: attach_help_data" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;

   PEL_Module* a_child = child( mod, keyword ) ;
   
   if( a_child==0 || !is_entry_description( a_child ) )
   {
      PEL_Error::object()->raise_plain(
         "When building pattern, you must first read data (" + keyword +
         ")\n then you can attach an help message to it." ) ;
   }
   set_property( a_child, help_keyword, help ) ;
}

//--------------------------------------------------------------------------
void
PEL_ModulePattern:: attach_file_extension( std::string const& filename,
                                           std::string const& mode )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: attach_file_extension" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !filename.empty() ) ;
   PEL_CHECK_PRE( mode == "read" || mode == "write" ) ;

   PEL_Module* a_child = child( mod, filename ) ;
   if( a_child==0 || !is_entry_description( a_child ) )
   {
      PEL_Error::object()->raise_plain(
         "When building pattern, you must first read data (" + filename +
         ")\n then you can attach an help message to it." ) ;
   }
   set_property( a_child, file_keyword, mode ) ;
}

//--------------------------------------------------------------------------
void
PEL_ModulePattern:: attach_default_data( std::string const& keyword,
                                         std::string const& value )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: attach_default_data" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !keyword.empty() ) ;
   PEL_CHECK_PRE( !value.empty() ) ;

   PEL_Module* a_child = child( mod, keyword ) ;
   
   if( a_child==0 || !is_entry_description( a_child ) )
   {
      if( ( keyword=="concrete_name" || keyword=="type"  ) )
      {            
         PEL_ASSERT( father!=0 ) ;
         if( father->has_module( keyword ) ) a_child = father->module( keyword ) ;
      }
      
      if( a_child==0 || !is_entry_description( a_child ) )
      {
         
         PEL_Error::object()->raise_plain(
            "When building pattern, you must first read data (" + keyword +
            ")\n then you can attach a default value to it." ) ;
      }
      
   }

   if( a_child==0 || !is_entry_description( a_child ) )
   {
      
      mod->print(PEL::out(),0) ;
      
      PEL_Error::object()->raise_plain(
         "When building pattern, you must first read data (" + keyword +
         ")\n then you can attach a default value to it." ) ;
   }
   if( ! a_child->has_entry( default_keyword ) )
   {
      std::string const& a_type = property( a_child, type_keyword ) ;
      if( a_type == "String" )
      {
         a_child->add_entry( default_keyword, PEL_String::create( a_child, value ) ) ;
      }
      else if( a_type == "Double" )
      {
         std::istringstream is( value ) ;
         double v ;
         is >> v ;
         a_child->add_entry( default_keyword, PEL_Double::create( a_child, v ) ) ;
      }
      else if( a_type == "Int" )
      {
         std::istringstream is( value ) ;
         int v ;
         is >> v ;
         a_child->add_entry( default_keyword, PEL_Int::create( a_child, v ) ) ;
      }
      else if( a_type == "Bool" )
      {
         bool v = value == "true" ;
         a_child->add_entry( default_keyword, PEL_Bool::create( a_child, v ) ) ;
      }
      else
      {
         PEL_Error::object()->raise_plain(
            "When setting default value for data (" + keyword +
            ")\n type " + a_type+ " is not supported." ) ;
      }
   }
}

//--------------------------------------------------------------------------
void
PEL_ModulePattern:: attach_verify_data( std::string const& keyword,
                                        std::string const& expression )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: attach_verify_data" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !keyword.empty() ) ;
   PEL_CHECK_PRE( !expression.empty() ) ;
   
   PEL_Module* a_child = child( mod, keyword ) ;

   if( a_child==0 || !is_entry_description( a_child ) )
   {
      PEL_Error::object()->raise_plain(
         "When building pattern, you must first read data (" + keyword +
         ")\n then you can test its value." ) ;
   }
   
   set_property( a_child, test_keyword, expression ) ;
}

//--------------------------------------------------------------------------
void
PEL_ModulePattern:: attach_list_of_valid_choices(
                                       std::string const& keyword,
                                       stringVector const& valid_choices )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: attach_list_of_valid_choices" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;
   PEL_CHECK_PRE( !keyword.empty() ) ;

   PEL_Module* a_child = child( mod, keyword ) ;
   if( a_child==0
       || !is_entry_description( a_child )
       || property( a_child, type_keyword ).empty() )
   {
      PEL_Error::object()->raise_plain(
         "When building pattern, you must first read data (" + keyword +
         ")\n then you can test its value." ) ;
   }
   
   if( !a_child->has_entry( in_keyword ) )
   {
      PEL_Data const* in_val = 0 ;
      
      std::string dtype = property( a_child, type_keyword ) ;
      if( dtype=="String" || dtype=="StringVector" )
      {
         in_val = PEL_StringVector::create( a_child, valid_choices ) ;
      }
      else if( dtype=="Double" || dtype=="DoubleVector" )
      {
         doubleVector tmp( valid_choices.size() ) ;
         for( size_t i=0 ; i<valid_choices.size() ; i++ )
         {
            std::istringstream is(valid_choices(i)) ;
            is >> tmp(i) ;
         }
         in_val = PEL_DoubleVector::create( a_child, tmp ) ;
      }
      else if( dtype=="Int" || dtype=="IntVector" )
      {
         intVector tmp( valid_choices.size() ) ;
         for( size_t i=0 ; i<valid_choices.size() ; i++ )
         {
            std::istringstream is(valid_choices(i)) ;
            is >> tmp(i) ;
         }
         in_val = PEL_IntVector::create( a_child, tmp ) ;
      }
      else if( dtype=="String" || dtype=="StringVector" )
      {
         stringVector tmp( valid_choices.size() ) ;
         for( size_t i=0 ; i<valid_choices.size() ; i++ )
         {
            std::istringstream is(valid_choices(i)) ;
            is >> tmp(i) ;
         }
         in_val = PEL_StringVector::create( a_child, tmp ) ;
      }
      else if( dtype=="Bool" || dtype=="BoolVector" )
      {
         boolVector tmp( valid_choices.size() ) ;
         for( size_t i=0 ; i<valid_choices.size() ; i++ )
         {
            std::istringstream is(valid_choices(i)) ;
            is >> tmp(i) ;
         }
         in_val = PEL_BoolVector::create( a_child, tmp ) ;
      }
      else
      {
         PEL_Error::object()->raise_plain(
            dtype + " : unsupported data type for attach_list_of_valid_choices facility" ) ;
      }
      PEL_ASSERT( in_val!=0 ) ;
      a_child->add_entry( in_keyword, in_val ) ;
   }
}

//--------------------------------------------------------------------------
void
PEL_ModulePattern:: attach_list_of_dynamic_choices(
                                       std::string const& keyword,
                                       std::string const& regexp,
                                       std::string const& where )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: attach_list_of_dynamic_choices" ) ;
   PEL_CHECK_PRE( has_pattern() ) ;

   PEL_Module* a_child = child( mod, keyword ) ;
   if( a_child==0 || !is_entry_description( a_child ) )
   {
      PEL_Error::object()->raise_plain(
         "When building pattern, you must first read data (" + keyword +
         ")\n then you can test its value." ) ;
   }
   
   if( !a_child->has_entry( select_in_keyword ) )
   {
      a_child->add_entry( select_in_keyword,
                        PEL_String::create( a_child, regexp ) ) ;
   }
   if( !a_child->has_entry( where_keyword ) )
   {
      a_child->add_entry( where_keyword,
                        PEL_String::create( a_child, where ) ) ;
   }
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: expand_then_simplify_one( void ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: expand_then_simplify_one" ) ;
   static size_t level = 0 ;
   static const size_t max_level = 100 ;
   
   level++ ;
   
   PEL_ASSERT( mod!=0 ) ;
   if( level > max_level )
   {
      PEL_Error::object()->raise_plain( "Maximum recursion level reached in expand_then_simplify_one on module of name " + mod->name() ) ;
   }
   
   if( is_mutable(mod) )
   {
      // Recover Pattern definition in main control description
      
      std::string const& to = indirection_to(mod) ;
   
      if( !to.empty() )
      {
         
         PEL_Module const* to_clone = MP_Pattern->module( to ) ;
         mod->remove_entry( indirection_to_keyword ) ;
         PEL_ModuleIterator* it = to_clone->create_module_iterator( 0 ) ;
         for( it->start() ; it->is_valid() ; it->go_next() )
            mod->add_module( it->item()->create_clone( mod ) ) ;
         it->destroy() ;
      }
      
   }   

   PEL_ModuleIterator* it = mod->create_module_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PEL_ModulePattern* pat = new PEL_ModulePattern( 0, it->item() ) ;
      pat->expand_then_simplify_one() ;
      pat->destroy() ;
   }
   it->destroy() ;

   level-- ; 
}

//--------------------------------------------------------------------------
stringVector const&
PEL_ModulePattern:: access_name( void )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: access_name" ) ;
   
   static stringVector result(0) ;
   if( result.size()==0 )
   {
      result.resize( 4 ) ;
      
      result(mandatory) = "mandatory" ;
      result(optional) = "optional" ;
      result(generic) = "generic" ;
      result(unspecified) = "" ;
   }
   return result ;
}

//---------------------------------------------------------------------------
bool
PEL_ModulePattern:: invariant( void ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( mod!=0 ) ;
   PEL_ASSERT( mod==MP_Base || mod->is_under_ownership_of( MP_Base ) ) ;
   return true ;
}

//---------------------------------------------------------------------------
PEL_Module*
PEL_ModulePattern:: child( PEL_Module const* module, std::string const& name )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: child" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( module != 0 ) ;
   PEL_CHECK( !name.empty() ) ;
   
   PEL_Module* result = 0 ;
   if( name.find( "/" ) < name.length() ) 
   {
      result = module->module( name ) ;
   }
   else
   {      
      PEL_ModuleIterator * pel_it = module->create_module_iterator( 0 ) ;
      for( pel_it->start() ; pel_it->is_valid() ; pel_it->go_next() )
      {
         PEL_Module* sub_mod = pel_it->item() ;
         if( property( sub_mod, name_keyword ) == name ) 
         {
            result = sub_mod ;
            break ;
         }
      }
      pel_it->destroy() ;
   }
   return(  result ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: extend_entries( PEL_Module* module,
                                    std::string const entrie_name,
                                    std::string const& value ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: extend_entries" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( module != 0 ) ;
   PEL_CHECK( !entrie_name.empty() ) ;
   PEL_CHECK( !value.empty() ) ;
   
   if( !module->has_entry( entrie_name ) )
   {
      module->add_entry( entrie_name,
                         PEL_StringVector::create( module, "" ) ) ;
   }
   PEL_StringVector const* vec =
      static_cast<PEL_StringVector const*>(
                           module->data_of_entry( entrie_name ) ) ;
   PEL_CHECK( dynamic_cast<PEL_StringVector const*>( vec ) !=0 ) ;
   check_evaluable( true, module, entrie_name, PEL_Data::StringVector ) ;
   stringVector tmp = vec->to_string_vector( module->context() ) ;
   tmp.extend( value ) ;
   tmp.sort() ;
   const_cast<PEL_StringVector*>( vec )->set( tmp ) ;
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: complete_indirection_choices( PEL_Module * root ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: complete_indirection_choices_stage" ) ;
   PEL_CHECK( has_pattern() ) ;
   PEL_CHECK( root!=0 ) ;
   
   if( is_mutable(root) )
   {
      std::string const& to = indirection_to(root) ;
      if( !to.empty() )
      {
         PEL_List*  list = list_of_conditional_modules(0,root) ;
         PEL_Iterator* it = list->create_iterator(list) ;
         
         for( it->start() ; it->is_valid() ; it->go_next() )
         {
            PEL_Module* item = static_cast<PEL_Module*>( it->item() ) ;
            std::string expression = property( item, condition_keyword ) ;
            // Searching for : #(key)='value'
            size_t const i1 = expression.find( "#(" ) ;
            size_t const i2 = expression.find( ")" ) ;
            size_t const i3 = expression.find( "=" ) ;
            size_t const i4 = expression.find( "\"" ) ;
            size_t const i5 = expression.find( "\"", i4+1 ) ;
            if( i1==0 && i1<i2 && i2+1==i3 && i3+1==i4 && i5<expression.length() )
            {
               std::string const key = expression.substr( i1+2, i2-i1-2 ) ;
               std::string const value = expression.substr( i4+1, i5-i4-1 ) ;
               if( root->has_module( key ) ) 
               {
                  extend_entries( root->module( key ), in_keyword, value ) ;
               }
               
            }
         }
         list->destroy() ;
      }
   }
   
   PEL_ModuleIterator* it = root->create_module_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      complete_indirection_choices( it->item() ) ;
   }
   it->destroy() ;   
}

//---------------------------------------------------------------------------
void
PEL_ModulePattern:: check_evaluable( bool pattern_mod,
                                     PEL_Module const* mod,
                                     std::string const& keyword,
                                     PEL_Data::Type data_type  )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: check_evaluable" ) ;
   PEL_CHECK( mod != 0 ) ;
   PEL_CHECK( !keyword.empty() ) ;
   std::ostringstream msg ;
   std::string s = (pattern_mod?"    ":"") ;

   if( !mod->has_entry( keyword ) )
   {
      msg << s << "In module: " << mod->absolute_path_name() << std::endl ;
      msg << s << "missing entry: " << keyword << std::endl ;      
   }
   else
   {
      PEL_Data const* d = mod->data_of_entry( keyword ) ;
      if( !d->value_can_be_evaluated( mod->context() ) )
      {
         msg << s << "In module: " << mod->absolute_path_name() << std::endl ;
         msg << s << "the data of keyword: " << keyword << std::endl ;
         msg << s << "cannot be evaluated" << std::endl ;
         stringVector const& undefined_variables =
            d->undefined_variables( mod->context() ) ;
         if( undefined_variables.size() > 0 )
         {
            msg << s << "undefined variable(s): " << std::endl ;
            for( size_t i=0 ; i<undefined_variables.size() ; ++i )
            {
               msg << s << "   - \"" << undefined_variables(i) << "\"" << std::endl ;
            }
         }
      }
      if( data_type!=PEL_Data::Undefined && d->data_type()!=data_type )
      {
         msg << s << "In module: " << mod->absolute_path_name() << std::endl ;
         msg << s << " entry: " << keyword << std::endl ;      
         msg << s << " should be of type: " << PEL_Data::type_name(data_type) << std::endl ;      
         msg << s << " but is of type: " <<  PEL_Data::type_name(d->data_type()) << std::endl ;      
      }
   }
   if( !msg.str().empty() )
   {
      if( pattern_mod )
      {
         PEL_ModulePattern_ERROR:: n0( pattern_filename, msg.str() ) ;
      }
      else
      {
         PEL_Error::object()->raise_plain( msg.str() ) ;
      }
   }

   PEL_CHECK_POST( mod->has_entry( keyword ) ) ;
   PEL_CHECK_POST(
      mod->data_of_entry( keyword )->value_can_be_evaluated(
                                                        mod->context() ) ) ;
   PEL_CHECK_POST(
      IMPLIES( data_type != PEL_Data::Undefined,
               mod->data_of_entry( keyword )->data_type() == data_type ) ) ;
}

//--------------------------------------------------------------------------
void
PEL_ModulePattern:: merge( PEL_Module* current, PEL_Module const* to_add )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ModulePattern:: merge" ) ;
   PEL_KeywordDataIterator* it_k = to_add->create_entry_iterator(0) ;
   for( it_k->start() ; it_k->is_valid() ; it_k->go_next() )
   {
      PEL_KeywordDataPair const* pair=it_k->item() ;
      std::string const& key = pair->keyword();
      
      if( current->has_entry(key) )
      {         
         PEL_Data const* newv = current->data_of_entry(key ) ;
         
         if( newv->value_as_string() != pair->data()->value_as_string() )
         {
            if( newv->data_type() != pair->data()->data_type() ||
                key!=in_keyword )
            {
               PEL_ModulePattern_ERROR::n1(key,
                                           pair->data()->value_as_string() ,
                                           newv->value_as_string() ) ;
            }
            if(newv->data_type()==PEL_Data::StringVector)
            {
               stringVector const& v_to_add =  pair->data()->to_string_vector() ;
               stringVector v_current =  newv->to_string_vector() ;
               for( size_t i=0 ; i<v_to_add.size() ; i++ )
               {
                  v_current.extend(v_to_add(i)) ;
               }
               current->replace_data_of_entry(key,
                                              PEL_StringVector::create(current,v_current)) ;
            }
            else if(newv->data_type()==PEL_Data::IntVector)
            {
               intVector const& v_to_add =  pair->data()->to_int_vector() ;
               intVector v_current =  newv->to_int_vector() ;
               for( size_t i=0 ; i<v_to_add.size() ; i++ )
               {
                  v_current.extend(v_to_add(i)) ;
               }
               current->replace_data_of_entry(key,
                                              PEL_IntVector::create(current,v_current)) ;
            }
            else if(newv->data_type()==PEL_Data::DoubleVector)
            {
               doubleVector const& v_to_add =  pair->data()->to_double_vector() ;
               doubleVector v_current =  newv->to_double_vector() ;
               for( size_t i=0 ; i<v_to_add.size() ; i++ )
               {
                  if(!v_current.has(v_to_add(i))) v_current.append(v_to_add(i)) ;
               }
               current->replace_data_of_entry(key,
                                              PEL_DoubleVector::create(current,v_current)) ;
            }
            else if(newv->data_type()==PEL_Data::BoolVector)
            {
               boolVector const& v_to_add =  pair->data()->to_bool_vector() ;
               boolVector v_current =  newv->to_bool_vector() ;
               for( size_t i=0 ; i<v_to_add.size() ; i++ )
               {
                  if(!v_current.has(v_to_add(i))) v_current.append(v_to_add(i)) ;
               }
               current->replace_data_of_entry(key,
                                              PEL_BoolVector::create(current,v_current)) ;
            }
         }
      }
      else
      {
         current->add_entry(key,pair->data()->create_clone(current)) ;
      }
   }
   it_k->destroy() ;
   
   PEL_ModuleIterator* it_m = to_add->create_module_iterator(0) ;
   for( it_m->start() ; it_m->is_valid() ; it_m->go_next() )
   {
      PEL_Module const* module=it_m->item() ;
      if( current->has_module(module->name()) )
      {
         merge(current->module(module->name()), module ) ;
      }
      else
      {
         current->add_module(module->create_clone(current)) ;
      }
   }
   it_m->destroy() ;
}

//internal---------------------------------------------------------------
void
PEL_ModulePattern_ERROR:: n0( std::string const& pattern,
                              std::string const& error )
//internal---------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** PEL_ModulePattern error\n" ;
   msg << "    bad pattern file: "
       << "\"" <<  pattern << "\"" << std::endl ;
   msg << error ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal---------------------------------------------------------------
void
PEL_ModulePattern_ERROR:: n1( std::string const& keyword,
                              std::string const& old,
                              std::string const& newv )
//internal---------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** PEL_ModulePattern error" << std::endl ;
   msg << "    problem when synchronizing pattern in parallel mode: "<< std::endl ;
   msg << "    keyword      : \"" <<  keyword << "\"" << std::endl ;
   msg << "    older value  : \"" <<  old << "\"" << std::endl ;
   msg << "    current value: \"" <<  newv << "\"" << std::endl ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
