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

#include <PEL_Module.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>

#include <PEL_BinStored.hh>
#include <PEL_Context.hh>
#include <PEL_ContextPair.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_Root.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_ModuleComparator.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_List.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>

#include <fstream>
#include <iomanip>
#include <sstream>

// Interface with Yacc
bool PEL_readFile( PEL_Module * top,
                   std::istream* input_stream,
                   std::string const& name,
                   bool debug ) ;
std::string PEL_current_parsed_module_path_name( void ) ;

//----------------------------------------------------------------------
PEL_Module*
PEL_Module:: create( PEL_Object* a_owner,
                     std::string const& a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create(a_owner,a_name)" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   PEL_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;

   PEL_Module* result = new PEL_Module( a_owner, a_name, 0 ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( result->is_empty() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_Module:: create( PEL_Object* a_owner, 
                     std::string const& a_name,
                     std::string const& file_name,
                     PEL_Context const* ct )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create(a_owner,a_name,file_name)" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   PEL_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   
   PEL_Module* result = new PEL_Module( a_owner, a_name, ct ) ;
   bool ok = PEL_readFile( result, 0, file_name, false ) ;
   
   if( !ok )
   {
      PEL_Error::object()->raise_plain(
         "*** PEL_Module error:\n"
         "    unable to complete data deck reading \""+file_name+"\"" ) ;
   }
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( !result->is_empty() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_Module:: create( PEL_Object* a_owner, 
                     std::string const& a_name,
                     std::istream& input_stream,
                     PEL_Context const* ct )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create(a_owner,a_name,input_stream)" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   PEL_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;

   PEL_Module* result = new PEL_Module( a_owner, a_name, ct ) ;
   bool ok = PEL_readFile( result, &input_stream, "", false ) ;

   if( !ok )
   {
      PEL_Error::object()->raise_plain(
         "*** PEL_Module error:\n"
         "    unable to complete data deck reading" ) ;
   }
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( !result->is_empty() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_Module:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create_clone" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* result = create_copy( a_owner ) ;
   if( FATHER != 0 )
   {
      complete_context( result, result, FATHER->context() ) ;
   }  

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   PEL_CHECK_POST( result->name() == name() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_Module:: create_as_difference( PEL_Object* a_owner,
				   std::string const& a_name,
                                   PEL_Module const* m1,
                                   PEL_Module const* m2,
				   PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create_as_difference" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   PEL_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;
   PEL_CHECK_PRE( m1 != 0 ) ;
   PEL_CHECK_PRE( m2 != 0 ) ;
   
   PEL_Module* result = PEL_Module::create( a_owner, a_name ) ;
   PEL_ModuleComparator* cmp = PEL_ModuleComparator::create( a_owner, exp ) ;
   int nb_err = cmp->compare( m1, m2, result ) ;
   result->add_entry( "nb_differences", PEL_Int::create( result, nb_err ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( result->has_entry( "nb_differences" ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module:: PEL_Module( PEL_Object* a_owner,
                         std::string const& a_name,
                         PEL_Context const* ct )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NAME( PEL_String::create( this, a_name ) )
   , FATHER( 0 )
   , MODS( PEL_List::create( this ) )
   , ENTRIES( PEL_List::create( this ) )
   , CTX( PEL_ContextSimple::create( this ) )
   , TMP_CTX( PEL_ContextPair::create( this, 0, 0) )
{
   PEL_LABEL( "PEL_Module:: PEL_Module" ) ;
   
   if( ct!=0 ) CTX->extend( ct ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Module:: ~PEL_Module( void ) 
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_Module:: modify_module_name( std::string const& a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: modify_module_name" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   PEL_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;
   PEL_CHECK_PRE( IMPLIES( a_name != name() && father() != 0,
                           !father()->has_module( a_name ) ) ) ;
   PEL_CHECK_PRE( IMPLIES( a_name != name() && father() != 0,
                           !father()->has_entry( a_name ) ) ) ;
   
   NAME->set( a_name ) ;

   PEL_CHECK_POST( name() == a_name ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: add_entry( std::string const& keyword,
                        PEL_Data const* data )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: add_entry" ) ;
   PEL_CHECK_PRE( !keyword.empty() ) ;
   PEL_CHECK_PRE( !has_module( keyword ) ) ;
   PEL_CHECK_PRE( keyword.find("/") >= keyword.size() ) ;
   PEL_CHECK_PRE( keyword.find(" ") >= keyword.size() ) ;
   PEL_CHECK_PRE( data != 0 ) ;
   PEL_CHECK_PRE( data->is_under_ownership_of(this) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( has_entry( keyword ) )
   {
      std::string mess = keyword ;
      mess += " is already the keyword of a data in Module " ;
      mess += name() ;
      PEL_Error::object()->raise_plain( mess ) ;
   }

   PEL_String* key_string = PEL_String::create( 0, keyword ) ;
   PEL_KeywordDataPair* a_entry = PEL_KeywordDataPair::create( this,
                                                               key_string,
                                                               data ) ;
   key_string->set_owner( a_entry ) ;

   if( keyword == "type" || keyword == "concrete_name" )
   {
      ENTRIES->prepend( a_entry ) ;
   }
   else
   {
      ENTRIES->append( a_entry ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( has_entry( keyword ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: add_module( PEL_Module* a_module )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: add_module" ) ;
   PEL_CHECK_PRE( a_module!=0 ) ;
   PEL_CHECK_PRE( !has_module( a_module->name() ) ) ;
   PEL_CHECK_PRE( !has_entry( a_module->name() ) ) ;
   PEL_CHECK_PRE( a_module->is_under_ownership_of(this) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   MODS->append( a_module ) ;
   a_module->FATHER = this ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( has_module( a_module->name() ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: merge_module( PEL_Module* a_module )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: merge_module" ) ;
   PEL_CHECK_PRE( a_module!=0 ) ;
   PEL_CHECK_PRE( a_module->name()==name() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   CTX->extend( a_module->CTX ) ;
   
   PEL_ModuleIterator* iteratorOtherMod =
         a_module->create_module_iterator( 0 ) ;
   for( iteratorOtherMod->start() ;
        iteratorOtherMod->is_valid() ;
        iteratorOtherMod->go_next() )
   {
      if( has_module( iteratorOtherMod->item()->name() ) )
      {
         PEL_Module* m = module( iteratorOtherMod->item()->name() ) ;
         m->merge_module( iteratorOtherMod->item() ) ;
      }
      else
      {
         add_module( iteratorOtherMod->item()->create_clone( this ) ) ;
      }
   }
   iteratorOtherMod->destroy() ;
   PEL_KeywordDataIterator* iteratorOtherAss =
         a_module->create_entry_iterator( 0 ) ;
   for( iteratorOtherAss->start() ;
        iteratorOtherAss->is_valid() ;
        iteratorOtherAss->go_next() )
   {
      if( !has_entry( iteratorOtherAss->item()->keyword() ) )
      {
         add_entry( iteratorOtherAss->item()->keyword(),
                    iteratorOtherAss->item()->data()->create_clone( this ) ) ;
      }
      else
      {
         replace_data_of_entry(
                    iteratorOtherAss->item()->keyword(),
                    iteratorOtherAss->item()->data()->create_clone( this ) ) ;
      }
   }
   iteratorOtherAss->destroy() ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: remove_module( std::string const& path_and_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: remove_module" ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   PEL_CHECK_PRE( has_module( path_and_name ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* m = module( path_and_name ) ;
   std::string dir = dirname( path_and_name ) ;
   if( dir=="" )
   {
      PEL_CHECK( MODS->has( m ) ) ;
      MODS->remove( m ) ;
   }
   else
   {
      PEL_Module* mod = 0 ;
      PEL_KeywordDataPair* assign = 0 ;
      bool ok = find( dirname( path_and_name ), mod, assign ) ;
      PEL_CHECK( ok && mod!=0 && assign==0 ) ;
      PEL_CHECK( mod->MODS->has( m ) ) ;
      mod->MODS->remove( m ) ;
   }
   if( m->FATHER!=0 )
   {
      m->CTX->extend( m->FATHER->context() ) ;
   }
   m->FATHER=0 ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !has_module( path_and_name ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: remove_entry( std::string const& path_and_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: remove_entry" ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   PEL_CHECK_PRE( has_entry( path_and_name ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* a_mod = 0 ;
   PEL_KeywordDataPair* key_data = 0 ;
   bool ok = find( path_and_name, a_mod, key_data ) ;
   PEL_CHECK( ok && a_mod==0 && key_data!=0 ) ;
   std::string dir = dirname( path_and_name ) ;
   if( dir=="" )
   {
      PEL_CHECK( ENTRIES->has( key_data ) ) ;
      ENTRIES->remove( key_data ) ;
   }
   else
   {
      PEL_Module* mod = 0 ;
      PEL_KeywordDataPair* assign = 0 ;
      ok = find( dirname( path_and_name ), mod, assign ) ;
      PEL_CHECK( ok && mod!=0 && assign==0 ) ;
      PEL_CHECK( mod->ENTRIES->has( key_data ) ) ;
      mod->ENTRIES->remove( key_data ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !has_entry( path_and_name ) ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Module:: name( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   std::string const& result = NAME->to_string() ;

   PEL_CHECK_POST( !result.empty() ) ;
   PEL_CHECK_POST( result.find( "/" ) >= result.size() ) ;
   PEL_CHECK_POST( result.find( " " ) >= result.size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Module:: absolute_path_name( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: absolute_path_name" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   static std::string result ;
   result = "" ;
   PEL_Module const* chain = this ;
   while( chain!=0 )
   {
      if( !result.empty() ) result = "/" + result ;
      
      result = chain->name() + result ;
      chain = chain->FATHER ;
   }
   result = "/" + result ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Module:: is_empty( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: is_empty" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   bool result = !has_module() && !has_entry() ;
   
   PEL_CHECK_POST( EQUIVALENT( result, !has_module() && !has_entry() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Module:: has_module( void ) const
//----------------------------------------------------------------------
{
   return( MODS->count()>0 ) ;
}

//----------------------------------------------------------------------
bool
PEL_Module:: has_module( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: has_module" ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Module* a_module = 0 ;
   PEL_KeywordDataPair* assign = 0 ;
   bool result = find( path_and_name, a_module, assign ) && a_module!=0 ;
   
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module const*
PEL_Module:: father( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: father" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Module const* result = FATHER ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( result!=0, result->module( name() )==this ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_Module:: module( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: module" ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Module* result = 0 ;
   PEL_KeywordDataPair* assign = 0 ;
   bool ok = find( path_and_name, result, assign ) && result!=0 ;
   if( !ok )
   {
      std::string mess = "Can't find module " ;
      mess += path_and_name ;
      mess += " in module " ;
      mess += name() ;
      PEL_Error::object()->raise_plain( mess ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST( result->name() == basename(path_and_name) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Module:: has_entry( void ) const
//----------------------------------------------------------------------
{
   return( ENTRIES->count()>0 ) ;
}

//----------------------------------------------------------------------
bool
PEL_Module:: has_entry( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: has_entry" ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_Module* a_module = 0 ;
   PEL_KeywordDataPair* assign = 0 ;
   bool result = find( path_and_name, a_module, assign ) && assign!=0 ;
   
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data const* 
PEL_Module:: data_of_entry( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: data_of_entry" ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Module* a_module = 0 ;
   PEL_KeywordDataPair* key_data = 0 ;
   bool ok = find( path_and_name, a_module, key_data ) && key_data!=0 ;
   if( !ok )
   {
      std::string mess = "Can't find entry " ;
      mess += path_and_name ;
      mess += " in module " ;
      mess += name() ;
      PEL_Error::object()->raise_plain( mess ) ;
   }
   PEL_CHECK( key_data->keyword() == basename(path_and_name) ) ;
   PEL_Data const* result = key_data->data() ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: replace_data_of_entry( std::string const& path_and_name,
                                    PEL_Data const* new_data ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: replace_data_of_entry" ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   PEL_CHECK_PRE( new_data != 0 ) ;
   PEL_CHECK_PRE( new_data->is_under_ownership_of( this ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* a_mod = 0 ;
   PEL_KeywordDataPair* key_data = 0 ;
   bool ok = find( path_and_name, a_mod, key_data ) && key_data!=0 ;
   if( !ok )
   {
      std::string mess = "Can't find entry " ;
      mess += path_and_name ;
      mess += " in module " ;
      mess += name() ;
      PEL_Error::object()->raise_plain( mess ) ;
   }
   PEL_CHECK( key_data->keyword() == basename(path_and_name) ) ;
   key_data->replace_data( new_data ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( data_of_entry( path_and_name ) == new_data ) ;
}

//----------------------------------------------------------------------
PEL_List*
PEL_Module:: create_data_selection( PEL_Object* a_owner,
                                    std::string const& regexp,
                                    PEL_List* result,
                                    std::string const& where ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create_data_selection" ) ;
   PEL_CHECK_PRE( !regexp.empty() ) ;
   PEL_CHECK_PRE( result==0 || result->owner()==a_owner ) ;
   
   std::string where_cp = where ;
   bool failed = substitute_variables( where_cp, 0 ) ;
   if( failed )
   {
      PEL_Error::object()->raise_plain(
         "Bad subsitution in chain "+where_cp ) ;
   }
   if(result==0)
   {
      result = PEL_List::create(a_owner) ;
      if( regexp[0]=='/' ) 
      {
         PEL_Module const* root = this ;
         while( root->FATHER!=0 ) root = root->FATHER ;
         return root->create_data_selection( a_owner, regexp, result, where_cp ) ;
      }
   }
   
   std::string token, otherToken ;
   split(regexp, token, otherToken) ;
   
   if( token==".." && FATHER!=0 )
   {
      FATHER->create_data_selection( a_owner, otherToken, result, where_cp ) ;
   }
   else if( has_module(token) )
   {
      module(token)->create_data_selection( a_owner, otherToken, result, where_cp ) ;
   }
   else if( token=="*" && !otherToken.empty() )
   {
      PEL_ModuleIterator* iteratorOtherMod = create_module_iterator( 0 ) ;
      for( iteratorOtherMod->start() ;
           iteratorOtherMod->is_valid() ;
           iteratorOtherMod->go_next() )
      {
         iteratorOtherMod->item()->create_data_selection( a_owner,
                                                          otherToken,
                                                          result, where_cp ) ;
      }
      iteratorOtherMod->destroy() ;
   }
   else if( token=="*" && otherToken.empty() ) 
   {
      PEL_KeywordDataIterator* iteratorOtherAss =
         create_entry_iterator( 0 ) ;
      for( iteratorOtherAss->start() ;
           iteratorOtherAss->is_valid() ;
           iteratorOtherAss->go_next() )
      {
         PEL_Data* clone_data =
            static_cast<PEL_Data*>
            ( iteratorOtherAss->item()->data()->create_clone( result ) ) ;
         result->append( clone_data ) ;
      }
      iteratorOtherAss->destroy() ;
   }
   else if( token=="$key" && otherToken.empty() ) 
   {
      PEL_KeywordDataIterator* iteratorOtherAss =
         create_entry_iterator( 0 ) ;
      for( iteratorOtherAss->start() ;
           iteratorOtherAss->is_valid() ;
           iteratorOtherAss->go_next() )
      {
         result->append( PEL_String::create( result,
                                             iteratorOtherAss->item()->keyword() ) ) ;
      }
      iteratorOtherAss->destroy() ;
   }
   else if( has_entry(token) && otherToken.empty() )
   {
      size_t idx ;
      while( ( idx = where_cp.find( "##(" ) ) < where_cp.length() )
      {
         where_cp.replace( idx, 3, "#(" ) ;
      }
      
      PEL_Data const* eval = create_evaluation( 0, where_cp, 0 ) ;
      if( eval->to_bool() )
      {
         result->append( PEL_DataWithContext::create( result,
                                                      data_of_entry(token),
                                                      context() ) ) ;
      }
      eval->destroy() ;
   }

   PEL_CHECK_POST( result!=0 && result->owner()==a_owner ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->index_limit() ; i++ ),
                           dynamic_cast<PEL_Data *>(result->at(i))!=0 ) ) ;
   return( result );
   
}

//----------------------------------------------------------------------
PEL_Context const*
PEL_Module:: context( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: context" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Context const* result = CTX ;
   if( FATHER!=0 )
   {
      TMP_CTX->re_initialize( FATHER->context(), CTX ) ;
      result = TMP_CTX ;
   }
   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: add_variable( PEL_Variable const* variable,
                           PEL_Data* value ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: add_variable" ) ;
   PEL_CHECK_PRE( !context()->has_variable( variable ) ) ;
   PEL_CHECK_PRE( variable->data_type()==value->data_type() ) ;
   PEL_CHECK_PRE( value->owner()==0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   value->set_owner( CTX ) ;
   CTX->extend( variable, value ) ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( context()->has_variable( variable ) ) ;
   PEL_CHECK_POST( context()->value( variable ) == value ) ;
   PEL_CHECK_POST( value->is_under_ownership_of( this ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: modify_variable( PEL_Variable const* variable,
                              PEL_Data* value ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: modify_variable" ) ;
   PEL_CHECK_PRE( context()->has_variable( variable ) ) ;
   PEL_CHECK_PRE( variable->data_type()==value->data_type() ) ;
   PEL_CHECK_PRE( value->owner()==0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   value->set_owner( CTX ) ;
   if( CTX->has_variable( variable ) )
   {
      CTX->set_value_of( variable, value ) ;
   }
   else
   {
      CTX->extend( variable, value ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( context()->has_variable( variable ) ) ;
   PEL_CHECK_POST( context()->value( variable ) == value ) ;
   PEL_CHECK_POST( value->is_under_ownership_of( this ) ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: write( std::string const& file,
                    std::string const& format ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: write" ) ;
   PEL_CHECK_PRE( !file.empty() ) ;
   PEL_CHECK_PRE( format=="text" || format=="hybrid" ) ;   
   PEL_CHECK_INV( invariant() ) ;
   
   std::ofstream out( file.c_str(), std::ios::out | std::ios::app ) ;
   if( !out )
   {
      std::string mess = "PEL_Module \"" ;
      mess += name() ;
      mess += "\" writing failure : \n   Unable to open file \"" ;
      mess += file ;
      mess += "\"" ;
      PEL_Error::object()->raise_plain( mess ) ;
   }
   bool hybrid = format=="hybrid" ;
   std::string const binary_file = file + ".bin" ;
   if( hybrid &&
       !PEL_BinStored::is_valid_binary_file( binary_file ) )
   {
      PEL_BinStored::init_binary_file( binary_file ) ;
   }
   recursive_print( out, 0, hybrid, binary_file ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: print( std::ostream& os, size_t indent_width ) const 
//----------------------------------------------------------------------
{	
   PEL_LABEL( "PEL_Module:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   recursive_print( os, (int)indent_width, false, "" ) ;
}

//-------------------------------------------------------------------------
void
PEL_Module:: display_info( std::ostream& os, size_t indent_width ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: display_info" ) ;
   PEL_Object::display_info( os, indent_width ) ;
   std::string const s( indent_width, ' ' ) ;
   os << s << "module name : " << NAME->to_string() << std::endl ;
}

//----------------------------------------------------------------------
std::string const&
PEL_Module:: current_parsed_module_path_name( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: current_parsed_module_path_name" ) ;
   static std::string result = PEL_current_parsed_module_path_name() ;
   return( result ) ;
}

//---------------------------------------------------------------------------
std::string
PEL_Module:: data_as_string( std::string const& path_and_name,
                             PEL_Context const* ct,
                             bool& failed ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: data_as_string" ) ;
   PEL_CHECK_PRE( !path_and_name.empty() ) ;
   
   if( ct==0 ) ct=context() ;
   
   std::string result ;
   
   size_t idx = path_and_name.find("/") ;
   
   if( idx < path_and_name.size() ) 
   {
      PEL_Module const* root = 0 ;
      std::string path = path_and_name.substr(idx+1) ; ;
      
      if( idx==0 ) // Absolute path
      {
         root = this ;
         
         while( root->FATHER!=0 ) root=root->FATHER ;
      }
      else if( idx==2 && path_and_name[0]=='.' && path_and_name[1]=='.' ) // ../path
      {
         root = FATHER ;
      }
      else 
      {
         std::string sub_mod = path_and_name.substr(0,idx) ;
         
         if( has_module( sub_mod ) ) 
         {
            root = module( sub_mod ) ;
         }
      }
      if( root!=0 ) 
      {         
         result = root->data_as_string( path, ct, failed ) ;
      }
      
   }
   else
   {
      failed = !has_entry( path_and_name ) ;
      if( !failed ) 
      {         
         PEL_Data const* data = data_of_entry( path_and_name ) ;
         failed = !data->value_can_be_evaluated( ct ) ;
         
         std::ostringstream str ;
         if( !failed )
         {
            result = data->value_as_string( ct ) ;
         }
      }
   }

   PEL_CHECK_POST( IMPLIES( !failed, !result.empty() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleIterator*
PEL_Module:: create_module_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create_module_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModuleIterator* result =
                           PEL_ModuleIterator::create( a_owner, MODS ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_KeywordDataIterator*
PEL_Module:: create_entry_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create_entry_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_KeywordDataIterator* result =
                  PEL_KeywordDataIterator::create( a_owner, ENTRIES ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string
PEL_Module:: basename( std::string const& path_and_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: basename" ) ;
   const char separator = '/' ;
   size_t idx = path_and_name.find_last_of( separator ) ;
   std::string result ;
   if( idx<path_and_name.length() )
   {
      result = path_and_name.substr( idx+1, path_and_name.length()-idx-1 ) ;
   }
   else
   {
      result = path_and_name ;
   }
   PEL_CHECK_POST( result.find( separator ) >= result.size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string
PEL_Module:: dirname( std::string const& path_and_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: dirname" ) ;
   const char separator = '/' ;
   size_t idx = path_and_name.find_last_of( separator ) ;
   std::string result ;
   if( idx>0 && idx<path_and_name.length() )
   {
      result = path_and_name.substr( 0, idx ) ;
   }
   else
   {
      result = "" ;
   }
   return( result ) ;
}

//---------------------------------------------------------------------------
bool
PEL_Module:: substitute_variables( std::string& replaced,
                                   PEL_Context const* ct ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: substitute_variables" ) ;
   PEL_CHECK_PRE( !replaced.empty() ) ;
   bool failed = false ;

   PEL_Context const* ctx = ( ct==0 ? context() :
                              PEL_ContextPair::create(0,context(),ct) ) ;
   // Adressage indirect des sous-modules
   size_t start = 0 ;
   size_t varidx = replaced.find("#(", start ) ;
   bool new_notation =  varidx < replaced.length() ;
   
   if( new_notation )
   {
      while( !failed &&
             varidx< replaced.length() )
      {
         if( varidx==0 || replaced[varidx-1] != '#' )
         {
            size_t endvaridx = replaced.find(")",varidx+2) ;
            if( endvaridx < replaced.length())
            {
               size_t i1 = varidx+2 ;
               size_t i2 = endvaridx-varidx-2 ;
               std::string default_value ;
            
               std::string varname =  replaced.substr( i1, i2 ) ;
               size_t i3 = varname.find( "," ) ;
               if( i3 < varname.length() ) 
               {
                  default_value = varname.substr( i3+1, varname.length()-i3 ) ;
                  varname = varname.substr( 0, i3 ) ;
               }

               std::string val = data_as_string(varname,ctx,failed) ;
               
               if( val.empty() )
               {
                  if( !default_value.empty() )
                  {
                     failed = false ;                  
                     replaced.replace( varidx, endvaridx-varidx+1, default_value ) ;
                  }
               
               }
               else
               {            
                  std::string valuestr = " " + val + " " ;
            
                  replaced.replace( varidx, endvaridx-varidx+1, valuestr ) ;
               }
            }
            else
            {
               PEL_Error::object()->raise_plain(
                  "Bad syntax in variable substitution in chain "+replaced ) ;
            }
            start = 0 ;
         }
         else
         {
            start = varidx+1 ; 
         }
         varidx = replaced.find("#(", start ) ;
      
      }
   }
   else
   {
      
      PEL_KeywordDataIterator* it = create_entry_iterator( 0 ) ;
   
      for( it->start() ; !failed && it->is_valid() ; it->go_next() )
      {
         std::string prefix ;
         PEL_KeywordDataPair const* pair = it->item() ;
         std::string a_keyword = pair->keyword() ;

         std::string valuestr ;
      
         size_t idx=0 ;
      
         while( !failed && (idx=replaced.find(a_keyword,idx))<replaced.length())
         {
            bool to_be_replaced = true ;
            if( idx > 0 )
            {
               char c = replaced[idx-1] ;
               to_be_replaced = !( ( ('a'<=c) && (c<='z') ) ||
                                   ( ('A'<=c) && (c<='Z') ) ||
                                   ( ('0'<=c) && (c<='9') ) ||
                                   ( c=='_' ) ||
                                   ( c=='"' ) ||
                                   ( c=='#' )
                  ) ; 
               
            }
            if( idx > 1 )
            {
               to_be_replaced = to_be_replaced &&
                  !( replaced[idx-2]=='#' && replaced[idx-1]=='(' ) ;
               
            }
            if( to_be_replaced && idx+a_keyword.length()<replaced.length() )
            {
               char c = replaced[idx+a_keyword.length()] ;
               to_be_replaced = !( ( ('a'<=c) && (c<='z') ) ||
                                   ( ('A'<=c) && (c<='Z') ) ||
                                   ( ('0'<=c) && (c<='9') ) ||
                                   ( c=='_' ) ||
                                   ( c=='"' ) ||
                                   ( c=='#' )
                  ) ; 
            }
         
            if( to_be_replaced )
            {
               if( valuestr.empty() )
               {
                  valuestr = " " ;
                  valuestr += data_as_string(a_keyword,ctx,failed) ;
                  valuestr += " " ;
               }
            
               if( !failed )
               {
                  replaced.replace( idx, a_keyword.length(), valuestr ) ;
               }
            }
            else idx++ ;
      
         }
         
      }
      it->destroy() ;
   }
   
   if( ct!=0 ) ctx->destroy() ;
   
   return( failed ) ;
}

//---------------------------------------------------------------------------
PEL_DataWithContext const*
PEL_Module:: create_evaluation( PEL_Object * a_owner,
                                std::string const& expression,
                                PEL_Context const* ct ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create_evaluation" ) ;
   PEL_CHECK_PRE( !expression.empty() ) ;
   
   std::ostringstream full_str ;
   full_str << "MODULE Root" << std::endl ;   
   std::string replaced = expression ;
   bool failed = substitute_variables( replaced, ct ) ;
   

   PEL_DataWithContext * result = 0 ;
   if( !failed )
   {
      full_str << "_expression_result = ( " << replaced << " )" << std::endl ;
      full_str << "END MODULE Root" << std::endl ;
   
      std::istringstream is( full_str.str() ) ;
   
      PEL_Module * built = PEL_Module::create( 0, "Test", is ) ;
      PEL_Module const* root = built->module( "Root" ) ;
      result = PEL_DataWithContext::create(
                     0,
                     root->data_of_entry( "_expression_result" ),
                     root->context() ) ;
      built->set_owner( result ) ;
      
      if( !result->value_can_be_evaluated( ct ) )
      {
         result->destroy() ; result = 0 ;
      }
      else if( a_owner != 0 )
      {
         result->set_owner( a_owner ) ;
      }
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( result != 0, result->owner() == a_owner ) ) ;
   PEL_CHECK_POST( IMPLIES( result != 0, result->value_can_be_evaluated( ct ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
PEL_Module:: comparable( PEL_Object const* other ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: comparable" ) ;
   return( PEL_Object::comparable( other ) ||
                    dynamic_cast<PEL_String const*>( other ) != 0 ) ;
}

//----------------------------------------------------------------------------
bool
PEL_Module:: is_equal( PEL_Object const* other ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: is_equal" ) ;
   // less restrictive than PEL_Object::is_equal_PRE
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;

   PEL_Module const* lex = dynamic_cast<PEL_Module const* >( other ) ;
   bool result ;
   
   if( lex!=0 )
   {
      result = NAME->is_equal( lex->NAME ) ;
   }
   else
   {
      PEL_String const* ss = dynamic_cast<PEL_String const* >( other ) ;
      PEL_ASSERT( ss != 0 ) ;
      result = NAME->is_equal( ss ) ;
   }
   PEL_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
int
PEL_Module:: three_way_comparison( PEL_Object const* other ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;
   PEL_Module const* lex = dynamic_cast<PEL_Module const* >( other ) ;
   int result ;
   
   if( lex!=0 )
   {
      result = NAME->three_way_comparison( lex->NAME ) ;
   }
   else
   {
      PEL_String const* ss = dynamic_cast<PEL_String const* >( other ) ;
      PEL_ASSERT( ss != 0 ) ;
      result = NAME->three_way_comparison( ss ) ;
   }
   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PEL_Module:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: hash_code" ) ;
   return( NAME->hash_code() ) ;
}
//----------------------------------------------------------------------
PEL_Module*
PEL_Module:: create_copy( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: create_copy" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* result = new PEL_Module( a_owner, name(), CTX ) ;

   PEL_ModuleIterator* iteratorOtherMod = create_module_iterator( 0 ) ;
   for( iteratorOtherMod->start() ;
        iteratorOtherMod->is_valid() ;
        iteratorOtherMod->go_next() )
   {
      PEL_Module* clone_mod =
                        iteratorOtherMod->item()->create_copy( result ) ;
      result->add_module( clone_mod ) ;
   }
   iteratorOtherMod->destroy() ;
   PEL_KeywordDataIterator* iteratorOtherAss = create_entry_iterator( 0 ) ;
   for( iteratorOtherAss->start() ;
        iteratorOtherAss->is_valid() ;
        iteratorOtherAss->go_next() )
   {
      PEL_Data* clone_data =
         static_cast<PEL_Data*>
            ( iteratorOtherAss->item()->data()->create_clone( result ) ) ;
      result->add_entry( iteratorOtherAss->item()->keyword(),
                         clone_data ) ;
   }
   iteratorOtherAss->destroy() ;   

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   PEL_CHECK_POST( result->name() == name() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: complete_context( PEL_Module* root,
                               PEL_Module* dup,
                               PEL_Context const* ref_ctx ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: complete_context" ) ;
   PEL_CHECK( root != 0 ) ;
   PEL_CHECK( dup != 0 ) ;
   PEL_CHECK( ref_ctx != 0 ) ;

   bool has_modif = false ;
   PEL_KeywordDataIterator* assIt = dup->create_entry_iterator( 0 ) ;
   for( assIt->start() ; assIt->is_valid() ; assIt->go_next() )
   {
      PEL_KeywordDataPair const* ass = assIt->item() ;
      PEL_Data const* dat = ass->data() ;
      complete_context( root, dat, ref_ctx, has_modif ) ;
      
   }
   assIt->destroy() ; assIt = 0 ;

   PEL_Context const* root_ctx = root->context() ;
   for( size_t i=0 ; i<root_ctx->nb_variables() ; i++ )
   {
      PEL_Variable const* var = root_ctx->variable(i) ;
      PEL_Data const* dat = root_ctx->value(var) ;
      complete_context( root, dat, ref_ctx, has_modif ) ;
   }
   
   PEL_Context const* dup_ctx = dup->context() ;
   for( size_t i=0 ; i<dup_ctx->nb_variables() ; i++ )
   {
      PEL_Variable const* var = dup_ctx->variable(i) ;
      PEL_Data const* dat = dup_ctx->value(var) ;
      complete_context( root, dat, ref_ctx, has_modif ) ;
   }
   
   if( has_modif )
   {
      complete_context( root, dup, ref_ctx ) ;
   }
   else
   {
      PEL_ModuleIterator* modIt = dup->create_module_iterator( 0 ) ;
      for( modIt->start() ; modIt->is_valid() ; modIt->go_next() )
      {
         PEL_Module* a_module = modIt->item() ;
         complete_context( root, a_module, ref_ctx ) ;
      }
      modIt->destroy() ; modIt = 0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_Module:: complete_context( PEL_Module* root,
                               PEL_Data const* dat,
                               PEL_Context const* ref_ctx,
                               bool & has_modif ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: complete_context" ) ;
   PEL_CHECK( root != 0 ) ;
   PEL_CHECK( dat != 0 ) ;
   PEL_CHECK( ref_ctx != 0 ) ;

   PEL_Context const* ctx = root->context() ;
   if( !dat->context_has_required_variables(ctx) )
   { 
      PEL_List* needed = PEL_List::create( 0 ) ;
      dat->declare( needed ) ;
      PEL_ListIterator* it = needed->create_iterator( needed ) ;
      for( it->start() ;
           it->is_valid() ;
           it->go_next() )
      {
         PEL_Variable const* var =
            static_cast<PEL_Variable const*>( it->item() ) ;
         if( ref_ctx->has_variable( var ) && ! ctx->has_variable( var ) )
         {
            has_modif = true ;
            root->add_variable( var, ref_ctx->value( var )->create_clone( 0 ) ) ;
         }   
      }
      needed->destroy() ; needed = 0 ;
   }
}

//----------------------------------------------------------------------
bool
PEL_Module:: find( std::string const& nom,
                   PEL_Module*& theModule,
                   PEL_KeywordDataPair*& theAssignment ) const
//----------------------------------------------------------------------
{
   PEL_CHECK( theModule==0 && theAssignment==0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   bool result = false ;
   
   // The name is null : we have found the name
   if( nom.empty() )
   {
      theModule = const_cast<PEL_Module*>( this ) ;
      result = true ;
   }
   else
   {
      std::string token, otherToken ;
      split(nom, token, otherToken) ;
      
      // The name begins with "/" : we must find the child
      // root = ( nom[0]=='/' ) ;
      PEL_String * str = PEL_String::create( 0, token ) ;
      PEL_Object* obj = MODS->item( str ) ;
      if( obj!=0 )
      {
         // Find a module
         PEL_Module* a_module = static_cast<PEL_Module* >(obj) ;
         if( otherToken.empty() )
         {
            theModule = a_module ;
            result = true ;
         }
         else
         {
            result = a_module->find( otherToken, theModule, theAssignment ) ;
         }
      }
      else
      {
         // Searching for an assignment
         obj = ENTRIES->item( str ) ;
         theAssignment = static_cast<PEL_KeywordDataPair* >(obj) ;
         PEL_CHECK( IMPLIES( obj!=0, dynamic_cast<PEL_KeywordDataPair* >(obj) != 0 ) ) ;
         
         result = theAssignment!=0 && otherToken.empty() ;
      }
      str->destroy() ;

   }
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !result || theModule!=0 || theAssignment!=0 ) ;
   PEL_CHECK_POST( !( theModule!=0 && theAssignment!=0 ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_Module:: split( std::string const& nom,
                    std::string& token,
                    std::string& otherToken ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: split" ) ;
   
   // si nom=="/toto/titi/tutu" ou nom=="toto/titi/tutu" alors 
   //     token="toto" et otherToken = "/titi/tutu"
   // si nom=="/toto" ou nom = "toto" alors
   //     token="toto" et otherToken=""
   int n = (int)nom.length() ;
      
   int start = (int)nom.find_first_not_of( "/" ) ;
   if( ( start >=0 ) && ( start < n ) )
   {
      int stop = (int)nom.find_first_of( "/", start ) ;
      if( ( stop >= 0 ) && ( stop < n ) ) 
      {
         token      = nom.substr( start, stop-start ) ;
         otherToken = nom.substr( stop, n ) ;
      }
      else
      {
         token = nom.substr( start, n-start ) ;
         otherToken = "" ;
      }
   }
   else
   {
      token = nom ;
      otherToken = "" ;
   }
}

//----------------------------------------------------------------------
std::ostream& 
PEL_Module:: recursive_print( std::ostream& s,
                              int n,
                              bool hybrid,
                              std::string const& bin_file) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Module:: recursive_print" ) ;
   PEL_CHECK( IMPLIES( hybrid, !bin_file.empty() ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::string bl( n, ' ' ) ;
   
   s << bl << "MODULE " << NAME->to_string()  << std::endl ;

   // print contexts:
   {
      PEL_Context const* father_ctx = 0 ;
      if( FATHER!=0 )
      {
         father_ctx = FATHER->context() ;
      }
      CTX->print( s, n+2, father_ctx ) ;
   }
   
   // print entries:
   {
      PEL_KeywordDataIterator* assIt = create_entry_iterator( 0 ) ;
      for( assIt->start() ; assIt->is_valid() ; assIt->go_next() )
      {
         PEL_KeywordDataPair const* ass = assIt->item() ;
      
         s << bl << "  " << ass->keyword() << " = " ;
         if( !hybrid ||
             !PEL_BinStored::is_type_supported( ass->data()->data_type() ) ||
             !ass->data()->is_constant() )
         {
            ass->data()->print( s, 0 ) ;
         }
         else
         {
            PEL_Data const* to_print =
               PEL_BinStored::create_reference( 0, ass->data(),
                                                bin_file, true ) ;
            to_print->print( s, 0 ) ;
            to_print->destroy() ;
         }
         s << std::endl ;
      }
      assIt->destroy() ;
   }

   // print modules:
   {
      PEL_ModuleIterator* modIt = create_module_iterator( 0 ) ;
      for( modIt->start() ; modIt->is_valid() ; modIt->go_next() )
      {
         PEL_Module const* a_module = modIt->item() ;
         a_module->recursive_print( s, n+2, hybrid, bin_file ) ;
      }
      modIt->destroy() ;
   }

   
   s << bl << "END MODULE " << name() << std::endl ;
   PEL_CHECK_INV( invariant() ) ;		
   return(s);
}

//----------------------------------------------------------------------
bool
PEL_Module:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( NAME!=0 ) ;
   PEL_ASSERT( MODS!=0 ) ;
   PEL_ASSERT( ENTRIES!=0 ) ;
   // To assume that FATHER hasn't been destroy before self
   PEL_ASSERT( FATHER==0 || FATHER->MODS!=0 ) ;
   return( true ) ;
}
