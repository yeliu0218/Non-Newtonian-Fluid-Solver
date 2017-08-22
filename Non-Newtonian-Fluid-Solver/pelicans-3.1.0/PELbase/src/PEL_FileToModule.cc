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

#include <PEL_FileToModule.hh>

#include <PEL_Error.hh>
#include <PEL_Iterator.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <sstream>

using std::endl ;

struct PEL_FileToModule_ERROR
{
   static void n0( std::string const& a_format,
                   std::string const& a_default_motif,
                   PEL_FileToModule const* other ) ;
} ;

//----------------------------------------------------------------------------
PEL_FileToModule const*
PEL_FileToModule:: object( std::string const& format )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FileToModule:: object" ) ;

   PEL_FileToModule const* result =
      static_cast<PEL_FileToModule const*>(
                                    plugins_map()->item( format ) ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( PEL_Root::object() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_FileToModule:: PEL_FileToModule( std::string const& a_format,
                                     std::string const& a_default_motif )
//----------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , MY_FORMAT( a_format )
   , MY_MOTIF( a_default_motif )
{
   PEL_LABEL( "PEL_FileToModule:: PEL_FileToModule" ) ;
   
   PEL_Iterator* it = plugins_map()->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PEL_FileToModule const* oo = 
                       static_cast< PEL_FileToModule const* >( it->item() ) ;
      if( oo->default_motif() == a_default_motif )
         PEL_FileToModule_ERROR::n0( a_format, a_default_motif, oo ) ;
   }
   it->destroy() ; it = 0 ;

   plugins_map()->register_item( a_format, this ) ;
   if( formats().empty() )
   {
      formats() = a_format ;
   }
   else
   {
      formats() += "," + a_format ;
   }
}

//----------------------------------------------------------------------------
PEL_FileToModule:: ~PEL_FileToModule( void  )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
bool
PEL_FileToModule:: has( std::string const& format )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FileToModule:: has" ) ;
   
   bool result = plugins_map()->has( format ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PEL_FileToModule:: find_file_format( std::string const& a_file_name,
                                     std::string& a_format )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_FileToModule:: find_file_format" ) ;
   
   PEL_Iterator* it = plugins_map()->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PEL_FileToModule const* oo = 
                       static_cast< PEL_FileToModule const* >( it->item() ) ;
      if( a_file_name.find( oo->default_motif() ) < a_file_name.length() )
      {
         a_format = oo->format() ;
         break ;
      }
   }
   it->destroy() ; it = 0 ;
}

//----------------------------------------------------------------------------
std::string const&
PEL_FileToModule:: list_of_formats( void )
//----------------------------------------------------------------------------
{
   return( formats() ) ; 
}

//----------------------------------------------------------------------------
std::string const&
PEL_FileToModule:: format( void ) const
//----------------------------------------------------------------------------
{
   return( MY_FORMAT ) ;
}

//----------------------------------------------------------------------------
std::string const&
PEL_FileToModule:: default_motif( void ) const
//----------------------------------------------------------------------------
{
   return( MY_MOTIF ) ;
}

//----------------------------------------------------------------------------
PEL_ObjectRegister*
PEL_FileToModule:: plugins_map( void )
//----------------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "PEL_FileToModule descendant" ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
std::string&
PEL_FileToModule:: formats( void )
//----------------------------------------------------------------------------
{
   static std::string result ;
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
PEL_FileToModule:: create_from_file_PRE( PEL_Object* a_owner,
                                         std::string const& module_name,
                                         std::string const& file_name ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( !module_name.empty() ) ;
   PEL_ASSERT( !file_name.empty() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------------
bool
PEL_FileToModule:: create_from_file_POST( PEL_Module const* result,
                                          PEL_Object* a_owner,
                                          std::string const& module_name,
                                          std::string const& file_name ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->name() == module_name ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------------
void 
PEL_FileToModule_ERROR:: n0( std::string const& a_format,
                             std::string const& a_default_motif,
                             PEL_FileToModule const* other )
//internal--------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** PEL_FileToModule error:" << endl << endl ;
   mesg << "    Attempt to register the format: \"" << a_format 
        << "\"" << endl ;
   mesg << "    associated to the filename motif: \"" 
        << a_default_motif << "\"" << endl << endl ;
   mesg << "    The filename motif \""
        << a_default_motif << "\" is already handled" << endl ;
   mesg << "    by the format: \"" << other->format() << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

