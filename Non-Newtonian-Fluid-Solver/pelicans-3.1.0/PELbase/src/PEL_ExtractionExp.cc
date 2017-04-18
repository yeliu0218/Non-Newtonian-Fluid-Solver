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

#include <PEL_ExtractionExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_Communicator.hh>
#include <PEL_Context.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_List.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Sequence.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>
#include <stringVector.hh>

#include <fstream>
#include <iostream>
#include <sstream>

PEL_Module const* PEL_ExtractionExp::DB_MOD = 0 ;
PEL_ExtractionExp* PEL_ExtractionExp::PROTO_HAS_DATA =
   new PEL_ExtractionExp( "has_data", PEL_ExtractionExp::has_data ) ;
PEL_ExtractionExp* PEL_ExtractionExp::PROTO_DATA =
   new PEL_ExtractionExp( "extracted_data", PEL_ExtractionExp::ext_data ) ;
PEL_ExtractionExp* PEL_ExtractionExp::PROTO_HAS_MOD =
   new PEL_ExtractionExp( "has_module", PEL_ExtractionExp::has_mod ) ;
PEL_ExtractionExp* PEL_ExtractionExp::PROTO_EXTRACTED_MODULE =
   new PEL_ExtractionExp( "extracted_module", PEL_ExtractionExp::ext_mod ) ;

//----------------------------------------------------------------------
PEL_ExtractionExp:: PEL_ExtractionExp( std::string const& a_name,
                                       ExtractionExp op ) 
//----------------------------------------------------------------------
   : PEL_TransferExp( a_name )
   , OP( op )
   , TEMP_FILE_NAME( "" )
   , SRC( 0 )
{
   PEL_LABEL( "PEL_ExtractionExp:: PEL_ExtractionExp" ) ;
}

//----------------------------------------------------------------------
PEL_ExtractionExp:: PEL_ExtractionExp( PEL_Object* a_owner,
                                       ExtractionExp op,
                                       std::string const& a_name,
                                       PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_TransferExp( a_owner, a_name, argument_list )
   , OP( op )
   , TEMP_FILE_NAME( "" )
   , SRC( 0 )
{
   PEL_LABEL( "PEL_ExtractionExp:: PEL_ExtractionExp" ) ;
   PEL_ASSERT( is_initialized() ) ;
   
   // Data name:
   std::string const& n = data_name( name(), arg(0) ) ;

   if( OP == ext_data )
   {
      PEL_ASSERT( name() == "extracted_data" ) ;
      if( DB_MOD->has_entry( n ) )
      {
         PEL_Data const* a_data = DB_MOD->data_of_entry( n ) ;
         PEL_List* l = PEL_List::create( 0 ) ;
         a_data->declare( l ) ;
         if( l->count() == 0 )
         {
            SRC = a_data ;
         }
         else
         {
            std::string const& dirname = PEL_Module::dirname( n )  ;
            PEL_Module const* m =
               ( dirname.empty() ? DB_MOD : DB_MOD->module( dirname ) ) ;
            PEL_Data* d =
                 PEL_DataWithContext::create( 0, a_data, m->context() ) ;
            SRC = d->create_simplification( this ) ;
            d->destroy() ; d = 0 ;
         }
         l->destroy() ; l = 0 ;
      }
      else if( argument_list->count() == 2 )
      {
         SRC = arg(1) ;
      }
      else
      {
         raise_error( "    missing entry: "+n ) ;
      }
   }
   else if( OP == has_data ) 
   {
      PEL_ASSERT( name() == "has_data" ) ;
      SRC = PEL_Bool::create( this, DB_MOD->has_entry( n ) ) ;
   }
   else if( OP == has_mod ) 
   {
      PEL_ASSERT( name() == "has_module" ) ;
      SRC = PEL_Bool::create( this, DB_MOD->has_module( n ) ) ;
   }
   else if( OP == ext_mod ) 
   {
      PEL_ASSERT( name() == "extracted_module" ) ;
      if( DB_MOD->has_module( n ) )
      {
         TEMP_FILE_NAME = temporary_file() ;
         extract_module( TEMP_FILE_NAME, n, data_name( name(), arg(1) ) ) ;
         SRC = PEL_String::create( this, TEMP_FILE_NAME ) ;
      }
      else if( argument_list->count() == 3 )
      {
         SRC = arg(2) ;
      }
      else
      {
         raise_error( "    missing module: "+n ) ;
      }
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ExtractionExp:: ~PEL_ExtractionExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: ~PEL_ExtractionExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( OP == ext_mod )
   {
      PEL_System::erase( TEMP_FILE_NAME ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_ExtractionExp:: initialize( PEL_Module const* mod )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: initialize" ) ;
   PEL_CHECK_PRE( mod != 0 ) ;
   
   DB_MOD = mod ;
   
   PEL_CHECK_POST( is_initialized() ) ;
   PEL_CHECK_POST( data_base() == mod ) ;
}

//----------------------------------------------------------------------
void
PEL_ExtractionExp:: reset( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: reset" ) ;
   
   DB_MOD = 0 ;
   
   PEL_CHECK_POST( !is_initialized() ) ;
}

//----------------------------------------------------------------------
bool
PEL_ExtractionExp:: is_initialized( void )
//----------------------------------------------------------------------
{
   return( DB_MOD != 0 ) ;
}

//----------------------------------------------------------------------
PEL_Module const*
PEL_ExtractionExp:: data_base( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: data_base" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;
   PEL_Module const* result = DB_MOD ;
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ExtractionExp*
PEL_ExtractionExp:: create_replica( PEL_Object* a_owner, 
                                    PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   if( !is_initialized() )
   {
      raise_error(
         "*** PEL_ExtractionExp: error\n"
         "    Can't evaluate expression before initialize self" ) ;
   }
   
   PEL_ExtractionExp* result =
           new PEL_ExtractionExp( a_owner, OP, name(), argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_ExtractionExp:: declare( PEL_List* lst ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: declare" ) ;
   PEL_CHECK_PRE( declare_PRE( lst ) ) ;
   
   SRC->declare( lst ) ;
   
   PEL_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_ExtractionExp:: context_has_required_variables(
                                           PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: context_has_required_variables" ) ;
   PEL_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;
   
   return( SRC->context_has_required_variables( ct ) ) ;
}

//----------------------------------------------------------------------
stringVector const&
PEL_ExtractionExp:: undefined_variables( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: undefined_variables" ) ;
   
   return( SRC->undefined_variables( ct ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_ExtractionExp:: value_can_be_evaluated( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: value_can_be_evaluated" ) ;
   
   return( SRC->value_can_be_evaluated( ct ) ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_ExtractionExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: data_type" ) ;
   return( SRC->data_type() ) ;
}

//----------------------------------------------------------------------
bool
PEL_ExtractionExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = true ;
   if( OP == has_data || OP == has_mod )
   {
      result = some_arguments->count()==1 &&
               extract_arg( some_arguments, 0 )->data_type() == String ;
   }
   else if( OP == ext_mod )
   {
      result = ( some_arguments->count()==2 || some_arguments->count()==3 ) &&
               extract_arg( some_arguments, 0 )->data_type() == String &&
               extract_arg( some_arguments, 1 )->data_type() == String ;
      if( result && some_arguments->count()==3  )
      {
         result = ( extract_arg( some_arguments, 2 )->data_type() == String ) ;
      }
   }
   else
   {
      result = ( some_arguments->count()==1 || some_arguments->count()==2 ) &&
               extract_arg( some_arguments, 0 )->data_type() == String ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_ExtractionExp:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: usage" ) ;
   static std::string result ;
   if( OP == has_data )
   {
      result = name() + "(SS)" ;
   }
   else if( OP == has_mod )
   {
      result = name() + "(SS)" ;
   }
   else if( OP == ext_mod )
   {
      result = name() + "(SS,SS,[,SS])" ;
   }
   else
   {
      result = name() + "(SS[,SS])" ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_ExtractionExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: print" ) ;

   if( !is_a_prototype() )
   {
      SRC->print( os, indent_width ) ;
   }
   else
   {
      PEL_Expression::print( os, indent_width ) ;
   }
}

//----------------------------------------------------------------------
PEL_Data const*
PEL_ExtractionExp:: data( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: data" ) ;
   PEL_CHECK( data_PRE( ct ) ) ;

   PEL_Data const* result = SRC ;

   PEL_CHECK_POST( data_POST( result, ct ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_ExtractionExp:: temporary_file( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: temporary_file" ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   
   static int EXTRACTED_IDX = 0 ;
   std::ostringstream ss ;
   ss << PEL_System::working_directory()
      << PEL_System::path_name_separator()
      << "temporary_" << EXTRACTED_IDX++ ;
   if( com->nb_ranks() > 1 )
   {
      ss << "#" << com->rank() ;
   }
   ss << ".pel" ;
   static std::string result = "" ;
   result = ss.str() ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_ExtractionExp:: data_name( std::string const& exp_name,
                               PEL_Data const* d )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: data_name" ) ;
   PEL_CHECK( !exp_name.empty() ) ;
   PEL_CHECK( d != 0 ) ;
   
   if( ! d->value_can_be_evaluated( 0 ) )
   {
      std::ostringstream msg ;
      msg << "*** " << exp_name << " expression error:" << std::endl
          << "    the entry name cannot be defined from variables " << std::endl ;
      stringVector const& undef_vars = d->undefined_variables( 0 ) ;
      if( undef_vars.size() > 0 )
      {
         msg << "    unexpected variable(s): " << std::endl ;
         for( size_t i=0 ; i<undef_vars.size() ; ++i )
         {
            msg << "        - \"" << undef_vars(i) << "\"" << std::endl ;
         }
      }
      PEL_Error::object()->raise_plain( msg.str() ) ;
   }   
   return( d->to_string() ) ;
}

//----------------------------------------------------------------------
void
PEL_ExtractionExp:: extract_module( std::string const& file_name,
                                    std::string const& d_name,
                                    std::string const& m_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ExtractionExp:: extract_module" ) ;
   PEL_CHECK( OP == ext_mod ) ;
   PEL_CHECK( ! file_name.empty() ) ;
   PEL_CHECK( ! d_name.empty() ) ;
   PEL_CHECK( ! m_name.empty() ) ;
   PEL_CHECK( DB_MOD != 0 ) ;
   
   std::ofstream out( file_name.c_str(),
                      std::ios::out | std::ios::trunc ) ;
   if( !out )
   {
     raise_error( "   unable to create temporary file \""+file_name+"\"" ) ;
   }
   PEL_Module* mod = DB_MOD->module( d_name ) ;
   PEL_Module* dup = mod->create_clone(0) ;
   dup->modify_module_name( m_name ) ;
   dup->print( out, 0 ) ;
   out.close() ;
   dup->destroy() ; dup = 0 ;
}
