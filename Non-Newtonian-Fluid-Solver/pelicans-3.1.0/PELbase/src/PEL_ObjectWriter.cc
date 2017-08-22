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

#include <PEL_ObjectWriter.hh>

#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_Exec.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_Int.hh>
#include <PEL_Error.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_ListIterator.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_String.hh>

#include <fstream>
#include <sstream>

struct PEL_ObjectWriter_ERROR
{
   static void n0( void ) ;
   static void n1( PEL_ModuleExplorer const* exp, std::string const& n ) ;
   static void n2( std::string const& n ) ;
} ;

//---------------------------------------------------------------------------
PEL_ObjectWriter*
PEL_ObjectWriter:: create( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp,
                           PEL_ModuleExplorer const* header_exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   std::string const& w_type = exp->string_data( "type" ) ;
   PEL_ObjectWriter::PEL_ObjectWriterType writer_type =
                                        PEL_ObjectWriter::last_two_cycles ;
   
   if( w_type == "all_cycles_in_one_file" )
   {
      writer_type = PEL_ObjectWriter::all_cycles ;
   }
   else if( w_type == "cycles_in_separate_files" )
   {
      writer_type = PEL_ObjectWriter::per_one_cycle ;
   }
   else if( w_type == "last_two_cycles" )
   {
      writer_type = PEL_ObjectWriter::last_two_cycles ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value(
         exp,
         "type",
         "   - \"all_cycles_in_one_file\"\n"
         "   - \"cycles_in_separate_files\"\n"
         "   - \"last_two_cycles\"" ) ;
   }
   
   PEL_ObjectWriter* result = new PEL_ObjectWriter( a_owner, writer_type,
                                                    exp, header_exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->has_an_opened_cycle() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_ObjectWriter:: PEL_ObjectWriter( PEL_Object* a_owner,
                                     PEL_ObjectWriterType const writer_type,
                                     PEL_ModuleExplorer const* exp,
                                     PEL_ModuleExplorer const* header_exp )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , TYPE( writer_type )
   , OFILE_FORMAT( exp->has_entry( "output_format" ) ?
                         exp->string_data( "output_format" ) : "hybrid" )
   , HEADER_EXP( header_exp != 0 ? header_exp->create_clone( this ) : 0 )
   , OFILE_NAME( )
   , OFILE_NAME1( )
   , OFILE_NAME2( )
   , iCYCLE( 0 )
   , NB_OBJECTS( 0 )
{
   if( exp->has_entry( "output_format" ) )
   {
      if( OFILE_FORMAT!="text" && OFILE_FORMAT!="hybrid" )
      {
	 PEL_Error::object()->raise_bad_data_value(
	    exp,
	    "output_format",
	    "   - \"text\"\n   - \"hybrid\"" ) ;
      }
   }

   if( TYPE == PEL_ObjectWriter::all_cycles )
   {
      OFILE_NAME = exp->string_data( "file_name" ) ;
      if( OFILE_NAME.empty() ) PEL_ObjectWriter_ERROR:: n1( exp, "file_name" ) ;
      PEL_Communicator const* com = PEL_Exec::communicator() ;
      if( com->nb_ranks()>1 )
      {
         std::ostringstream rank ;
         rank << "." << com->rank() ;
         OFILE_NAME += rank.str() ;
      }
   }
   else if( TYPE == PEL_ObjectWriter::per_one_cycle )
   {
      OFILE_NAME1 = exp->string_data( "files_basename" ) ;
      if( OFILE_NAME1.empty() ) PEL_ObjectWriter_ERROR:: n1( exp, "files_basename" ) ;
   }
   else if( TYPE == PEL_ObjectWriter::last_two_cycles )
   {
      OFILE_NAME1 = exp->string_data( "file_name_0" ) ;
      OFILE_NAME2 = exp->string_data( "file_name_1" ) ;
      if( OFILE_NAME1.empty() ) PEL_ObjectWriter_ERROR:: n1( exp, "file_name_0" ) ;
      if( OFILE_NAME2.empty() ) PEL_ObjectWriter_ERROR:: n1( exp, "file_name_1" ) ;
      PEL_Communicator const* com = PEL_Exec::communicator() ;
      if( com->nb_ranks()>1 )
      {
         std::ostringstream rank ;
         rank << "." << com->rank() ;
         OFILE_NAME1 += rank.str() ;
         OFILE_NAME2 += rank.str() ;
      }
      OFILE_NAME = OFILE_NAME1 ;
   }
}

//---------------------------------------------------------------------------
PEL_ObjectWriter:: ~PEL_ObjectWriter( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PEL_ObjectWriter:: start_cycle( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: start_cycle" ) ;
   PEL_CHECK_PRE( !has_an_opened_cycle() ) ;
   PEL_SAVEOLD( size_t, cycle_number, cycle_number() ) ;
   
   ++iCYCLE ;
   
   initialize_saving_file() ;

   std::ostringstream mn ;
   mn << "cycle#" << iCYCLE ;
   PEL_Module* mod = PEL_Module::create( this, mn.str() ) ;

   mod->add_entry( "cycle_number", PEL_Int::create( mod, iCYCLE ) ) ;

   MODS.push( mod ) ;

   PEL_CHECK_POST( has_an_opened_cycle() ) ;
   PEL_CHECK_POST( cycle_number() == OLD( cycle_number ) + 1 ) ;
   PEL_CHECK_POST( current_object_number() == 0 ) ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectWriter:: terminate_cycle( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: terminate_cycle" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;

   PEL_Module* mod = MODS.top() ;
   MODS.pop() ;

   if( !MODS.empty() ) PEL_ObjectWriter_ERROR::n0() ;

   mod->write( OFILE_NAME, OFILE_FORMAT ) ;

   destroy_possession( mod ) ;

   NB_OBJECTS=0 ;

   PEL_CHECK_POST( !has_an_opened_cycle() ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ObjectWriter:: has_an_opened_cycle( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: has_an_opened_cycle" ) ;
   return( !MODS.empty() ) ;
}

//---------------------------------------------------------------------------
size_t
PEL_ObjectWriter:: cycle_number( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: cycle_number" ) ;
   return( iCYCLE ) ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectWriter:: start_new_object( std::string const& class_name )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: start_new_object" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_SAVEOLD( size_t, current_object_number, current_object_number() ) ;

   ++NB_OBJECTS ;

   std::ostringstream name ;
   name << "object#" << NB_OBJECTS ;

   PEL_Module* mod = PEL_Module::create( MODS.top(), name.str() ) ;
   mod->add_entry( "class", PEL_String::create( mod, class_name ) ) ;
   mod->add_entry( "object_number", PEL_Int::create( mod, (int)NB_OBJECTS ) ) ;
   MODS.top()->add_module( mod ) ;
   MODS.push( mod ) ;

   PEL_CHECK_POST( current_object_number()==OLD(current_object_number) + 1 ) ;
}

//---------------------------------------------------------------------------
size_t
PEL_ObjectWriter:: current_object_number( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: current_object_number" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;

   return( MODS.size()-1 ) ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectWriter:: add_entry( std::string const& keyword, PEL_Data* data )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: add_entry" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_CHECK_PRE( current_object_number() != 0 ) ;
   PEL_CHECK_PRE( data->owner() == 0 ) ;

   data->set_owner( MODS.top() ) ;
   MODS.top()->add_entry( keyword, data ) ;

   PEL_CHECK_POST( data->is_under_ownership_of( this ) ) ;
}

//---------------------------------------------------------------------------
void PEL_ObjectWriter:: finalize_object( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "void PEL_ObjectWriter:: finalize_object" ) ;
   PEL_CHECK_PRE( has_an_opened_cycle() ) ;
   PEL_SAVEOLD( size_t, current_object_number, current_object_number() ) ;

   MODS.pop() ;

   PEL_CHECK_POST( current_object_number() == OLD(current_object_number)-1 ) ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectWriter:: initialize_saving_file( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: initialize_saving_file" ) ;

   static bool first = true ;
   
   if( first || TYPE != all_cycles )
   {
      set_file_name() ;
      
      std::ofstream file( OFILE_NAME.c_str(), std::ios::out | std::ios::trunc ) ;
      if( !file ) PEL_ObjectWriter_ERROR:: n2( OFILE_NAME ) ;
      file.close() ;
      if( OFILE_FORMAT=="hybrid" )
      {
         std::string const bin_file_name =  OFILE_NAME+ ".bin" ;
         std::ofstream file_bin( bin_file_name.c_str(),
                                 std::ios::out |
                                 std::ios::binary | 
                                 std::ios::trunc ) ;
         if( !file_bin ) PEL_ObjectWriter_ERROR:: n2( bin_file_name ) ;
         file_bin.close() ;
      }

      std::ofstream os( OFILE_NAME.c_str(), std::ios::out | std::ios::app ) ;
      os.close() ;
      write_communicator() ;
      write_header() ;
   }
   first = false ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectWriter:: set_file_name( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectWriter:: set_file_name" ) ;

   if( TYPE == all_cycles )
   {
      // Nothing to do
   }
   else if( TYPE == per_one_cycle )
   {
      PEL_ASSERT( iCYCLE<99999 ) ;
      std::ostringstream tmp ;
      tmp << iCYCLE ;
      std::string nb_string = tmp.str() ;
      OFILE_NAME = OFILE_NAME1+".00000";
      OFILE_NAME.replace( OFILE_NAME.length()-nb_string.length(),
                          nb_string.length(), nb_string ) ;
      OFILE_NAME += ".pel" ;
      PEL_Communicator const* com = PEL_Exec::communicator() ;
      if( com->nb_ranks()>1 )
      {
         std::ostringstream rank ;
         rank << "." << com->rank() ;
         OFILE_NAME += rank.str() ;
      }
   }
   else if( TYPE == last_two_cycles )
   {
      if( OFILE_NAME == OFILE_NAME1 )
      {
         OFILE_NAME = OFILE_NAME2 ;
      }
      else if( OFILE_NAME == OFILE_NAME2 )
      {
         OFILE_NAME = OFILE_NAME1 ;
      }
   }
}

//---------------------------------------------------------------------------
void
PEL_ObjectWriter:: write_header( void ) const
//---------------------------------------------------------------------------
{
   std::ofstream os( OFILE_NAME.c_str(), std::ios::out | std::ios::app ) ;
   os << "MODULE header" << std::endl ;
   os.close() ;

   if( HEADER_EXP != 0 )  HEADER_EXP->write( OFILE_NAME, OFILE_FORMAT ) ;

   os.open( OFILE_NAME.c_str(), std::ios::out | std::ios::app ) ;
   os << "END MODULE header" << std::endl ;
   os.close() ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectWriter:: write_communicator( void ) const
//---------------------------------------------------------------------------
{
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   std::ofstream os( OFILE_NAME.c_str(), std::ios::out | std::ios::app ) ;
   os << "MODULE communicator" << std::endl ;
   os << "   nb_ranks = " << com->nb_ranks() << std::endl ;
   os << "   rank = " << com->rank() << std::endl ;
   os << "END MODULE communicator" << std::endl ;
   os.close() ;
}

//---------------------------------------------------------------------------
bool
PEL_ObjectWriter:: invariant( void ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( MODS.size() <= NB_OBJECTS ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectWriter_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "PEL_ObjectWriter :" << std::endl
        << "   impossible to terminate a cycle when an " << std::endl 
        << "   object storage is in progress" << std::endl
        << "   (\"start_new_object\" and \"finalize_object\"" 
        << std::endl
        << "    should be called the same number of times between two calls to"
	<< std::endl
	<< "    \"start_cycle\" and \"terminate_cycle\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectWriter_ERROR:: n1( PEL_ModuleExplorer const* exp,
                             std::string const& n )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_bad_data_value(
      exp, n, "a no empty string is expected" ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectWriter_ERROR:: n2( std::string const& n )
//internal--------------------------------------------------------------
{
   std::string mess ;
   mess += "PEL_ObjectWriter :\n" ;
   mess += "   Saving failure : unable to open file \"" ;
   mess += n ;
   mess += "\" for writing" ;
   PEL_Error::object()->raise_plain( mess ) ;
}

