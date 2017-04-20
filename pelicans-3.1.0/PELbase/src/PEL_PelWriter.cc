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

#include <PEL_PelWriter.hh>

#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <sstream>
#include <fstream>

bool PEL_PelWriter::APPEND = false ;

PEL_PelWriter const* PEL_PelWriter::PROTOTYPE = new PEL_PelWriter() ;

//----------------------------------------------------------------------
PEL_PelWriter:: PEL_PelWriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "PEL_PelWriter" )
   , FILENAME()
   , FORMAT()
   , WRITER_NAME()
{
}

//----------------------------------------------------------------------
PEL_PelWriter*
PEL_PelWriter:: create_replica( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelWriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_PelWriter* result = new PEL_PelWriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PEL_PelWriter:: PEL_PelWriter( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , FILENAME( exp->string_data( "files_basename" )+".pel" )
   , FORMAT( exp->string_data( "writing_mode" ) )
   , WRITER_NAME()
{
   // Meshing name :
   static size_t writer_count = 0 ;
   size_t writer_id = writer_count++ ;
   std::ostringstream mm ;
   mm << "PEL_PelWriter#" << writer_id ;
   WRITER_NAME = mm.str() ;
   
   // Format :
   if( !( FORMAT=="text" || FORMAT=="binary" ) )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "writing_mode", "\"text\" or \"binary\"" ) ;
   }
   if( FORMAT=="binary" )
   {
      FORMAT = "hybrid" ;
   }
   
   bool append_mode = false ;
   if( exp->has_entry( "append_mode" ) )
   {
      append_mode = exp->bool_data( "append_mode" ) ;
      APPEND = true ;
   }
   if( append_mode )
   {
      std::ofstream file( FILENAME.c_str(),
                          std::ios::out | std::ios::app ) ;
      if( !file )
      {
         std::string mess = "PEL_PelWriter : unable to open file \"" ;
         mess += FILENAME ;
         mess += "\"" ;
         PEL_Error::object()->raise_plain( mess ) ;
      }
      file.close() ;
   }
   else
   {
      std::ofstream file( FILENAME.c_str(),
                          std::ios::out | std::ios::trunc ) ;
      if( !file )
      {
         std::string mess = "PEL_PelWriter : unable to open file \"" ;
         mess += FILENAME ;
         mess += "\"" ;
         PEL_Error::object()->raise_plain( mess ) ;
      }
      file.close() ;
   }

   PEL_CHECK_INV( invariant() ) ;
}


//----------------------------------------------------------------------
PEL_PelWriter:: ~PEL_PelWriter( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_PelWriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelWriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;

   PEL_ModuleExplorer const* dummy_exp = exp ;
   PEL_Module* dummy_m = 0 ;
   if( APPEND )
   {
      dummy_m = PEL_Module::create( 0, WRITER_NAME ) ;
      PEL_Module* m = PEL_Module::create( dummy_m, exp->name() ) ;
      PEL_ModuleExplorer* e = exp->create_clone( dummy_m ) ; 
      add_modules( e, m ) ;
      add_entries( e, m ) ;
      dummy_m->add_module( m ) ;
      dummy_exp = PEL_ModuleExplorer::create( dummy_m, dummy_m ) ;
   }
   
   dummy_exp->write( FILENAME, FORMAT ) ;

   if( dummy_m!=0 )
   {
      dummy_m->destroy() ; dummy_m = 0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_PelWriter:: add_modules( PEL_ModuleExplorer* exp,
                             PEL_Module* m ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelWriter:: add_modules" ) ;

   exp->start_module_iterator() ;   
   for( ; exp->is_valid_module() ; exp->go_next_module() )
   {
      PEL_ModuleExplorer* e = exp->create_subexplorer( 0 ) ;
      PEL_Module* mm = PEL_Module::create( m, e->name() ) ;
      add_entries( e, mm ) ;
      add_modules( e, mm ) ;
      m->add_module( mm ) ;
      e->destroy() ; e = 0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_PelWriter:: add_entries( PEL_ModuleExplorer* exp,
                             PEL_Module* m ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_PelWriter:: add_entries" ) ;

   exp->start_entry_iterator() ;   
   for( ; exp->is_valid_entry() ; exp->go_next_entry() )
   {
      m->add_entry( exp->keyword(), exp->data( m ) ) ;
   }
}
