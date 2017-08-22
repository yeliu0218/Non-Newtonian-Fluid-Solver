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

#include <PEL_MemoryTracer.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Root.hh>
#include <PEL_System.hh>

#include <iostream>
#include <fstream>
#include <sstream>

struct PEL_MemoryTracer_ERROR
{
   static void n0( std::string const& file_name ) ;
   static void n1( void ) ;
} ;

//----------------------------------------------------------------------
PEL_MemoryTracer*
PEL_MemoryTracer:: object( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: object" ) ;

   static PEL_MemoryTracer* result = new PEL_MemoryTracer() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_MemoryTracer:: PEL_MemoryTracer( void )
//----------------------------------------------------------------------
   : PEL_Object( PEL_Root::object() )
   , MEMORY_TRACE( false )
   , OFILENAME( "" )
   , INDENT( "" )
   , EVENTS( 0 )
   , MEM0( 0 )
   , OBJ0( 0 )
{
}

//----------------------------------------------------------------------
PEL_MemoryTracer:: ~PEL_MemoryTracer( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_MemoryTracer:: enable_memory_trace( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: enable_memory_trace" ) ;
   PEL_CHECK_PRE( !memory_trace_enabled() ) ;

   init_trace_file() ;
   PEL::out() << "*** PEL_MemoryTracer: memory trace enabled" << std::endl
              << "    trace_file: " << OFILENAME  << std::endl ;
   
   MEMORY_TRACE = true ;
   
   PEL_CHECK_POST( memory_trace_enabled() ) ;
}

//----------------------------------------------------------------------
void
PEL_MemoryTracer:: disable_memory_trace( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: disable_memory_trace" ) ;
   PEL_CHECK_PRE( memory_trace_enabled() ) ;

   PEL::out() << "*** PEL_MemoryTracer: memory trace disabled" << std::endl ;
   
   MEMORY_TRACE = false ;
   
   PEL_CHECK_POST( !memory_trace_enabled() ) ;
}

//----------------------------------------------------------------------
bool
PEL_MemoryTracer:: memory_trace_enabled( void ) const
//----------------------------------------------------------------------
{
   return( MEMORY_TRACE ) ;
}

//----------------------------------------------------------------------
size_t
PEL_MemoryTracer:: used_memory( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: used_memory" ) ;
   return( PEL_System::used_memory() ) ;
}

//----------------------------------------------------------------------
void
PEL_MemoryTracer:: display_memory( std::ostream& os, size_t memory )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: display_memory" ) ;
   PEL_CHECK_PRE( os ) ;

   static size_t const mo = 1024*1024 ;
   static size_t const go = 1024*1024*1024 ;

   if( memory > go )
   {
      os << ( (double) memory )/go << " Go" ;
   }
   else if( memory > mo )
   {
      os << ( (double) memory )/mo << " Mo" ;
   }
   else
   {
      os << memory << " octets" ;
   }
}

//----------------------------------------------------------------------
void
PEL_MemoryTracer:: start_event( std::string const& label )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: start_event" ) ;
   PEL_CHECK_PRE( !label.empty() ) ;

   if( MEMORY_TRACE )
   {
      std::ofstream file( OFILENAME.c_str(), std::ios::out | std::ios::app ) ;
      if( !file )
      {
         PEL_MemoryTracer_ERROR:: n0( OFILENAME ) ;
      }
      size_t const mem = used_memory() ;
      size_t const nb_objs = PEL_Object::GetNumberOf_PEL_objects() ;
      file << INDENT << "### Start: " << label << " (memory: " ;
      display_memory( file, mem ) ;
      file << ", objects: " << nb_objs << ")" << std::endl ;
      file << INDENT << "#" << std::endl ;
      EVENTS.append( label ) ;
      MEM0.append( mem ) ;
      OBJ0.append( nb_objs ) ;
      INDENT += "   " ;
   }
}

//----------------------------------------------------------------------
void
PEL_MemoryTracer:: stop_event( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: stop_event" ) ;
   
   if( MEMORY_TRACE )
   {
      if( INDENT.size() < 3 ) PEL_MemoryTracer_ERROR:: n1() ;
      INDENT.erase( INDENT.length()-3 ) ;
      std::ofstream file( OFILENAME.c_str(), std::ios::out | std::ios::app ) ;
      if( !file )
      {
         PEL_MemoryTracer_ERROR:: n0( OFILENAME ) ;
      }
      size_t const mem = used_memory() ;
      size_t const nb_objs = PEL_Object::GetNumberOf_PEL_objects() ;
      size_t const last = EVENTS.size()-1 ;
      file << INDENT << "#          diff memory: " ;
      if( mem>=MEM0(last) )
      {
         display_memory( file, mem-MEM0(last) ) ;
      }
      else
      {
         file << "-" ;
         display_memory( file, MEM0(last)-mem ) ;
      }
      file << ", diff objects: " ;
      if( nb_objs>=OBJ0(last) )
      {
         file << nb_objs-OBJ0(last) ;
      }
      else
      {
         file << "-" << OBJ0(last)-nb_objs ;
      }
      file << std::endl ;
      file << INDENT << "### Stop:  " ;
      file << EVENTS(last) << " (memory: " ;
      display_memory( file, mem ) ;
      file << ", objects: " << nb_objs << ")" << std::endl ;
      file << std::endl ;
      EVENTS.resize( last ) ;
      MEM0.resize( last ) ;
      OBJ0.resize( last ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_MemoryTracer:: trace( std::string const& a_message )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: trace" ) ;
   PEL_CHECK_PRE( !a_message.empty() ) ;

   if( MEMORY_TRACE )
   {
      std::ofstream file( OFILENAME.c_str(), std::ios::out | std::ios::app ) ;
      if( !file )
      {
         PEL_MemoryTracer_ERROR:: n0( OFILENAME ) ;
      }
      file << INDENT << a_message << std::endl ;
      file << std::endl ;
      file.close() ;
   }
}

//----------------------------------------------------------------------
std::string const&
PEL_MemoryTracer:: message( std::string const& label ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: message" ) ;
   static std::string result ;
   std::ostringstream msg ;
   msg << INDENT << "# " << label << " (memory_used: " ;
   display_memory( msg, used_memory() ) ;
   msg << ", "
       << PEL_Object::GetNumberOf_PEL_objects() << " objects)" ;
   result = msg.str() ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_MemoryTracer:: indent( void ) const
//----------------------------------------------------------------------
{
   return( INDENT ) ;
}

//----------------------------------------------------------------------
void
PEL_MemoryTracer:: init_trace_file( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MemoryTracer:: init_trace_file" ) ;

   if( OFILENAME.empty() )
   {
      std::stringstream m ;
      m << "memory" ;
      PEL_Communicator const* com = PEL_Exec::communicator() ;
      if( com->nb_ranks() > 1 )
      {
         m << "#" << com->rank() ;
      }
      m << ".txt" ;
      OFILENAME = m.str() ;
      std::ofstream file( OFILENAME.c_str(),
                          std::ios::out | std::ios::trunc ) ;
      if( !file )
      {
         PEL_MemoryTracer_ERROR:: n0( OFILENAME ) ;
      }
      else
      {
         std::string const s = std::string( 60, '#' ) ;
         file << s << std::endl ;
         file << "#" << std::endl ;
         file << "# PEL_MemoryTracer generated file" << std::endl ;
         file << "#" << std::endl ;
         file << s << std::endl ;
         file << std::endl ;
      }
      file.close() ;
   }
}

//internal--------------------------------------------------------------
void
PEL_MemoryTracer_ERROR:: n0( std::string const& file_name )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** PEL_MemoryTracer error:\n"
      "    Unable to open file \""+file_name+"\" for writing" ) ;
}

//internal--------------------------------------------------------------
void
PEL_MemoryTracer_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** PEL_MemoryTracer error:\n"
      "    attempt to decrease indentation below the zero limit" ) ;
}
