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

#include <EXT_MPI_API.hh>

#include <PEL_Bool.hh>
#include <PEL_Exec.hh>
#include <PEL_Error.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>
#include <PEL_assertions.hh>

// Both stdio.h and the MPI C++ interface use SEEK_SET, SEEK_CUR, SEEK_END.
// This is really a bug in the MPI-2 standard.
// A possibility would be to undefine the 3 names SEEK_SET, SEEK_CUR, SEEK_END
//    #undef SEEK_SET
//    #undef SEEK_CUR
//    #undef SEEK_END
// Our solution is to define MPICH_IGNORE_CXX_SEEK which works at least 
// with MPICH2
#define MPICH_IGNORE_CXX_SEEK 1
#include <mpi.h>

#include <unistd.h>

#include <iostream>
#include <sstream>

using std::ostringstream ;
using std::endl ;

EXT_MPI_API* EXT_MPI_API:: SINGLETON = new EXT_MPI_API() ;

struct EXT_MPI_API_ERROR
{
   static void n0( std::string const& func ) ;
} ;

//----------------------------------------------------------------------
EXT_MPI_API:: EXT_MPI_API( void )
//----------------------------------------------------------------------
   : PEL_ExternalAPI( "EXT_MPI_API", 8 )
{   
}

//----------------------------------------------------------------------
EXT_MPI_API:: ~EXT_MPI_API( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPI_API:: ~EXT_MPI_API" ) ;
   int mpierr ;

   int size ;
   mpierr = MPI_Comm_size( MPI_COMM_WORLD, &size ) ;
   if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Comm_size" ) ;

   int err = PEL_Exec::exit_code() ;
   if( size>1 && err!=0 )
   {
      int rank ;

      mpierr = MPI_Comm_rank( MPI_COMM_WORLD, &rank ) ;
      if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Comm_rank" ) ;

      std::cout << "Parallel execution interrupted by processor #" << rank 
                << "." << std::endl ;

      mpierr = MPI_Abort( MPI_COMM_WORLD, err ) ;
      if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Abort" ) ;
   }

   mpierr = MPI_Finalize() ;
   if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Finalize" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPI_API:: initialize( int& argc, char **& argv )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPI_API:: initialize" ) ;
   
   int mpierr ;

   std::string const cwd = PEL_System::working_directory() ;

   mpierr = MPI_Init( &argc, &argv ) ;
   if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Init" ) ;

   int size ;
   mpierr = MPI_Comm_size( MPI_COMM_WORLD, &size ) ;
   if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Comm_size" ) ;

   if( size == 1 )
   {
      PEL_System::changedir( cwd ) ;
   }
  
   PEL_Exec::add_variable_to_execution_context(
                      PEL_Variable::object( "BS_with_MPI" ),
                      PEL_Bool::create( 0, true ) ) ;
   
#ifndef MPIRUN
# error \
Macro MPIRUN must be set when compiling Pelicans (opt: -DMPIRUN=<value>).
#else
   if( !PEL_System::can_read( MPIRUN ) )
   {
      ostringstream mesg ;
      mesg << "*** EXT_MPI_API:" << endl ;
      mesg << "    unable to read the file called" << endl ;
      mesg << "       \"" << MPIRUN << "\""  << endl << endl ; 
      mesg << "This file was specified by the macro MPIRUN when" << endl ;
      mesg << "compiling PELICANS, eg via the compiler option:" << endl ;
      mesg << "  -DMPIRUN=" << MPIRUN << endl ;
      mesg << "or in the extra-makefile with an instruction such as:" << endl ;
      mesg << "  MPIRUN=" << MPIRUN << endl << endl ;
      mesg << "The macro MPIRUN should contain the full path of the" << endl ;
      mesg << "command used to run an MPI application" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }
   PEL_Exec::add_variable_to_execution_context(
                      PEL_Variable::object( "SS_MPI_RUN" ),
                      PEL_String::create( 0, MPIRUN ) ) ;
#undef MPIRUN
#endif
}

//internal---------------------------------------------------------------
void
EXT_MPI_API_ERROR:: n0( std::string const& func )
//internal---------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** EXT_MPI_API:" << endl ;
   mesg << "    call to " << func << " failed" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
