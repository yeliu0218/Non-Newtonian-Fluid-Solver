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

#include <EXT_MPIcommunicator.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_assertions.hh>
#include <doubleVector.hh>
#include <intVector.hh>

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

#include <iostream>
#include <sstream>

using std::ostringstream ;
using std::endl ;

EXT_MPIcommunicator*
EXT_MPIcommunicator:: SINGLETON = new EXT_MPIcommunicator() ;
   
struct EXT_MPIcommunicator_ERROR
{
   static void n0( std::string const& func ) ;
} ;

#include <PEL_Timer.hh>
#include <PEL_Map.hh>
#include <PEL_Root.hh>
#include <PEL_Root.hh>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


//----------------------------------------------------------------------
EXT_MPIcommunicator:: EXT_MPIcommunicator( void )
//----------------------------------------------------------------------
   : PEL_Communicator( "EXT_MPIcommunicator" )
{
}

//----------------------------------------------------------------------
EXT_MPIcommunicator:: ~EXT_MPIcommunicator( void )
//----------------------------------------------------------------------
{
   SINGLETON = 0 ;
   
}

//----------------------------------------------------------------------
size_t
EXT_MPIcommunicator:: rank( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: rank" ) ;
   
   static size_t result = PEL::bad_index() ;
   if( result==PEL::bad_index() )
   {
      int tmp ;
      int  mpierr = MPI_Comm_rank( MPI_COMM_WORLD, &tmp ) ;
      if( mpierr != MPI_SUCCESS ) 
         EXT_MPIcommunicator_ERROR::n0( "MPI_Comm_rank" ) ;
      result = (size_t) tmp ;
   }

   PEL_CHECK_POST( rank_POST( result ) ) ;   
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
EXT_MPIcommunicator:: nb_ranks( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: nb_ranks" ) ;

   static size_t result = PEL::bad_index() ;
   if( result==PEL::bad_index() )
   {
      int tmp ;
      int  mpierr = MPI_Comm_size( MPI_COMM_WORLD, &tmp ) ;
      if( mpierr != MPI_SUCCESS ) 
         EXT_MPIcommunicator_ERROR::n0( "MPI_Comm_size" ) ;
      result = (size_t) tmp ;
   }

   PEL_CHECK_POST( nb_ranks_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: send( size_t dest, int const* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: send( int const* )" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;

   int mpierr =  MPI_Send( (void*)const_cast<int*>(value),
                           nb, MPI_INT, dest, TAG_INT, MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS )
      EXT_MPIcommunicator_ERROR::n0( "MPI_Send(int*)" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: send( size_t dest, double const* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: send( double const* )" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;
   
   int mpierr = MPI_Send( (void*)const_cast<double*>(value),
                          nb, MPI_DOUBLE, dest, TAG_DOUBLE, MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Send(double*)" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: send( size_t dest, char const* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: send( char const* )" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;
   
   int mpierr = MPI_Send( (void*)const_cast<char*>(value),
                          nb, MPI_CHAR, dest, TAG_CHAR, MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Send(char*)" ) ;
}

//----------------------------------------------------------------------
void*
EXT_MPIcommunicator:: Isend( size_t dest, int const* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: Isend( int const* )" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;
   
   MPI_Request* request = new MPI_Request ;
   
   int mpierr =  MPI_Isend( (void*)const_cast<int*>(value),
                            nb, MPI_INT, dest, TAG_INT, MPI_COMM_WORLD,
                            request ) ;
   if( mpierr != MPI_SUCCESS )
      EXT_MPIcommunicator_ERROR::n0( "MPI_Isend(int*)" ) ;
   
   return( request ) ;
   
}

//----------------------------------------------------------------------
void*
EXT_MPIcommunicator:: Ireceive( size_t src, int* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: Ireceive( int* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;
   MPI_Request* request = new MPI_Request ;
   
   int mpierr =  MPI_Irecv( (void*)value,
                            nb, MPI_INT, src, TAG_INT, MPI_COMM_WORLD,
                            request ) ;
   if( mpierr != MPI_SUCCESS )
      EXT_MPIcommunicator_ERROR::n0( "MPI_Irecv(int*)" ) ;
   
   return( request ) ;
   
}

//----------------------------------------------------------------------
void*
EXT_MPIcommunicator:: Isend( size_t dest, double const* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: Isend( double const* )" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;
   MPI_Request* request = new MPI_Request ;
   
   int mpierr = MPI_Isend( (void*)const_cast<double*>(value),
                           nb, MPI_DOUBLE, dest, TAG_DOUBLE, MPI_COMM_WORLD,
                           request ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Isend(double*)" ) ;
   
   return( request ) ;
}

//----------------------------------------------------------------------
void*
EXT_MPIcommunicator:: Ireceive( size_t src, double* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: Ireceive( double* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;
   MPI_Request* request = new MPI_Request ;
   
   int mpierr = MPI_Irecv( (void*)value,
                           nb, MPI_DOUBLE, src, TAG_DOUBLE, MPI_COMM_WORLD,
                           request ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Irecv(double*)" ) ;

   return( request ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: wait( void* request  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: wait" ) ;
   
   MPI_Status status ;
   MPI_Request* r = (MPI_Request*) request ;
   
   int mpierr = MPI_Wait( r, &status ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_wait" ) ;

   delete (MPI_Request*)r ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: receive( size_t src, int* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: receive( int const* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;

   static MPI_Status status ;
   int mpierr = MPI_Recv( (void*)value, nb, MPI_INT, src, TAG_INT,
                          MPI_COMM_WORLD, &status ) ;
   if( mpierr!=MPI_SUCCESS || status.MPI_ERROR!=MPI_SUCCESS )
      EXT_MPIcommunicator_ERROR::n0( "MPI_Recv(int*)" ) ;

   PEL_CHECK( MPI_Get_count( &status, MPI_INT, &count )==MPI_SUCCESS &&
              count==nb ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: receive( size_t src, double* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: receive( double* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;

   static MPI_Status status ;
   int mpierr = MPI_Recv( (void*)value, nb, MPI_DOUBLE, src, TAG_DOUBLE,
                          MPI_COMM_WORLD, &status ) ;
   if( mpierr!=MPI_SUCCESS || status.MPI_ERROR!=MPI_SUCCESS )
      EXT_MPIcommunicator_ERROR::n0( "MPI_Recv(double*)" ) ;

   PEL_CHECK( MPI_Get_count( &status, MPI_DOUBLE, &count )==MPI_SUCCESS &&
              count==nb ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: receive( size_t src, char* value, int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: receive( char* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;

   static MPI_Status status ;
   int mpierr = MPI_Recv( (void*)value, nb, MPI_CHAR, src, TAG_CHAR,
                          MPI_COMM_WORLD, &status ) ;
   if( mpierr!=MPI_SUCCESS || status.MPI_ERROR!=MPI_SUCCESS )
      EXT_MPIcommunicator_ERROR::n0( "MPI_Recv" ) ;

   PEL_CHECK( MPI_Get_count( &status, MPI_CHAR, &count )==MPI_SUCCESS &&
              count==nb ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: all_gather( int const* value,
                                  size_t nb,
                                  int* result  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: all_gather(int)" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   int mpierr =  MPI_Allgather( const_cast<void*>((void const*)value), 
                                nb, MPI_INT,
                                result, nb, MPI_INT, MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Allgather" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator::  all_gather_v( double const* values,
                                     size_t nb,
                                     double* result,
                                     intVector const& partition,
                                     intVector const& start ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: all_gather_v(double vector)" ) ;
   PEL_CHECK( partition.size()==nb_ranks() ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   int mpierr =  MPI_Allgatherv( const_cast<void*>((void const*)values), 
                                 nb, MPI_DOUBLE,
                                 result, const_cast<int*>(partition.data()),
                                 const_cast<int*>(start.data()),
                                 MPI_DOUBLE, MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Allgatherv" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: all_to_all( int const* value, size_t nb,
                                  int* result  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: all_to_all(int)" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   int mpierr =  MPI_Alltoall( const_cast<void*>((void const*)value), 
                               nb, MPI_INT,
                               result, nb, MPI_INT, MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Alltoall" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: gather( double const* value,size_t nb,
                              double* result, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: gather(double)" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   int mpierr = MPI_Gather( const_cast<void*>((void const*)value), 
                            nb, MPI_DOUBLE,
                            result, nb, MPI_DOUBLE, root, MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Gather" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: gather( int const* value,size_t nb,
                              int* result, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: gather(int)" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   int mpierr = MPI_Gather( const_cast<void*>((void const*)value), 
                            nb, MPI_INT,
                            result, nb, MPI_INT, root, MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Gather" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: broadcast( int* value, size_t nb, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: broadcast(int*)" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   int mpierr =  MPI_Bcast( (void*)(value), nb, MPI_INT, root,
                            MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_broadcast(int*)" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: broadcast( double* value, size_t nb, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: broadcast(double*)" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   int mpierr =  MPI_Bcast( (void*)(value), nb, MPI_DOUBLE, root,
                            MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_broadcast(double*)" ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: broadcast( char* value, size_t nb, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: broadcast(char*)" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   int mpierr =  MPI_Bcast( (void*)(value), nb, MPI_CHAR, root,
                            MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_broadcast(char*)" ) ;
}

//----------------------------------------------------------------------
bool
EXT_MPIcommunicator:: same_value_everywhere( int val ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: same_value_everywhere(int)" ) ;
   // PEL_CHECK_COLLECTIVE(true) ; //??? Recursive call of
                                   //??? PEL_Maker::is_collective()

   intVector other( nb_ranks() ) ;
   int aval=val ;
   
   int mpierr = MPI_Gather(
                    const_cast<void*>( (void const*) &aval ),
                    1, MPI_INT,
                    (void*) other.data(), 1, MPI_INT, 0,
                    MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Gather" ) ;
   
   int ok = 1 ;
   if( rank()==0 )
   {
      
      for( size_t i=0 ; i<nb_ranks() ; i++ )
      {
         if( other(i)!=val )
         {
            ok = 0 ;
         }
      }
   }   
   mpierr =  MPI_Bcast( (void*)(&ok), 1, MPI_INT, 0,
                        MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Bcast" ) ;
   
   return( ok==1 ) ;
}

//----------------------------------------------------------------------
bool
EXT_MPIcommunicator:: same_value_everywhere( double val ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: same_value_everywhere(double)" ) ;
   PEL_CHECK_COLLECTIVE(true) ;

   doubleVector other( nb_ranks() ) ;
   double aval=val ;
   
   int mpierr = MPI_Gather(
                    const_cast<void*>( (void const*) &aval ),
                    1, MPI_DOUBLE,
                    (void*) other.data(), 1, MPI_DOUBLE, 0,
                    MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Gather" ) ;
   
   int ok = 1 ;
   if( rank()==0 )
   {
      
      for( size_t i=0 ; i<nb_ranks() ; i++ )
      {
         if( other(i)!=val )
         {
            ok = 0 ;
         }
      }
   }   
   mpierr =  MPI_Bcast( (void*)(&ok), 1, MPI_INT, 0,
                        MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Bcast" ) ;
   
   return( ok==1 ) ;
}

//----------------------------------------------------------------------
bool
EXT_MPIcommunicator:: boolean_and( bool value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: boolean_and" ) ;  
   PEL_CHECK_COLLECTIVE(true) ;
   
   int input = ( int ) value ;
   int result ;

   int mpierr = MPI_Allreduce( (void*)(&input), (void*)(&result),
                               1, MPI_INT, MPI_LAND, MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Allreduce and" ) ;
   
   return( (bool) result ) ;
}

//----------------------------------------------------------------------
bool
EXT_MPIcommunicator:: boolean_or( bool value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: boolean_or" ) ;   
   PEL_CHECK_COLLECTIVE(true) ;
   
   int input = (int) value ;
   int result ;
   
   int mpierr = MPI_Allreduce( (void*)(&input), (void*)(&result),
                               1, MPI_INT, MPI_LOR, MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Allreduce or" ) ;
   
   return( (bool) result ) ;
}

//----------------------------------------------------------------------
double
EXT_MPIcommunicator:: sum( double value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: sum" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   double result ;
   
   int mpierr = MPI_Allreduce( (void*)(&value), (void*)(&result),
                               1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Allreduce sum" ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
double
EXT_MPIcommunicator:: max( double value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: max" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   double result ;
   
   int mpierr = MPI_Allreduce( (void*)(&value), (void*)(&result),
                               1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Allreduce max" ) ;
   
   PEL_CHECK_POST( max_POST( result, value ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
EXT_MPIcommunicator:: min( double value ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: min" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   double result ;
   
   int mpierr = MPI_Allreduce( (void*)(&value), (void*)(&result),
                               1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Allreduce min" ) ;
   
   PEL_CHECK_POST( min_POST( result, value ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: sum_vector( double* values, int nb ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: sum_vector(double*)" ) ;
   PEL_CHECK( sum_vector_PRE( values, nb ) ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   double* result = new double[nb] ;
   
   int mpierr = MPI_Allreduce( values, result,
                               nb, MPI_DOUBLE, MPI_SUM,
                               MPI_COMM_WORLD ) ;
   
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Allreduce sum_vector" ) ;

   for( size_t i=0 ; i<nb ; ++i ) values[i] = result[i] ;  
   delete [] result ;
}

//----------------------------------------------------------------------
void
EXT_MPIcommunicator:: barrier( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_MPIcommunicator:: barrier" ) ;
   PEL_CHECK_COLLECTIVE(true) ;
   
   int mpierr = MPI_Barrier( MPI_COMM_WORLD ) ;
   if( mpierr != MPI_SUCCESS ) 
      EXT_MPIcommunicator_ERROR::n0( "MPI_Barrier" ) ;
}

//internal---------------------------------------------------------------
void
EXT_MPIcommunicator_ERROR:: n0( std::string const& func )
//internal---------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** EXT_MPIcommunicator : " << endl ;
   mesg << "    call to " << func << " failed" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
