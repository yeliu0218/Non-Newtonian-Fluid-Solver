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

#include <PEL_SequentialCommunicator.hh>

#include <PEL_assertions.hh>


PEL_SequentialCommunicator*
PEL_SequentialCommunicator:: SINGLETON = new PEL_SequentialCommunicator() ;
   
//----------------------------------------------------------------------
PEL_SequentialCommunicator:: PEL_SequentialCommunicator( void )
//----------------------------------------------------------------------
   : PEL_Communicator( "PEL_SequentialCommunicator" )
{
}

//----------------------------------------------------------------------
PEL_SequentialCommunicator:: ~PEL_SequentialCommunicator( void )
//----------------------------------------------------------------------
{
   SINGLETON = 0 ;
}

//----------------------------------------------------------------------
size_t
PEL_SequentialCommunicator:: nb_ranks( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: nb_ranks" ) ;

   static size_t result = 1 ;
   
   PEL_CHECK_POST( nb_ranks_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PEL_SequentialCommunicator:: rank( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: rank" ) ;

   static size_t result = 0 ;

   PEL_CHECK_POST( rank_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: send( size_t dest, int const* value, 
                                   int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: send( int const* )" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: send( size_t dest, double const* value, 
                                   int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: send( double const* )" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: send( size_t dest, char const* value, 
                                   int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: send( char const* )" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;
}

//----------------------------------------------------------------------
void*
PEL_SequentialCommunicator:: Isend( size_t dest, int const* value,
                                    int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: Isend( int const*)" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
void*
PEL_SequentialCommunicator:: Isend( size_t dest, double const* value,
                                    int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: Isend( double const* )" ) ;
   PEL_CHECK( send_PRE( dest, value, nb ) ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: receive( size_t src, int* value,
                                      int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: receive( int* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: receive( size_t src, double* value, 
                                      int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: receive( double* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: receive( size_t src, char* value, 
                                      int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: receive( char* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;
}

//----------------------------------------------------------------------
void*
PEL_SequentialCommunicator:: Ireceive( size_t src, int* value,
                                       int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: Ireceive( int* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
void*
PEL_SequentialCommunicator:: Ireceive( size_t src, double* value,
                                       int nb  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: Ireceive( double* )" ) ;
   PEL_CHECK( receive_PRE( src, value, nb ) ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: all_gather_v( double const* values,
                                           size_t nb,
                                           double* result,
                                           intVector const& partition,
                                           intVector const& start ) const 
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<nb ; ++i ) result[i] = values[i] ;   
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: all_gather(
                       int const* value, size_t nb, int* result  ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<nb ; ++i ) result[i] = value[i] ;   
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: all_to_all( int const* value,
                                         size_t nb,
                                         int* result  ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<nb ; ++i ) result[i] = value[i] ;   
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: gather( double const* value, size_t nb, 
                                     double* result, size_t root ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<nb ; ++i ) result[i] = value[i] ;
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: gather( int const* value, size_t nb, 
                                     int* result, size_t root ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SequentialCommunicator:: gather(int)" ) ;

   for( size_t i=0 ; i<nb ; ++i ) result[i] = value[i] ;
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: broadcast(
                              int* value,size_t nb, size_t root ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: broadcast(
                            double* value,size_t nb, size_t root ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: broadcast(
                              char* value,size_t nb, size_t root ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_SequentialCommunicator:: barrier( void ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
PEL_SequentialCommunicator:: same_value_everywhere( double val ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_SequentialCommunicator:: same_value_everywhere( int val ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}
