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

#include <PEL_DistributedPartition.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_Exec.hh>

#include <doubleVector.hh>

#ifdef OUTLINE
   #define inline
   #include <PEL_DistributedPartition.icc>
   #undef inline
#endif

//----------------------------------------------------------------------
PEL_DistributedPartition*
PEL_DistributedPartition:: create( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DistributedPartition:: create" ) ;

   PEL_DistributedPartition* result =
                               new PEL_DistributedPartition( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->communicator() == PEL_Exec::communicator() ) ;
   PEL_CHECK_POST( result->global_number() == 0 ) ;
   PEL_CHECK_POST( result->local_number() == 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition:: PEL_DistributedPartition( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , COMM( PEL_Exec::communicator() )
   , SIZE( PEL_Exec::communicator()->nb_ranks() )
   , RANK( PEL_Exec::communicator()->rank() )
   , FIRST( 0 )
   , LAST( 0 )
   , GLOBAL_NB( 0 )
   , LOCAL_NB( 0 )
   , PARTITION( 0 )
   , START( 0 )
{
   PEL_LABEL( "PEL_DistributedPartition:: PEL_DistributedPartition" ) ;

   PARTITION.re_initialize( SIZE ) ;
   START.re_initialize( SIZE ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition:: ~PEL_DistributedPartition( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_DistributedPartition:: set( PEL_DistributedPartition const* other ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DistributedPartition:: set" ) ;
   PEL_CHECK_PRE( other!=0 ) ;

   PARTITION = other->PARTITION ;
   START = other->START ;
   FIRST = other->FIRST ;
   LAST = other->LAST ;
   LOCAL_NB = other->LOCAL_NB ;
   GLOBAL_NB = other->GLOBAL_NB ;

   PEL_CHECK_POST( first_local_index() == other->first_local_index() ) ;
   PEL_CHECK_POST( local_index_limit() == other->local_index_limit() ) ;
   PEL_CHECK_POST( global_number() == other->global_number() ) ;
   PEL_CHECK_POST( local_number() == other->local_number() ) ;
   PEL_CHECK_POST( partitioning() == other->partitioning() ) ;
   PEL_CHECK_POST( start_of_partition() ==  other->start_of_partition() ) ;
}

//----------------------------------------------------------------------
bool
PEL_DistributedPartition:: is_compatible( 
                              PEL_DistributedPartition const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DistributedPartition:: is_compatible" ) ;
   PEL_CHECK_PRE( other!=0 ) ;
   
   bool result = ( PARTITION == other->PARTITION ) ;

   PEL_CHECK_POST(
             EQUIVALENT( result, other->partitioning()==partitioning() ) ) ; 
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_DistributedPartition:: set_local_number( size_t a_local_number )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DistributedPartition:: set_local_number" ) ;
   PEL_CHECK_COLLECTIVE( true ) ;

   COMM->all_gather( a_local_number, PARTITION ) ;
   GLOBAL_NB = PARTITION.sum() ;
   START( 0 ) = 0 ;
   for( size_t i=1 ; i<SIZE ; ++i )
   {
      START( i ) = START( i-1 ) + PARTITION( i-1 ) ;
   }   
   FIRST = START( RANK ) ;
   LAST = FIRST + PARTITION( RANK ) ;
   LOCAL_NB = PARTITION( RANK ) ;

   PEL_CHECK_POST( local_number() == a_local_number ) ;
   PEL_CHECK_POST( global_number() ==
                   (  size_t) communicator()->sum( (double) a_local_number ) ) ;
   PEL_CHECK_POST( partitioning()( communicator()->rank() ) == (int) a_local_number ) ;
   PEL_CHECK_POST( partitioning().sum() == (int) global_number() ) ;
}

//----------------------------------------------------------------------
void
PEL_DistributedPartition:: distribute_global_number( size_t a_global_number )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DistributedPartition:: distribute_global_number" ) ;
   PEL_CHECK_COLLECTIVE( false ) ;

   GLOBAL_NB = a_global_number ;
   size_t const dim = GLOBAL_NB / SIZE ;
   size_t const r   = GLOBAL_NB % SIZE ;
   for( size_t i=0 ; i<SIZE ; ++i )
   {
      PARTITION(i) = ( i<r ? dim+1 : dim ) ;
   }
   START(0) = 0 ;
   for( size_t i=1 ; i<SIZE ; ++i )
   {
      START(i) = START(i-1)+PARTITION(i-1) ;
   }   
   FIRST = START( RANK ) ;
   LAST = FIRST + PARTITION( RANK ) ;
   LOCAL_NB = PARTITION( RANK ) ;

   PEL_CHECK_POST( global_number() == a_global_number ) ;
   PEL_CHECK_POST( partitioning().sum() == (int) global_number() ) ;
   PEL_CHECK_POST(
      local_number() == a_global_number/communicator()->nb_ranks() ||
      local_number() == a_global_number/communicator()->nb_ranks()+1 ) ;
}

//----------------------------------------------------------------------
void
PEL_DistributedPartition:: set_global_number( size_t a_global_number )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DistributedPartition:: set_global_number" ) ;
   PEL_CHECK_COLLECTIVE( false ) ;

   GLOBAL_NB = a_global_number ;
   PARTITION.set( 0 ) ;
   PARTITION( RANK ) = GLOBAL_NB ;
   START(0) = 0 ;
   for( size_t i=1 ; i<SIZE ; ++i )
   {
      START(i) = START(i-1)+PARTITION(i-1) ;
   }   
   FIRST = 0 ;
   LAST =  GLOBAL_NB ;
   LOCAL_NB = GLOBAL_NB ;

   PEL_CHECK_POST( global_number() == a_global_number ) ;
   PEL_CHECK_POST( partitioning().sum() == (int) global_number() ) ;
   PEL_CHECK_POST( partitioning()( communicator()->rank() ) == a_global_number ) ;
   PEL_CHECK_POST( local_number() == a_global_number ) ;
   PEL_CHECK_POST( first_local_index() == 0 ) ;
   PEL_CHECK_POST( local_index_limit() == a_global_number ) ;
}
