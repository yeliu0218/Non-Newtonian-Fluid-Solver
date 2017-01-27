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

#ifndef PEL_DISTRIBUTED_PARTITION_HH
#define PEL_DISTRIBUTED_PARTITION_HH

#include <PEL_Object.hh>
#include <intVector.hh>

class PEL_Communicator ;

/* 
Distribution of a set of indices on processes.

PUBLISHED
*/

class PEL_EXPORT PEL_DistributedPartition : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_DistributedPartition* create( PEL_Object* a_owner ) ;
      
   //-- Distribution built

      void set( PEL_DistributedPartition const* other ) ;

      void set_local_number( size_t a_local_number ) ;
      
      void distribute_global_number( size_t a_global_number ) ;

      void set_global_number( size_t a_global_number ) ;

   //-- Distribution comparison
     
      // Has `other' the same distributed structure than `self' ?
      bool is_compatible( PEL_DistributedPartition const* other ) const ;
      
   //-- Access

      // communicator
      PEL_Communicator const* communicator( void ) const ;
      
      // global number of items owned by all processes
      size_t global_number( void ) const ;
      
      // local number of items owned by current process  
      size_t local_number( void ) const ;
      
      /* first index of local owned items:
           `::first_local_index()' = `::start_of_partition'( `::communicator'()->rank() ) */
      size_t first_local_index( void ) const ;
      
      /* last index + 1 of local owned items:
           `::local_index_limit()' = `::first_local_index()'+`::local_number'() */
      size_t local_index_limit( void ) const ;

      /* table of partition:
           `::partitioning'( `::communicator'()->rank() ) = `::local_number'()
           `::partitioning'().sum() = `::global_number'() */
      intVector const& partitioning( void ) const ;

      // table of the first index of owned by each process
      intVector const& start_of_partition( void ) const ;

      // rank of the process which owned the index `i'
      size_t rank_of( size_t i ) const ;
      
   protected: //--------------------------------------------------------
     
   private: //----------------------------------------------------------

      PEL_DistributedPartition( PEL_Object* a_owner ) ;
     ~PEL_DistributedPartition( void ) ;
      
      PEL_DistributedPartition( void ) ;
      PEL_DistributedPartition( PEL_DistributedPartition const& other ) ;
      PEL_DistributedPartition& operator=( 
                                PEL_DistributedPartition const& other ) ;

   //-- Attributes
      
      PEL_Communicator const* const COMM ;
      
      size_t const SIZE ;
      size_t const RANK ;      

      size_t FIRST ;
      size_t LAST ;
      size_t GLOBAL_NB ;
      size_t LOCAL_NB ;
      intVector PARTITION ;
      intVector START ;
} ;

#ifndef OUTLINE
   #include <PEL_DistributedPartition.icc>
#endif

#endif



