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

#include <LA_Scatter.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_DistributedPartition.hh>

#include <LA_SeqVector.hh>

#include <size_t_vector.hh>

//----------------------------------------------------------------------
LA_Scatter:: LA_Scatter( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
{
}

//----------------------------------------------------------------------
LA_Scatter:: ~LA_Scatter( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_Scatter:: implementation_POST( LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: repatriated_items_POST( size_t_vector const& result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result.size() == size() ) ;
   PEL_ASSERT(
      FORALL(
         ( size_t i=0 ; i<size() ; ++i ),
               result(i)<distribution()->global_number() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: local_indices_POST( size_t_vector const& result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result.size() == size() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: distribution_POST( PEL_DistributedPartition const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == this ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: get_PRE( LA_Vector const* source,
                      LA_SeqVector const* dest ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( source != 0 ) ;
   PEL_ASSERT( source->is_synchronized() ) ;
   PEL_ASSERT( source->implementation() == implementation() ) ;
   PEL_ASSERT( source->row_distribution()->is_compatible( distribution() ) ) ;
   PEL_ASSERT( dest!=0 ) ;
   PEL_ASSERT(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              repatriated_items()(i) < source->nb_rows() ) ) ;
   PEL_ASSERT(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              local_indices()(i) < dest->nb_rows() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: get_POST( LA_Vector const* source,
                       LA_SeqVector const* dest ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( dest->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: set_PRE( LA_SeqVector const* source,
                      LA_Vector const* dest ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( source!=0 ) ;
   PEL_ASSERT( source->is_synchronized() ) ;
   PEL_ASSERT( dest!=0 ) ;
   PEL_ASSERT( dest->implementation() == implementation() ) ;
   PEL_ASSERT( dest->row_distribution()->is_compatible( distribution() ) ) ;
   PEL_ASSERT(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              local_indices()(i) < source->nb_rows() ) ) ;
   PEL_ASSERT(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              repatriated_items()(i) < dest->nb_rows() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: set_POST( LA_SeqVector const* source,
                       LA_Vector const* dest ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( dest->is_synchronized() ) ;
   return( true ) ;
}


