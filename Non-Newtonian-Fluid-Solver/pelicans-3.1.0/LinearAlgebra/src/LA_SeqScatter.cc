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

#include <LA_SeqScatter.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_DistributedPartition.hh>

#include <LA_SeqImplementation.hh>
#include <LA_SeqVector.hh>

//----------------------------------------------------------------------
LA_SeqScatter*
LA_SeqScatter:: create( PEL_Object* a_owner,
                        size_t a_nb_rows,
                        size_t_vector const& a_repatriated_items_table,
                        size_t_vector const& a_local_indices_table )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqScatter:: create" ) ;
   PEL_CHECK_PRE( a_repatriated_items_table.size() ==
                                        a_local_indices_table.size() ) ;
   PEL_CHECK_PRE(
      FORALL(
         ( size_t i=0 ; i<a_repatriated_items_table.size() ; ++i ),
         a_repatriated_items_table(i)<a_nb_rows ) ) ;

   LA_SeqScatter* result =
      new LA_SeqScatter( a_owner,
                         a_nb_rows,
                         a_repatriated_items_table, a_local_indices_table ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->size() == a_repatriated_items_table.size() ) ;
   PEL_CHECK_POST( result->repatriated_items() == a_repatriated_items_table ) ;
   PEL_CHECK_POST( result->local_indices() == a_local_indices_table ) ;
   PEL_CHECK_POST( result->distribution()->global_number() == a_nb_rows ) ;
   PEL_CHECK_POST( result->distribution()->local_number() == a_nb_rows ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqScatter:: LA_SeqScatter(
                     PEL_Object* a_owner,
                     size_t a_nb_rows,
                     size_t_vector const& a_repatriated_items_table,
                     size_t_vector const& a_local_indices_table )
//----------------------------------------------------------------------
   : LA_Scatter( a_owner )
   , NB_ROWS( a_nb_rows )
   , NEEDED( a_repatriated_items_table )
   , LOCAL( a_local_indices_table )
   , DIST( 0 )
{
}

//----------------------------------------------------------------------
LA_SeqScatter:: ~LA_SeqScatter( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_SeqScatter:: implementation( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqScatter:: implementation" ) ;
   
   static LA_Implementation const* result = LA_SeqImplementation::object() ;

   PEL_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_SeqScatter:: size( void ) const
//----------------------------------------------------------------------
{
   return( NEEDED.size() ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
LA_SeqScatter:: repatriated_items( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqScatter:: repatriated_items" ) ;
   
   size_t_vector const& result = NEEDED ;

   PEL_CHECK_POST( repatriated_items_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
LA_SeqScatter:: local_indices( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqScatter:: local_indices" ) ;
   
   size_t_vector const& result = LOCAL ;

   PEL_CHECK_POST( local_indices_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition const*
LA_SeqScatter:: distribution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqScatter:: distribution" ) ;

   if( DIST == 0 )
   {
      DIST = PEL_DistributedPartition::create(
                                const_cast<LA_SeqScatter*>( this ) ) ;
      DIST->set_global_number( NB_ROWS ) ;
   }
   
   PEL_DistributedPartition const* result = DIST ;

   PEL_CHECK_POST( distribution_POST( result ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
void
LA_SeqScatter:: get( LA_Vector const* source,
                     LA_SeqVector* dest ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqScatter:: get" ) ;
   PEL_CHECK_PRE( get_PRE( source, dest) ) ;

   LA_SeqVector const* bsource = static_cast<LA_SeqVector const*>( source ) ;
   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( source ) != 0 ) ;

   size_t n = size() ;
   
   for( size_t i=0 ; i<n ; ++i )
   {
      dest->set_item( LOCAL(i), bsource->item( NEEDED(i) ) ) ;
   }

   dest->synchronize() ;
   PEL_CHECK_POST( get_POST( source, dest) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqScatter:: set( LA_SeqVector const* source,
                     LA_Vector* dest ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_SeqScatter:: set" ) ;
   PEL_CHECK_PRE( set_PRE( source, dest) ) ;

   PEL_CHECK( dynamic_cast<LA_SeqVector const*>( dest ) != 0 ) ;

   size_t n = size() ;
   
   for( size_t i=0 ; i<n ; ++i )
   {
      dest->set_item( NEEDED(i), source->item( LOCAL(i) ) ) ;
   }
   
   dest->synchronize() ;
   PEL_CHECK_POST( set_POST( source, dest) ) ;

}

//----------------------------------------------------------------------
bool
LA_SeqScatter:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Scatter::implementation_POST( result ) ) ;
   PEL_ASSERT( result == LA_SeqImplementation::object() ) ;
   return( true ) ;
}
