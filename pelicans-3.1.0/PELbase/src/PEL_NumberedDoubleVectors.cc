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

#include <PEL_NumberedDoubleVectors.hh>

#include <PEL_BalancedBinaryTree.hh>
#include <PEL_BalancedBinaryTreeIterator.hh>
#include <PEL_DoubleComparator.hh>
#include <PEL_Error.hh>

#include <iostream>
#include <sstream>

class PEL_NumberedDoubleVectorsItem : public PEL_Object
{
   public:
      
      PEL_NumberedDoubleVectorsItem( PEL_NumberedDoubleVectors* a_owner,
                                     doubleVector const& vector,
                                     size_t n ) ;
     ~PEL_NumberedDoubleVectorsItem( void ) ;

      void set( doubleVector const& vector, size_t n ) ;
      virtual bool is_equal( PEL_Object const* other ) const ;
      virtual int three_way_comparison( PEL_Object const* other ) const ;
      
      size_t N ;
      doubleVector VECTOR ;
} ;

//-----------------------------------------------------------------------------
PEL_NumberedDoubleVectors*
PEL_NumberedDoubleVectors:: create( PEL_Object* a_owner,
                                    PEL_DoubleComparator const* dbl_comp,
                                    size_t a_size_of_items )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NumberedDoubleVectors:: create" ) ;
   PEL_CHECK_PRE( dbl_comp != 0 ) ;
   
   PEL_NumberedDoubleVectors* result =
      new PEL_NumberedDoubleVectors( a_owner, dbl_comp, a_size_of_items ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->size_of_items() == a_size_of_items ) ;
   PEL_CHECK_POST( result->nb_items() == 0 ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
PEL_NumberedDoubleVectors:: PEL_NumberedDoubleVectors(
                                    PEL_Object* a_owner,
                                    PEL_DoubleComparator const* dbl_comp,
                                    size_t a_size_of_items ) 
//-----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , DBL_COMP( dbl_comp )
   , DIM( a_size_of_items )
   , NB_PTS( 0 )
   , PT_TREE( 0 )
   , PT_TREE_IT( 0 )
   , ALL_ITEMS( 0, 0 )
   , ORDER( 0 )
{
   PEL_LABEL( "PEL_NumberedDoubleVectors:: PEL_NumberedDoubleVectors" ) ;
   PT_TREE = PEL_BalancedBinaryTree::create( this ) ;
   PT_TREE_IT = PT_TREE->create_iterator( PT_TREE ) ;
   TMP = new PEL_NumberedDoubleVectorsItem( this, doubleVector( DIM ), 0 ) ;
}

//-----------------------------------------------------------------------------
PEL_NumberedDoubleVectors:: ~PEL_NumberedDoubleVectors( void ) 
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NumberedDoubleVectors:: ~PEL_NumberedDoubleVectors" ) ;
}

//-----------------------------------------------------------------------------
size_t
PEL_NumberedDoubleVectors:: nb_items( void ) const
//-----------------------------------------------------------------------------
{
   return( NB_PTS ) ;
}

//-----------------------------------------------------------------------------
size_t
PEL_NumberedDoubleVectors:: size_of_items( void ) const
//-----------------------------------------------------------------------------
{
   return( DIM ) ;
}

//-----------------------------------------------------------------------------
bool
PEL_NumberedDoubleVectors:: has( doubleVector const& item ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NumberedDoubleVectors:: has" ) ;
   PEL_CHECK_PRE( item.size()==size_of_items() ) ;

   TMP->set( item, 0 ) ;
   bool result = PT_TREE->has( TMP ) ;
   
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t
PEL_NumberedDoubleVectors:: index( doubleVector const& item ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NumberedDoubleVectors:: index" ) ;
   PEL_CHECK_PRE( item.size()==size_of_items() ) ;
   PEL_CHECK_PRE( has( item ) ) ;

   TMP->set( item, 0 ) ;
   PEL_NumberedDoubleVectorsItem const* vec =
      static_cast<PEL_NumberedDoubleVectorsItem*>( PT_TREE->item( TMP ) ) ;
   
   size_t result = vec->N ;

   PEL_CHECK_POST( result < nb_items() ) ;

   return( result ) ;
}

//-----------------------------------------------------------------------------
doubleArray2D const& 
PEL_NumberedDoubleVectors:: ordered_items( void ) const 
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NumberedDoubleVectors:: ordered_items" ) ;

   ALL_ITEMS.re_initialize( size_of_items(), nb_items() ) ;
   size_t cpt = 0 ;
   for( PT_TREE_IT->start() ; PT_TREE_IT->is_valid() ; PT_TREE_IT->go_next() )
   {
      PEL_NumberedDoubleVectorsItem const* item =
         static_cast<PEL_NumberedDoubleVectorsItem*>( PT_TREE_IT->item() ) ;
      for( size_t j=0 ; j<DIM ; j++ ) ALL_ITEMS( j, cpt ) = item->VECTOR( j ) ;
      cpt++ ;
   }
   PEL_CHECK( cpt==nb_items() ) ;
   doubleArray2D const& result = ALL_ITEMS ;
   
   PEL_CHECK_POST( result.index_bound(1) == nb_items() ) ;
   PEL_CHECK_POST( result.index_bound(0) == size_of_items() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
PEL_NumberedDoubleVectors:: order( void ) const 
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NumberedDoubleVectors:: order" ) ;

   ORDER.re_initialize( nb_items() ) ;
   size_t cpt=0 ;
   for( PT_TREE_IT->start() ; PT_TREE_IT->is_valid() ; PT_TREE_IT->go_next() )
   {
      PEL_NumberedDoubleVectorsItem const* item =
         static_cast<PEL_NumberedDoubleVectorsItem*>( PT_TREE_IT->item() ) ;
      size_t i = item->N ;
      ORDER( i ) = cpt++ ;
   }
   size_t_vector const& result = ORDER ;
   
   PEL_CHECK_POST( result.size() == nb_items() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
PEL_NumberedDoubleVectors:: extend( doubleVector const& item )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NumberedDoubleVectors:: extend" ) ;
   PEL_CHECK_PRE( item.size() == size_of_items() ) ;
   if( !has( item ) )
   {
      PEL_NumberedDoubleVectorsItem* new_one =
         new PEL_NumberedDoubleVectorsItem( this, item, NB_PTS++ ) ;
      PT_TREE->extend( new_one ) ;
      PEL_CHECK( PT_TREE->count()==NB_PTS ) ;
   }
   
   PEL_CHECK_POST( has( item ) ) ;
}

//internal--------------------------------------------------------------
PEL_NumberedDoubleVectorsItem:: PEL_NumberedDoubleVectorsItem(
                                       PEL_NumberedDoubleVectors* a_owner,
                                       doubleVector const& vector,
                                       size_t n )
//internal--------------------------------------------------------------
   : PEL_Object( a_owner )
   , N( n )
   , VECTOR( vector )
{
}

//internal--------------------------------------------------------------
PEL_NumberedDoubleVectorsItem:: ~PEL_NumberedDoubleVectorsItem( void )
//internal--------------------------------------------------------------
{
}

//internal--------------------------------------------------------------
void
PEL_NumberedDoubleVectorsItem:: set( doubleVector const& vector, size_t n )
//internal--------------------------------------------------------------
{
   VECTOR = vector ;
   N = n ;
}

//internal--------------------------------------------------------------
bool
PEL_NumberedDoubleVectorsItem:: is_equal( PEL_Object const* other ) const
//internal--------------------------------------------------------------
{
   PEL_CHECK( is_equal_PRE( other ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   bool result = three_way_comparison( other ) == 0 ;

   return( result ) ;
}

//internal--------------------------------------------------------------
int
PEL_NumberedDoubleVectorsItem:: three_way_comparison( 
                                              PEL_Object const* other ) const
//internal--------------------------------------------------------------
{
   PEL_CHECK( three_way_comparison_PRE( other ) ) ;
   PEL_CHECK( 
    dynamic_cast<PEL_NumberedDoubleVectorsItem const*>( other ) != 0 ) ;
   PEL_CHECK( VECTOR.size() == 
    static_cast<PEL_NumberedDoubleVectorsItem const*>(other)->VECTOR.size() ) ;
   PEL_CHECK( 
      dynamic_cast<PEL_NumberedDoubleVectors const*>( owner() ) != 0 ) ;

   PEL_DoubleComparator const* DBL_COMP =
      static_cast<PEL_NumberedDoubleVectors const*>( owner() )->DBL_COMP ;
   
   PEL_NumberedDoubleVectorsItem const* dbl_vec =
                   static_cast<PEL_NumberedDoubleVectorsItem const*>( other ) ;
   
   int result = 0 ;
   size_t const nb_coords = VECTOR.size() ;
   for( size_t i = 0 ; i<nb_coords && result==0 ; ++i )
   {
      result = DBL_COMP->three_way_comparison( VECTOR(i), dbl_vec->VECTOR(i) ) ;
   }

   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}


