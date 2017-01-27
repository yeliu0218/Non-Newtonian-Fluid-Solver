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

#include <GE_SetOfPoints.hh>

#include <PEL_BalancedBinaryTree.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_KeyItemPair.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_ListIterator.hh>
#include <PEL_Module.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_Root.hh>
#include <PEL_ObjectWriter.hh>
#include <PEL_Vector.hh>

#include <GE_Color.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <doubleArray2D.hh>

#include <iostream>
#include <sstream>

//-----------------------------------------------------------------------------
GE_SetOfPoints*
GE_SetOfPoints:: create( PEL_Object* a_owner,
                         size_t a_dimension )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: create" ) ;
   
   GE_SetOfPoints* result =  new GE_SetOfPoints( a_owner, a_dimension ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->dimension() == a_dimension ) ;
   PEL_CHECK_POST( result->nb_points() == 0 ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_SetOfPoints:: GE_SetOfPoints( PEL_Object* a_owner,
                                 size_t a_dimension ) 
//-----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , DIM( a_dimension )
   , NB_PTS( 0 )
   , PT_VECTOR(  PEL_Vector::create( this, 0 ) )
   , PT_COLOR_VECTOR( PEL_Vector::create( this, 0 ) )
   , PT_TREE( 0 )
   , PT_TREE_IT( 0 )
   , OBSERVERS( PEL_ListIdentity::create( this )  )
   , OBSERVERS_IT( 0 )
{
   PEL_LABEL( "GE_SetOfPoints:: GE_SetOfPoints" ) ;

   OBSERVERS_IT = OBSERVERS->create_iterator( this ) ;
   
   rebuild_points_tree() ;
}

//-----------------------------------------------------------------------------
GE_SetOfPoints:: ~GE_SetOfPoints( void ) 
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: ~GE_SetOfPoints" ) ;

   reset_points_tree() ;
}

//-----------------------------------------------------------------------------
size_t
GE_SetOfPoints:: nb_points( void ) const
//-----------------------------------------------------------------------------
{
   return( NB_PTS ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_SetOfPoints:: dimension( void ) const
//-----------------------------------------------------------------------------
{
   return( DIM ) ;
}

//-----------------------------------------------------------------------------
bool
GE_SetOfPoints:: has( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: has" ) ;
   PEL_CHECK_PRE( pt != 0 && pt->nb_coordinates()==dimension() ) ;

   bool result ;
   if( !is_ok_points_tree() )
   {
      result = PT_VECTOR->has( pt ) ;
   }
   else
   {
      result = PT_TREE->has( pt ) ;
   }
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_SetOfPoints:: index( GE_Point const* pt ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: index" ) ;
   PEL_CHECK_PRE( pt != 0 && pt->nb_coordinates()==dimension() ) ;
   PEL_CHECK_PRE( has( pt ) ) ;

   size_t result ;
   
   if( !is_ok_points_tree() )
   {
      result = PT_VECTOR->index_of( pt ) ;
   }
   else
   {
      PEL_KeyItemPair const* key_val =
         static_cast<PEL_KeyItemPair const*>( PT_TREE->item( pt ) ) ;
      result =
         (size_t)  static_cast<PEL_Int const*>( key_val->item() )->to_int() ;
   }

   PEL_CHECK_POST( result < nb_points() ) ;
   PEL_CHECK_POST( pt->is_equal( point( result ) ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_SetOfPoints:: point( size_t i ) const 
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: point" ) ;
   PEL_CHECK_PRE( i < nb_points() ) ;
   
   GE_Point const* result = static_cast<GE_Point const*>( PT_VECTOR->at(i) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_coordinates() == dimension() ) ;
   PEL_CHECK_POST( index(result) == i ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Color const*
GE_SetOfPoints:: color( size_t i ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: color" ) ;
   PEL_CHECK_PRE( i < nb_points() ) ;

   GE_Color const* result = 
                   static_cast<GE_Color const* >( PT_COLOR_VECTOR->at( i ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: start( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: start" ) ;
   PEL_CHECK( is_ok_points_tree() ) ;

   PT_TREE_IT->start() ;
}

//-----------------------------------------------------------------------------
bool
GE_SetOfPoints:: is_valid( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: is_valid" ) ;
   PEL_CHECK( is_ok_points_tree() ) ;

   return PT_TREE_IT->is_valid() ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: go_next( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: go_next" ) ;
   PEL_CHECK_PRE( is_valid() ) ;

   PT_TREE_IT->go_next() ;
}

//-----------------------------------------------------------------------------
GE_Point const*
GE_SetOfPoints:: point( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: point" ) ;
   PEL_CHECK_PRE( is_valid() ) ;

   PEL_KeyItemPair const* key_val =
                  static_cast<PEL_KeyItemPair const*>( PT_TREE_IT->item() ) ;
   GE_Point const* result =  static_cast<GE_Point const*>( key_val->key() ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_coordinates() == dimension() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_SetOfPoints:: index( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: index" ) ;
   PEL_CHECK_PRE( is_valid() ) ;

   PEL_KeyItemPair const* key_val =
      static_cast<PEL_KeyItemPair const*>( PT_TREE_IT->item() ) ;
   size_t result =
         (size_t)  static_cast<PEL_Int const*>( key_val->item() )->to_int() ;

   PEL_CHECK_POST( result < nb_points() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: extend( GE_Point const* pt )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: extend" ) ;
   PEL_CHECK_PRE( pt != 0 && pt->nb_coordinates()==dimension() ) ;
   PEL_SAVEOLD( bool, has, has(pt) ) ;
   PEL_SAVEOLD( size_t, nb_points, nb_points() ) ;

   if( !has( pt ) )
   {
      append( pt, GE_Color::null_color() ) ;
   }
   
   PEL_CHECK_POST( has( pt ) ) ;
   PEL_CHECK_POST( IMPLIES(  OLD(has), nb_points()==OLD(nb_points) ) ) ;
   PEL_CHECK_POST( IMPLIES( !OLD(has), nb_points()==OLD(nb_points)+1 ) ) ;
   PEL_CHECK_POST( IMPLIES( !OLD(has), index(pt)==OLD(nb_points) ) ) ;
   PEL_CHECK_POST( point( index(pt) )->owner()==this ) ;
   PEL_CHECK_POST( pt->is_equal( point( index(pt) ) ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: append( GE_Point const* pt, GE_Color const* pt_color )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: append" ) ;
   PEL_CHECK_PRE( pt != 0  && pt->nb_coordinates()==dimension() ) ;
   PEL_CHECK_PRE( !has( pt ) ) ;
   PEL_SAVEOLD( size_t, nb_points, nb_points() ) ;

   // Point to append :
   GE_Point* new_pt = pt->create_clone( this ) ;
   size_t const new_index = NB_PTS ;
   NB_PTS++ ;

   // Vector :
   if( new_index >= PT_VECTOR->index_limit() )
   {
      PT_VECTOR->resize( new_index+1 ) ;
      PT_COLOR_VECTOR->resize( new_index+1 ) ;
   }
   PT_VECTOR->set_at( new_index, new_pt ) ;
   if( pt_color!=0 )
   {
      PT_COLOR_VECTOR->set_at( new_index,
                          const_cast<GE_Color*>( pt_color ) ) ;
   }
   else
   {
      PT_COLOR_VECTOR->set_at( new_index,
                          const_cast<GE_Color*>( GE_Color::null_color() ) ) ;
   }
   
   // Binary tree :
   if( is_ok_points_tree() )
   {
      add_vertex_in_points_tree( new_index, new_pt ) ;
   }

   PEL_CHECK_POST( has( pt ) ) ;
   PEL_CHECK_POST( IMPLIES( ( pt_color == 0 ),
                      color( index( pt ) ) == GE_Color::null_color() ) ) ;
   PEL_CHECK_POST( IMPLIES( ( pt_color != 0 ),
                      color( index( pt ) ) == pt_color ) ) ;
   PEL_CHECK_POST( nb_points() == OLD( nb_points ) + 1 ) ;
   PEL_CHECK_POST( index( pt ) == OLD( nb_points ) ) ;
   PEL_CHECK_POST( point( index(pt) )->owner() == this ) ;
   PEL_CHECK_POST( pt->is_equal( point( index( pt ) ) ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: insert_at( size_t i, GE_Point const* pt, 
                            GE_Color const* pt_color )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: insert_at" ) ;
   PEL_CHECK_PRE( i <= nb_points() ) ;
   PEL_CHECK_PRE( pt != 0  && pt->nb_coordinates()==dimension() ) ;
   PEL_CHECK_PRE( !has( pt ) ) ;
   PEL_SAVEOLD( size_t, nb_points, nb_points() ) ;

   if( i==nb_points() )
   {
      append( pt, pt_color ) ;
   }
   else
   {
      // Vector :
      if( NB_PTS >= PT_VECTOR->index_limit() )
      {
	 PT_VECTOR->resize( NB_PTS+1 ) ;
	 PT_COLOR_VECTOR->resize( NB_PTS+1 ) ;
      }

      for( size_t j=NB_PTS; j>i; j-- )
      {
	 PT_VECTOR->set_at( j, PT_VECTOR->at( j-1 ) ) ; 
	 PT_COLOR_VECTOR->set_at( j, PT_COLOR_VECTOR->at( j-1 ) ) ; 
      }
      NB_PTS++ ;
 
      // Point to append :
      GE_Point* new_pt = pt->create_clone( this ) ;

      PT_VECTOR->set_at( i, new_pt ) ;

      if( pt_color!=0 )
      {
	 PT_COLOR_VECTOR->set_at( i,
			   const_cast<GE_Color*>( pt_color ) ) ;
      }
      else
      {
	 PT_COLOR_VECTOR->set_at( i,
			   const_cast<GE_Color*>( GE_Color::null_color() ) ) ;
      }
   
      // Binary tree :
      if( is_ok_points_tree() )
      {
	 reset_points_tree() ;
      }
   }

   PEL_CHECK_POST( has( pt ) ) ;
   PEL_CHECK_POST( IMPLIES( ( pt_color == 0 ),
                      color( index( pt ) ) == GE_Color::null_color() ) ) ;
   PEL_CHECK_POST( IMPLIES( ( pt_color != 0 ),
                      color( index( pt ) ) == pt_color ) ) ;
   PEL_CHECK_POST( nb_points() == OLD( nb_points ) + 1 ) ;
   PEL_CHECK_POST( index( pt ) == i ) ;
   PEL_CHECK_POST( point( index( pt ) )->owner() == this ) ;
   PEL_CHECK_POST( pt->is_equal( point( index( pt ) ) ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: remove_at( size_t i )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: remove_at" ) ;
   PEL_CHECK_PRE( i < nb_points() ) ;
   PEL_SAVEOLD( size_t, nb_points, nb_points() ) ;

   PEL_Object* old_pt = PT_VECTOR->at( i ) ;

   destroy_possession( old_pt) ;

   // Vector :
   PT_VECTOR->remove_at( i ) ;
   PT_COLOR_VECTOR->remove_at( i ) ;
   NB_PTS-- ;

   // Binary tree :
   if( is_ok_points_tree() )
   {
      reset_points_tree() ;
   }

   PEL_CHECK_POST( nb_points() == OLD( nb_points ) - 1 ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: move_and_update( PEL_Vector const* dM_table ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: move_and_update" ) ;
   PEL_CHECK_PRE( dM_table!=0 &&
                  dM_table->index_limit() == nb_points() ) ;
   PEL_CHECK_PRE(
      FORALL( ( size_t i=0 ; i<nb_points() ; ++i ),
         dynamic_cast<GE_Vector const*>( dM_table->at(i) )!=0 &&
         static_cast<GE_Vector const*>
                   ( dM_table->at(i) )->nb_components()==dimension() ) ) ;
   
   // Move :
   for( size_t i=0 ; i<nb_points() ; ++i )
   {
      GE_Vector const* dm = static_cast<GE_Vector const*>( dM_table->at(i) ) ;
      PEL_CHECK( dynamic_cast<GE_Point*>( PT_VECTOR->at(i) )!=0 ) ;
      GE_Point* pt = static_cast<GE_Point*>( PT_VECTOR->at(i) ) ;
      pt->translate( 1., dm ) ;
   }
   reset_points_tree() ;
   rebuild_points_tree() ;
   
   // And update :
   update_observers() ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: take_position_and_update( 
                                        PEL_Vector const* new_positions ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: take_position_and_update" ) ;
   PEL_CHECK_PRE( new_positions!=0 &&
                  new_positions->index_limit() == nb_points() ) ;
   PEL_CHECK_PRE(
      FORALL( ( size_t i=0 ; i<nb_points() ; ++i ),
         dynamic_cast<GE_Point const*>( new_positions->at(i) )!=0 &&
         static_cast<GE_Point const*>
                  ( new_positions->at(i) )->nb_coordinates()==dimension() ) ) ;
   
   // Move :
   for( size_t i=0 ; i<nb_points() ; ++i )
   {
      GE_Point const* new_pt = 
                      static_cast<GE_Point const*>( new_positions->at(i) ) ;
      PEL_CHECK( dynamic_cast<GE_Point*>( PT_VECTOR->at(i) )!=0 ) ;
      GE_Point* pt = static_cast<GE_Point*>( PT_VECTOR->at(i) ) ;
      pt->set( new_pt ) ;
   }
   reset_points_tree() ;
   rebuild_points_tree() ;
   
   // And update :
   update_observers() ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: modify_color( size_t i, GE_Color const* pt_color )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: modify_color" ) ;
   PEL_CHECK_PRE( i<nb_points() ) ;
   PEL_CHECK_PRE( pt_color!=0 ) ;

   PT_COLOR_VECTOR->set_at( i, const_cast<GE_Color*>( pt_color ) ) ;
   
   PEL_CHECK_POST( color(i)==pt_color ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: attach_observer( PEL_Object* observer )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: attach_observer" ) ;
   PEL_CHECK_PRE( observer!=0 ) ;
   PEL_CHECK_PRE( !has_as_observer( observer ) ) ;

   OBSERVERS->append( observer ) ;

   PEL_CHECK_POST( has_as_observer( observer ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: detach_observer( PEL_Object const* observer )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: detach_observer" ) ;
   PEL_CHECK_PRE( observer!=0 ) ;
   PEL_CHECK_PRE( has_as_observer( observer ) ) ;
   
   OBSERVERS->remove( observer ) ;
  
   PEL_CHECK_POST( !has_as_observer( observer ) ) ;
}

//-----------------------------------------------------------------------------
bool
GE_SetOfPoints:: has_as_observer( PEL_Object const* observer ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: has_as_observer" ) ;
   PEL_CHECK_PRE( observer!=0 ) ;

   return( OBSERVERS->has( observer ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: save_state( PEL_ObjectWriter* writer ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: save_state" ) ;
   PEL_CHECK_PRE( save_state_PRE( writer ) ) ;

   writer->start_new_object( "GE_SetOfPoints" ) ;

   writer->add_entry( "dimension", PEL_Int::create( 0, dimension() ) ) ;
   writer->add_entry( "nb_points", PEL_Int::create( 0, nb_points() ) ) ;
   
   // Coordinates :
   doubleArray2D coordinates( nb_points(), dimension() ) ;
   for( size_t i=0 ; i<nb_points() ; ++i )
   {
      GE_Point const* pt = static_cast<GE_Point const*>( PT_VECTOR->at(i) ) ;
      for( size_t d=0 ; d<dimension() ; ++d )
      {
         coordinates( i, d ) = pt->coordinate(d) ;
      }
   }
   writer->add_entry( "points_coordinates",
                      PEL_DoubleArray2D:: create( 0, coordinates ) ) ;

   writer->finalize_object() ;   

   PEL_CHECK_POST( save_state_POST( writer ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: restore_state( PEL_ObjectReader* reader )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: restore_state" ) ;
   PEL_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "GE_SetOfPoints" ) ;

   // Does some checks :
   PEL_ASSERT( reader->data_of_entry( "dimension" )->to_int() == 
               (int) dimension() ) ;
   PEL_ASSERT( reader->data_of_entry( "nb_points" )->to_int() ==
               (int) nb_points() ) ;

   // Retrieving coordinates :
   doubleArray2D const& coordinates =
         reader->data_of_entry( "points_coordinates" )->to_double_array2D() ;
   
   reader->end_object_retrieval() ;

   for( size_t i=0 ; i<nb_points() ; ++i )
   {
      GE_Point* pt = static_cast<GE_Point*>( PT_VECTOR->at( i ) ) ;
      for( size_t d=0 ; d<dimension() ; ++d )
      {
         pt->set_coordinate( d, coordinates( i, d ) ) ;
      }
   }
   reset_points_tree() ;
   rebuild_points_tree() ;

   // Update :
   update_observers_for_restore_state() ;

   PEL_CHECK_POST( restore_state_POST( reader ) ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: print( std::ostream& os, size_t indent_width ) const 
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << "GE_SetOfPoints :" << std::endl ;
   for( size_t i=0 ; i<nb_points() ; ++i )
   {
      os << space << "   n°" << i << " : " ;
      point(i)->print( os, 0 ) ;
      os << "   " ;
      color(i)->print( os, 0 ) ;
      os << std::endl ;
   }
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: update_observers( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: update_observers" ) ;

   OBSERVERS_IT->start() ;
   for( ; OBSERVERS_IT->is_valid() ; OBSERVERS_IT->go_next() )
   {
      PEL_CHECK( OBSERVERS_IT->item() != 0 ) ;
      OBSERVERS_IT->item()->update() ;
   }
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: update_observers_for_restore_state( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: update_observers_for_restore_state" ) ;

   OBSERVERS_IT->start() ;
   for( ; OBSERVERS_IT->is_valid() ; OBSERVERS_IT->go_next() )
   {
      PEL_CHECK( OBSERVERS_IT->item() != 0 ) ;
      OBSERVERS_IT->item()->update_for_restore_state() ;
   }
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: reset_points_tree( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: reset_points_tree" ) ;
   if( is_ok_points_tree() )
   {
      PT_TREE->destroy() ;
      PT_TREE = 0 ;
   }
   PEL_CHECK( !is_ok_points_tree() ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: rebuild_points_tree( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: rebuild_points_tree" ) ;
   PEL_CHECK( !is_ok_points_tree() ) ;

   PT_TREE = PEL_BalancedBinaryTree::create( 0 ) ;
   PT_TREE_IT = PT_TREE->create_iterator( PT_TREE ) ;
   
   for( size_t i=0 ; i<nb_points() ; ++i )
   {
      GE_Point* pt = static_cast<GE_Point*>( PT_VECTOR->at( i ) ) ;
      add_vertex_in_points_tree( i, pt ) ;
   }
   PEL_CHECK( is_ok_points_tree() && PT_TREE->count()==nb_points() ) ;
}

//-----------------------------------------------------------------------------
void
GE_SetOfPoints:: add_vertex_in_points_tree( size_t const pt_index,
                                            GE_Point* pt ) const 
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: add_vertex_in_points_tree" ) ;
   PEL_CHECK( is_ok_points_tree() ) ;
   PEL_CHECK( pt_index < nb_points() ) ;

   if( PT_TREE->has( pt ) )
   {
      PEL_KeyItemPair const* key_val =
         static_cast<PEL_KeyItemPair const*>( PT_TREE->item( pt ) ) ;
      size_t old_index =
         (size_t)  static_cast<PEL_Int const*>( key_val->item() )->to_int() ;

      std::ostringstream message ;
      message << std::endl ;
      message << "Invalid set of points : " << std::endl
              << "                      index : " << pt_index  << std::endl
              << "                  and index : " << old_index << std::endl
              << "    refer to the same point : " ;
      pt->print( message, 0 ) ;
      message << std::endl ;
      PEL_Error::object()->raise_plain( message.str() ) ;
   }
   PEL_Int* i_pt = PEL_Int::create( PT_TREE, (int) pt_index ) ;
   PEL_KeyItemPair* key_val = PEL_KeyItemPair::create( PT_TREE, pt, i_pt ) ;
   PT_TREE->extend( key_val ) ;
   PEL_CHECK( has( pt ) && index( pt )==pt_index ) ;
}

//-----------------------------------------------------------------------------
bool
GE_SetOfPoints:: is_ok_points_tree( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SetOfPoints:: is_ok_points_tree" ) ;
   return( PT_TREE!=0 ) ;
}


