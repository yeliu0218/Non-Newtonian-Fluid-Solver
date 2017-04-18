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

#ifndef GE_SET_OF_POINTS_HH
#define GE_SET_OF_POINTS_HH

#include <PEL_Object.hh>

class PEL_BalancedBinaryTree ;
class PEL_BalancedBinaryTreeIterator ;
class PEL_ListIdentity ;
class PEL_ListIterator ;
class PEL_Vector ;

class GE_Color ;
class GE_Point ;

/*
Sets of colorized points.

A set of point is composed of several `GE_Point::' objects of the same number
of coordinates, each point beeing colorized with a `GE_Color::' object.

The points might be moved, in which case `self' is reponsible for updating
all attached observers.

PUBLISHED
*/

class PEL_EXPORT GE_SetOfPoints : public PEL_Object
{
   public: //------------------------------------------------------------------
   
   //-- Instance delivery and initialization

      // Create, initialize and return an instance devoted to store points
      // with `a_dimension' nunmber of components.
      static GE_SetOfPoints* create( PEL_Object* a_owner,
                                     size_t a_dimension ) ;

   //-- Access

      // number of points
      size_t nb_points( void ) const ;

      // number of coordinates of the points (same for all of them)
      size_t dimension( void ) const ;

      // Does `self' contains a point comparing equal to `pt' ?
      bool has( GE_Point const* pt ) const ;

   //-- Access with indices(51)

      // index of the point matching `pt'
      size_t index( GE_Point const* pt ) const ;

      // `i'-th point
      GE_Point const* point( size_t i ) const ;    

      // color of the `i'-th point
      GE_Color const* color( size_t i ) const ;

   //-- Access with an iteration independent of the introduction order(52)
      
      // Move point-iterator to the first position (independent of the
      // introduction order).
      void start( void ) ;

      // Is point-iterator position valid ?
      bool is_valid( void ) const ;

      // Move point iterator one position (independent of the
      // introduction order).
      void go_next( void ) ;

      // point at current iterator position
      GE_Point const* point( void ) const ;
      
      // index of the point at current iterator position
      size_t index( void ) const ;
      
   //-- Element change

      // Ensure that self includes a point equal to `pt'
      // (if necessary, create a new point with `self' as owner).
      void extend( GE_Point const* pt ) ;

      // Append a new point equal to `pt' with `self' as owner.
      void append( GE_Point const* pt,
                   GE_Color const* pt_color = 0 ) ;

      // Increase by one the index of all items starting from the `i'-th 
      // position and add a new point equal to `pt' to the `i'-th position 
      // with `self' as owner.
      void insert_at( size_t i,
                      GE_Point const* pt,
                      GE_Color const* pt_color = 0 ) ;

      // Remove the point of index `i'.
      void remove_at( size_t i ) ;

      // Translate all the points according to the vectors of `dM_table',
      // and update the polyhedrons attached to `self'.
      void move_and_update( PEL_Vector const* dM_table ) const ;

      // Give to all the points the position specified by new_positions,
      // and update the polyhedrons attached to `self'.
      void take_position_and_update( PEL_Vector const* new_positions) const ;

      // Modify the color of the point of index `i'.
      void modify_color( size_t i, GE_Color const* pt_color ) ;

   //-- Observer(800)

      // Make `observer' subscribe to receive notification when `self' changes
      // its state.
      void attach_observer( PEL_Object* observer ) ;

      // Remove `observer' from the subscription list.
      void detach_observer( PEL_Object const* observer ) ;

      // Does `observer' belong to the subscribtion list ?
      bool has_as_observer( PEL_Object const* observer ) const ;

   //-- Persistence   

      virtual void save_state( PEL_ObjectWriter* writer ) const ;

      virtual void restore_state( PEL_ObjectReader* reader ) ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //---------------------------------------------------------------

   private: //-----------------------------------------------------------------

      GE_SetOfPoints( void ) ;
     ~GE_SetOfPoints( void ) ;
      GE_SetOfPoints( GE_SetOfPoints const& other ) ;
      GE_SetOfPoints& operator=( GE_SetOfPoints const& other ) ;
      
      GE_SetOfPoints( PEL_Object* a_owner, size_t a_dimension ) ;

   //-- Observer
      
      void update_observers( void ) const ;
      void update_observers_for_restore_state( void ) const ;

   //-- Points tree

      void reset_points_tree( void ) const ;

      void rebuild_points_tree( void ) const ;

      void add_vertex_in_points_tree( size_t const pt_index,
                                      GE_Point* pt ) const ;

      bool is_ok_points_tree( void ) const ;

   //-- Attributes

      // default increment for resizing the vector
      static size_t const BUCKET_SIZE ;

      // number of components of the points
      size_t const DIM ;

      // number of points
      size_t NB_PTS ;

      // optimized storage for fast searching from global number
      PEL_Vector* PT_VECTOR ;
      PEL_Vector* PT_COLOR_VECTOR ;

      // optimized storage for fast searching from GE_Point
      mutable PEL_BalancedBinaryTree* PT_TREE ;
      mutable PEL_BalancedBinaryTreeIterator* PT_TREE_IT ;

      // observer
      PEL_ListIdentity* OBSERVERS ;
      PEL_ListIterator* OBSERVERS_IT ;
} ;

#endif
