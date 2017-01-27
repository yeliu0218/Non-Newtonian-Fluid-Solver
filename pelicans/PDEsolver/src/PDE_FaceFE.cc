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

#include <PDE_FaceFE.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Vector.hh>
#include <size_t_vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Transform.hh>

#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>

#include <iostream>

using std::cout ; using std::endl ;
using std::string ;

//------------------------------------------------------------------------
PDE_FaceFE*
PDE_FaceFE:: create( PEL_Object* a_owner,
                     size_t a_number,
                     GE_Mpolyhedron* a_polyhedron,
                     GE_Color const* a_color,
                     size_t a_refinement_level )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: create" ) ;

   PEL_CHECK_PRE( a_polyhedron!=0 ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;
   
   PDE_FaceFE* result = new PDE_FaceFE( a_owner, 
                                        a_number, 
                                        a_polyhedron, 
                                        a_color,
                                        a_refinement_level ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->id_number()==a_number ) ;
   PEL_CHECK_POST( result->polyhedron()==a_polyhedron ) ;
   PEL_CHECK_POST( result->color()==a_color ) ;
   PEL_CHECK_POST( result->nb_adjacent_cells()==0 ) ;
   
   return( result ) ;
}

//------------------------------------------------------------------------
PDE_FaceFE:: PDE_FaceFE( PEL_Object* a_owner,
                         size_t a_number,
                         GE_Mpolyhedron* a_polyhedron,
                         GE_Color const* a_color,
                         size_t a_refinement_level )
//-------------------------------------------------------------------------
   : PDE_MeshFE( a_owner, a_number, a_polyhedron, a_color, a_refinement_level )
   , CELL_1( 0 )
   , CELL_2( 0 )
   , bound( 0 )
   , SIDE_ID( PEL::bad_index() )
   , PERIODIC_NEIGHBOUR( 0 )
   , PERIODIC_TRANSFORM( 0 )
   , PARENT( 0 )
   , CHILDS( 0 )
{
}

//------------------------------------------------------------------------
PDE_FaceFE:: ~PDE_FaceFE( void )
//------------------------------------------------------------------------
{
}

//------------------------------------------------------------------------
void
PDE_FaceFE:: insert_adjacent_cell( PDE_CellFE* cell )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: insert_adjacent_cell" ) ;
   PEL_CHECK_PRE( cell != 0 ) ;
   PEL_CHECK_PRE( cell->refinement_level() <= refinement_level() ) ;
   
   //???? preconditions
   PEL_ASSERT( CELL_1!=cell && CELL_2!=cell ) ;
   PEL_ASSERT( CELL_2==0 || (CELL_1==0 && CELL_2==0 ) ) ;

   if( CELL_1 == 0 )
   {
      CELL_1 = cell ;
   }
   else if( CELL_2 == 0 )
   {
      CELL_2 = cell ;
   }
}

//------------------------------------------------------------------------
PDE_BasisFunction*
PDE_FaceFE:: basis_function( size_t ee, size_t ln ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: basis_function" ) ;
   PEL_CHECK_PRE( basis_function_PRE( ee, ln ) ) ;

   PDE_BasisFunction* result = 0 ;

   PEL_CHECK_POST( basis_function_POST( result, ee, ln ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_FaceFE:: replace_adjacent_cell( PDE_CellFE* old_cell,
                                    PDE_CellFE* new_cell )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: replace_adjacent_cell" ) ;
   PEL_CHECK_PRE( old_cell != 0 && new_cell != 0 ) ;
   PEL_CHECK_PRE( old_cell != new_cell ) ;
   PEL_CHECK_PRE( old_cell == new_cell->parent() || 
                  new_cell == old_cell->parent() ) ;
   PEL_CHECK_PRE( 
      EXISTS( ( size_t i=0 ; i<nb_adjacent_cells() ; ++i ),
         adjacent_cell(i) == old_cell ) ) ;
   PEL_CHECK_PRE( new_cell->refinement_level() <= refinement_level() ) ;
   PEL_CHECK_PRE( refinement_level() == new_cell->refinement_level() ||
                  refinement_level() == 
                   adjacent_cell_other_than( old_cell )->refinement_level() ) ;

   //??????
   PEL_ASSERT( CELL_1==old_cell || CELL_2==old_cell ) ;
   
   if( CELL_1 == old_cell )
   {
      CELL_1 = new_cell ;
   }
   else
   {
      PEL_CHECK( CELL_2 == old_cell ) ;
      CELL_2 = new_cell ;
   }
}

//----------------------------------------------------------------------
void
PDE_FaceFE:: set_adjacent_cells( PDE_CellFE* cell_a, PDE_CellFE* cell_b )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: set_adjacent_cells" ) ;
   PEL_CHECK_PRE( cell_a != 0 && cell_b != 0 ) ;
   PEL_CHECK_PRE( ( cell_a->refinement_level() == refinement_level() &&
                    cell_b->refinement_level() <= refinement_level() )
                  ||
                  ( cell_b->refinement_level() == refinement_level() &&
                    cell_a->refinement_level() <= refinement_level() ) ) ;
                  
   CELL_1 = cell_a ;
   CELL_2 = cell_b ;
}

//----------------------------------------------------------------------
void
PDE_FaceFE:: insert_adjacent_bound( PDE_BoundFE* a_bound )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: insert_adjacent_bound" ) ;
   PEL_CHECK_PRE( side_id() == PEL::bad_index() ) ;
   PEL_CHECK_PRE( a_bound->refinement_level() == refinement_level() ) ;
   
   //???? 
   PEL_ASSERT( a_bound->refinement_level() == refinement_level() ) ;

   bound = a_bound ;
}

//----------------------------------------------------------------------
void 
PDE_FaceFE:: set_side_id( size_t an_id )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: set_side_id" ) ;
   
   SIDE_ID = an_id ;

   PEL_CHECK_POST( side_id() == an_id ) ;
}

//----------------------------------------------------------------------
size_t
PDE_FaceFE:: side_id( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: side_id" ) ;
   
   size_t result = SIDE_ID ;

   PEL_CHECK_POST( IMPLIES( ( side_id() == PEL::bad_index() ), 
                            ( has_adjacent_bound() || is_periodic() ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_FaceFE:: has_adjacent_bound( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: has_adjacent_bound" ) ;

   bool result = ( bound != 0 ) ;

   PEL_CHECK_POST( IMPLIES( result, side_id()==PEL::bad_index() ) ) ;
   return result ;
}

//----------------------------------------------------------------------
PDE_BoundFE*
PDE_FaceFE:: adjacent_bound( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: adjacent_bound" ) ;
   PEL_CHECK_PRE( has_adjacent_bound() ) ;

   PDE_BoundFE* result = bound ;

   PEL_CHECK_POST( result->refinement_level() == refinement_level() ) ;
   PEL_CHECK_POST( result->refinement_level() ==
                   adjacent_cell( 0 )->refinement_level() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_FaceFE:: nb_adjacent_cells( void ) const
//----------------------------------------------------------------------
{
   size_t result = 0 ;
   if( CELL_2 != 0 )
   {
      result = 2 ;
   }
   else if( CELL_1 != 0 )
   {
      result = 1 ;
   }

   PEL_CHECK_POST( result==1 || result==2 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_FaceFE:: has_adjacent_cell( PDE_CellFE const* m ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: has_adjacent_cell" ) ;
   PEL_CHECK_PRE( m!= 0 ) ;
   
   bool result = false ;
   if( CELL_1 != 0 && CELL_1 == m )
   {
      PEL_ASSERT( CELL_1->id_number() == m->id_number() ) ;
      result = true ;
   }
   else if( CELL_2 != 0 && CELL_2 == m )
   {
      PEL_ASSERT( CELL_2->id_number() == m->id_number() ) ;
      result = true ;
   }
   else
   {
      if( CELL_1 != 0 ) PEL_ASSERT( CELL_1->id_number() != m->id_number() ) ;
      if( CELL_2 != 0 ) PEL_ASSERT( CELL_2->id_number() != m->id_number() ) ;
   }
   
   //???? postcondition
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_CellFE*
PDE_FaceFE:: adjacent_cell( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: adjacent_cell" ) ;
   PEL_CHECK_PRE( i < nb_adjacent_cells() ) ;

   PDE_CellFE* result = ( i==0 ? CELL_1 : CELL_2 ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->refinement_level() <= refinement_level() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_CellFE*
PDE_FaceFE:: adjacent_cell_other_than( PDE_CellFE const* m ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: adjacent_cell_other_than" ) ;
   PEL_CHECK_PRE( !has_adjacent_bound() ) ;
   PEL_CHECK_PRE( !is_periodic() ) ;
   PEL_CHECK_PRE( adjacent_cell( 0 ) == m || adjacent_cell( 1 ) == m ) ;

   PDE_CellFE* result = 0 ;

   if( CELL_1 == m )
   {
      result = CELL_2 ;
   }
   else
   {
      PEL_CHECK( CELL_2 == m ) ;
      result = CELL_1 ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->id_number() != m->id_number() ) ;
   PEL_CHECK_POST( result->refinement_level() <= refinement_level() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_FaceFE:: set_periodicity( PDE_FaceFE* side, GE_Transform const* tr )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: set_periodicity" ) ;
   PEL_CHECK_PRE( side != 0 ) ;
   PEL_CHECK_PRE( side != this ) ;
   PEL_CHECK_PRE( tr != 0 ) ;
   
   PERIODIC_NEIGHBOUR = side ;
   PERIODIC_TRANSFORM = tr ;
   
   PEL_CHECK_POST( is_periodic() ) ;
   PEL_CHECK_POST( periodic_neighbour() == side ) ;
   PEL_CHECK_POST( periodic_transform() == tr ) ;
}

//----------------------------------------------------------------------
bool
PDE_FaceFE:: is_periodic( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: is_periodic" ) ;

   bool result = ( PERIODIC_NEIGHBOUR != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_FaceFE*
PDE_FaceFE:: periodic_neighbour( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: periodic_neighbour" ) ;
   PEL_CHECK_PRE( is_periodic() ) ;

   PDE_FaceFE* result = PERIODIC_NEIGHBOUR ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->refinement_level() == refinement_level() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Transform const*
PDE_FaceFE:: periodic_transform( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: periodic_transform" ) ;
   PEL_CHECK_PRE( is_periodic() ) ;

   GE_Transform const* result = PERIODIC_TRANSFORM ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_FaceFE:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: print" ) ;

   string space( indent_width, ' ' ) ;

   os << space << "side_id : " << side_id() << endl ;

   PDE_MeshFE::print( os, indent_width ) ;

   os << space << "adjacent cells :" ;
   if( nb_adjacent_cells() == 1 )
   {
      os << "  " << adjacent_cell( 0 )->id_number() << endl ;
   }
   else
   {
      size_t i0 = adjacent_cell( 0 )->id_number() ;
      size_t i1 = adjacent_cell( 1 )->id_number() ;
      if( i0 < i1 )
         os << "  " << i0 << "  " << i1 << endl ;
      else
         os << "  " << i1 << "  " << i0 << endl ;
   }
   if( nb_childs() != 0 )
   {
      os << space << "childs :" ;
      for( size_t i=0 ; i<nb_childs() ; ++i )
      {
         os << "  " << child( i )->side_id() 
            << "(" << child( i )->id_number() << ")" ; 
      }
      os << endl ;
   }
   else
   {
      os << space << "no child" << endl ;
   }
   if( has_adjacent_bound() )
   {
      os << space << "adjacent bound :  " 
         << adjacent_bound()->id_number() << endl ;   
   }
   else
   {
      os << space << "no adjacent bound" << endl ;   
   }
   if( is_periodic() )
   {
      PDE_FaceFE const* ps = periodic_neighbour() ;
      PEL_ASSERT( ps->nb_adjacent_cells() == 1 ) ;
      os << space << "periodic neighbour of color \""
                  << ps->color()->name() << "\" adjacent to cell " ;
      os << ps->adjacent_cell(0)->id_number() << endl ;
   }
}

//-----------------------------------------------------------------------
void
PDE_FaceFE:: set_parent( PDE_FaceFE* a_parent )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: set_parent" ) ; 
   PEL_CHECK_PRE( a_parent->color() == color() ) ;
   
   PARENT = a_parent ;
}

//----------------------------------------------------------------------
PDE_FaceFE*
PDE_FaceFE:: parent( void ) const
//----------------------------------------------------------------------
{
   return( PARENT ) ;
}

//----------------------------------------------------------------------
void
PDE_FaceFE:: set_nb_childs( size_t a_nb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: set_nb_childs" ) ;
   
   if( CHILDS == 0 )
   {
      CHILDS = PEL_Vector::create( this, a_nb ) ;      
   }
   else
   {
      PEL_ASSERT( CHILDS->index_limit() == a_nb ) ;
   }
}

//----------------------------------------------------------------------
size_t
PDE_FaceFE:: nb_childs( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: set_nb_childs" ) ;

   size_t result = ( CHILDS==0 ? 0 : CHILDS->index_limit() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_FaceFE:: append_child( PDE_FaceFE* a_child )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: append_child" ) ;
   PEL_CHECK_PRE( nb_childs() != 0 ) ;
   PEL_CHECK_PRE( a_child != 0 ) ;
   PEL_CHECK_PRE( a_child->color() == color() ) ;
   PEL_CHECK_PRE( a_child->refinement_level() == refinement_level()+1 ) ;
   
   size_t i=0 ;
   for( ; i<CHILDS->index_limit() ; ++i )
   {
      if( CHILDS->at( i ) == 0 )
         break ; 
      else
      {
         PDE_FaceFE* face = static_cast< PDE_FaceFE* >( CHILDS->at( i ) ) ;
         PEL_ASSERT( face != a_child ) ;
         PEL_ASSERT( face->id_number() != a_child->id_number() ) ;
      }
   }
   CHILDS->set_at( i, a_child ) ;
   
   PEL_CHECK_POST( child( i ) == a_child ) ;
}

//----------------------------------------------------------------------
PDE_FaceFE*
PDE_FaceFE:: child( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: child" ) ;
   PEL_CHECK_PRE( i < nb_childs() ) ;

   PDE_FaceFE* result = static_cast< PDE_FaceFE* >( CHILDS->at( i ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->refinement_level() == refinement_level()+1 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_FaceFE:: is_active( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FaceFE:: is_active" ) ;

   bool result = false ;
   if( CELL_1 == 0 )
   {
      PEL_ASSERT( CELL_2 == 0 ) ;
      result = false ;
   }
   else
   {
      result = CELL_1->is_active() ;
      if( result )
      {
         if( CELL_2 != 0 )
         {
            PEL_ASSERT( PERIODIC_NEIGHBOUR == 0 ) ;
            result = CELL_2->is_active() ;
         }
         else if( PERIODIC_NEIGHBOUR != 0 )
         {
            PEL_ASSERT( PERIODIC_NEIGHBOUR->CELL_2 == 0 ) ;
            if( PERIODIC_NEIGHBOUR->CELL_1 == 0 )
            {
               result = false ;
            }
            else
            {
               result = PERIODIC_NEIGHBOUR->CELL_1->is_active() ;
            }
         }
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DiscOnMeshFE const*
PDE_FaceFE:: disc( void ) const
//----------------------------------------------------------------------
{
   return( 0 ) ;
}

