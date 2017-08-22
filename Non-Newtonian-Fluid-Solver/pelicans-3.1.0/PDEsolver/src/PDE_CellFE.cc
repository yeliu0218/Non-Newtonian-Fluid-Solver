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

#include <PDE_CellFE.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PDE_DiscOnMeshFE.hh>
#include <PDE_FacesOfCellFE.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_FaceFE.hh>

#include <PEL.hh>
#include <PEL_List.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <iostream>

using std::cout ;
using std::string ;
using std::endl ;

//----------------------------------------------------------------------
PDE_CellFE*
PDE_CellFE:: create( PEL_Object* a_owner,
                     size_t a_number,
                     GE_Mpolyhedron* a_polyhedron,
                     GE_Color const* a_color,
                     size_t a_refinement_level,
                     PDE_DiscOnMeshFE const* a_disc )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: create" ) ;

   PEL_CHECK_PRE( a_polyhedron!=0 ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;
   
   PDE_CellFE* result = new PDE_CellFE( a_owner, 
                                        a_number, 
                                        a_polyhedron,
                                        a_color,
                                        a_refinement_level,
                                        a_disc ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->id_number()==a_number ) ;
   PEL_CHECK_POST( result->polyhedron()==a_polyhedron ) ;
   PEL_CHECK_POST( result->color()==a_color ) ;
   PEL_CHECK_POST( result->nb_faces()==0 ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_CellFE:: PDE_CellFE( PEL_Object* a_owner,
                         size_t a_number,
                         GE_Mpolyhedron* a_polyhedron,
                         GE_Color const* a_color,
                         size_t a_refinement_level,
                         PDE_DiscOnMeshFE const* a_disc )
//----------------------------------------------------------------------
   : PDE_MeshFE( a_owner, a_number, a_polyhedron, a_color, a_refinement_level )
   , DISC( a_disc )
   , FACES( PEL_Vector::create( this, 0 ) )
   , PARENT( 0 )
   , CHILDS( 0 )
   , ACTIVE( true )
{
   if( a_disc != 0 )
   {
      for( size_t e=0 ; e<a_disc->nb_reference_elements() ; ++e )
      {
         PDE_ReferenceElement const* elm = a_disc->reference_element( e );
         std::vector< PDE_BasisFunctionCell* > bf( elm->nb_nodes() ,
                                              (PDE_BasisFunctionCell*) 0 ) ;
         BFS.push_back( bf ) ;
      }
   }
}

//----------------------------------------------------------------------
PDE_CellFE*
PDE_CellFE:: create_with_discretizations_pattern( 
                     PEL_Object* a_owner,
                     size_t a_number,
                     GE_Mpolyhedron* a_polyhedron,
                     GE_Color const* a_color,
                     size_t a_refinement_level,
                     PDE_CellFE const* a_pattern_mesh )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: create_with_discretization_pattern" ) ;

   PEL_CHECK_PRE( a_polyhedron!=0 ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;
   
   PDE_CellFE* result = new PDE_CellFE( a_owner, 
                                        a_number, 
                                        a_polyhedron,
                                        a_color,
                                        a_refinement_level,
                                        a_pattern_mesh ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->id_number()==a_number ) ;
   PEL_CHECK_POST( result->polyhedron()==a_polyhedron ) ;
   PEL_CHECK_POST( result->color()==a_color ) ;
   PEL_CHECK_POST( result->nb_faces()==0 ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_CellFE:: PDE_CellFE( PEL_Object* a_owner,
                         size_t a_number,
                         GE_Mpolyhedron* a_polyhedron,
                         GE_Color const* a_color,
                         size_t a_refinement_level,
                         PDE_CellFE const* a_pattern_mesh )
//----------------------------------------------------------------------
   : PDE_MeshFE( a_owner, a_number, a_polyhedron, a_color, a_refinement_level,
                 a_pattern_mesh )
   , DISC( a_pattern_mesh->DISC )
   , FACES( PEL_Vector::create( this, 0 ) )
   , PARENT( 0 )
   , CHILDS( 0 )
   , ACTIVE( true )
{
   for( size_t i=0 ; i<a_pattern_mesh->nb_reference_elements() ; ++i )
   {
      PDE_ReferenceElement const* elm = a_pattern_mesh->reference_element( i );
      std::vector< PDE_BasisFunctionCell* > bd( elm->nb_nodes() ,
                                            (PDE_BasisFunctionCell*) 0 ) ;
      BFS.push_back( bd ) ;
   }
}

//----------------------------------------------------------------------
PDE_CellFE:: ~PDE_CellFE( void )
//----------------------------------------------------------------------
{
}

//------------------------------------------------------------------------
bool
PDE_CellFE:: all_basis_functions_are_dropped( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: all_basis_functions_are_dropped" ) ;
   
   bool result = true ;
   for( size_t i=0 ; result==true && i<nb_childs() ; ++i )
   {
      result = child( i )->all_basis_functions_are_dropped() ;
   }
   for( size_t ee=0 ; result==true && ee<nb_reference_elements() ; ++ee )
   {
      size_t nb_lns = reference_element( ee )->nb_nodes() ;
      for( size_t ln=0 ; result==true && ln<nb_lns ; ++ln )
      {
         PDE_BasisFunction const* bf = basis_function( ee, ln ) ;
         result = ( bf == 0 ) || 
                  ( !bf->is_attached_to_valid_DOFs() ) ;
      }
   }
   return( result ) ;
}

//------------------------------------------------------------------------
void
PDE_CellFE:: set_basis_function( size_t ee, size_t ln, 
                                 PDE_BasisFunctionCell* bf )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: set_basis_function" ) ;
   PEL_CHECK_PRE( ee < nb_reference_elements() ) ;
   PEL_CHECK_PRE( ln < nb_basis_functions(ee) ) ;
   PEL_CHECK_PRE( bf != 0 ) ;
   PEL_CHECK_PRE( ( bf->refinement_level() == PEL::bad_index() ) ||
                  ( bf->refinement_level() == refinement_level() ) ) ;
   PEL_CHECK_PRE( basis_function( ee, ln ) == 0 ) ;

   BFS[ee][ln] = bf ;

   PEL_CHECK_POST( basis_function( ee, ln ) == bf ) ;
}

//------------------------------------------------------------------------
void
PDE_CellFE:: remove_basis_function( size_t ee, size_t ln )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: remove_basis_function" ) ;
   PEL_CHECK_PRE( ee < nb_reference_elements() ) ;
   PEL_CHECK_PRE( ln < nb_basis_functions(ee) ) ;
   PEL_CHECK_PRE( basis_function( ee, ln ) != 0 ) ;

   BFS[ee][ln] = 0 ;

   PEL_CHECK_POST( basis_function( ee, ln ) == 0 ) ;
}

//------------------------------------------------------------------------
PDE_BasisFunctionCell*
PDE_CellFE:: basis_function( size_t ee, size_t ln ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: basis_function" ) ;
   PEL_CHECK_PRE( basis_function_PRE( ee, ln ) ) ;

   PDE_BasisFunctionCell* result = BFS[ee][ln] ;

   PEL_CHECK_POST( basis_function_POST( result, ee, ln ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CellFE:: append_face( PDE_FaceFE* a_face )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: append_face" ) ;
   PEL_CHECK_PRE( a_face->refinement_level() == refinement_level() ) ;

   //??????? on connait à l'avance la taille ?????????????????
   FACES->append( a_face ) ; 
}

//----------------------------------------------------------------------
size_t
PDE_CellFE:: nb_faces( void ) const
//----------------------------------------------------------------------
{
   return( FACES->count() ) ;
}

//----------------------------------------------------------------------
PEL_Vector const*
PDE_CellFE:: faces( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: side" ) ;

   PEL_Vector const* result = FACES ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->count() == nb_faces() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_FacesOfCellFE*
PDE_CellFE:: create_faces_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: create_faces_iterator" ) ;
   
   PDE_FacesOfCellFE* result = new PDE_FacesOfCellFE( a_owner, FACES ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_CellFE:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: print" ) ;

   string space( indent_width, ' ' ) ;
   
   PDE_MeshFE::print( os, indent_width ) ;

   os << space << "faces :" ;
   PEL_VectorIterator* it = PEL_VectorIterator::create( 0, FACES ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_FaceFE const* face = static_cast<PDE_FaceFE*>( it->item() ) ;
      os << "  " ;
      if( face->side_id() != PEL::bad_index() ) os << face->side_id() ;
      os << "(" << face->id_number() << ")" ;
   }
   os << endl ;
   if( nb_childs() != 0 )
   {
      os << space << "childs :" ;
      for( size_t i=0 ; i<nb_childs() ; ++i )
      {
         os << "  " << child( i )->id_number() ;
      }
      os << endl ;
   }
   else
   {
      os << space << "no child" << endl ;
   }
   it->destroy() ;
}

//----------------------------------------------------------------------
void
PDE_CellFE:: set_parent( PDE_CellFE* a_parent )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: set_parent" ) ;
   PEL_CHECK_PRE( a_parent != 0 ) ;
   PEL_CHECK_PRE( a_parent->refinement_level() == (refinement_level()-1) ) ;
   PEL_CHECK_PRE( a_parent->color() == color() ) ;
   
   PARENT = a_parent ;
}

//----------------------------------------------------------------------
PDE_CellFE*
PDE_CellFE:: parent( void ) const
//----------------------------------------------------------------------
{
   PDE_CellFE* result = PARENT ;

   PEL_ASSERT( 
      IMPLIES( result != 0, 
               result->refinement_level() == (refinement_level()-1) ) ) ;
   return( PARENT ) ;
}

//----------------------------------------------------------------------
void
PDE_CellFE:: set_nb_childs( size_t a_nb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: set_nb_childs" ) ;

   PEL_ASSERT( CHILDS == 0 ) ;
   CHILDS = PEL_Vector::create( this, a_nb ) ;
}

//----------------------------------------------------------------------
size_t
PDE_CellFE:: nb_childs( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: set_nb_childs" ) ;

   size_t result = ( CHILDS==0 ? 0 : CHILDS->index_limit() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_CellFE:: set_child( size_t i, PDE_CellFE* a_child )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: set_child" ) ;
   PEL_CHECK_PRE( nb_childs() != 0 ) ;
   PEL_CHECK_PRE( i < nb_childs() ) ;
   PEL_CHECK_PRE( a_child != 0 ) ;
   PEL_CHECK_PRE( a_child->color() == color() ) ;
   PEL_CHECK_PRE( a_child->refinement_level() == refinement_level()+1 ) ;

   CHILDS->set_at( i, a_child ) ;

   PEL_CHECK_POST( child( i ) == a_child ) ;
}

//-----------------------------------------------------------------------
PDE_CellFE*
PDE_CellFE:: child( size_t i ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: child" ) ;
   PEL_CHECK_PRE( i < nb_childs() ) ;

   PDE_CellFE* result = static_cast< PDE_CellFE* >( CHILDS->at( i ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_ASSERT( result->refinement_level() == refinement_level()+1 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_CellFE:: is_active( void ) const
//----------------------------------------------------------------------
{
   return( ACTIVE ) ;
}

//----------------------------------------------------------------------
void
PDE_CellFE:: set_active( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: set_active" ) ;

   ACTIVE = true ;

   PEL_CHECK_POST( is_active() ) ;
}

//----------------------------------------------------------------------
void
PDE_CellFE:: set_inactive( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CellFE:: set_inactive" ) ;

   ACTIVE = false ;

   PEL_CHECK_POST( !is_active() ) ;
}

//----------------------------------------------------------------------
PDE_DiscOnMeshFE const*
PDE_CellFE:: disc( void ) const
//----------------------------------------------------------------------
{
   return( DISC ) ;
}
