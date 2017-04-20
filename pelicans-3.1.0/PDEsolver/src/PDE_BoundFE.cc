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

#include <PDE_BoundFE.hh>

#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DiscOnMeshFE.hh>
#include <PDE_ReferenceElement.hh>

#include <GE_QuadratureRule.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Vector.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;

//----------------------------------------------------------------------
PDE_BoundFE*
PDE_BoundFE:: create( PEL_Object* a_owner,
                      size_t a_number,
                      GE_Mpolyhedron* a_polyhedron,
                      GE_Color const* a_color,
                      size_t a_refinement_level,
                      PDE_DiscOnMeshFE const* a_disc ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BoundFE:: create" ) ;
   PEL_CHECK_PRE( a_polyhedron!=0 ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;
   
   PDE_BoundFE* result = new PDE_BoundFE( a_owner, 
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
   PEL_CHECK_POST( !result->has_adjacent_cell() ) ;
   PEL_CHECK_POST( !result->has_adjacent_bound() ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_BoundFE:: PDE_BoundFE( PEL_Object* a_owner,
                           size_t a_number,
                           GE_Mpolyhedron* a_polyhedron,
                           GE_Color const* a_color,
                           size_t a_refinement_level,
                           PDE_DiscOnMeshFE const* a_disc ) 
//----------------------------------------------------------------------
   : PDE_MeshFE( a_owner, a_number,a_polyhedron, a_color, a_refinement_level )
   , DISC( a_disc )
   , adjacentMesh( 0 )
   , ADJ_BOUND( 0 )
   , outwardNormal( GE_Vector::create( this,
                                       polyhedron()->nb_space_dimensions() ) )
   , PARENT( 0 )
{
   if( a_disc != 0 )
   {
      for( size_t e=0 ; e<a_disc->nb_reference_elements() ; ++e )
      {
         PDE_ReferenceElement const* elm = a_disc->reference_element( e );
         std::vector< PDE_BasisFunctionBound* > bf( elm->nb_nodes() ,
                                                (PDE_BasisFunctionBound*) 0 ) ;
         BFS.push_back( bf ) ;
      }
   }
}

//----------------------------------------------------------------------
PDE_BoundFE:: ~PDE_BoundFE( void )
//----------------------------------------------------------------------
{
}

//------------------------------------------------------------------------
void
PDE_BoundFE:: set_basis_function( size_t ee, size_t ln, 
                                  PDE_BasisFunctionBound* bf )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: set_basis_function" ) ;
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
PDE_BasisFunctionBound*
PDE_BoundFE:: basis_function( size_t ee, size_t ln ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: basis_function" ) ;
   PEL_CHECK_PRE( basis_function_PRE( ee, ln ) ) ;

   PDE_BasisFunctionBound* result = BFS[ee][ln] ;

   PEL_CHECK_POST( basis_function_POST( result, ee, ln ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_BoundFE:: insert_adjacent_cell( PDE_CellFE* mesh )
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( !has_adjacent_cell() ) ;

   adjacentMesh = mesh ;

   PEL_CHECK_POST( has_adjacent_cell() ) ;
}

//----------------------------------------------------------------------
PDE_CellFE*
PDE_BoundFE:: adjacent_cell( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( has_adjacent_cell() ) ;

   return( adjacentMesh ) ;
}

//----------------------------------------------------------------------
GE_Vector const*
PDE_BoundFE:: unit_outward_normal( void ) const
//----------------------------------------------------------------------
{
   // Set the outward normal :
   {
      outwardNormal->set( polyhedron()->unit_normal() ) ;
      GE_Vector* uu = GE_Vector::create( 0,
                                         adjacentMesh->polyhedron()->center(),
                                         polyhedron()->vertex(0) ) ;
      if( outwardNormal->dot_product(uu)>0.0 )
      {
         outwardNormal->scale( -1.0 ) ;
      }
      uu->destroy() ;      
   }
   
   GE_Vector const* result = outwardNormal ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->nb_components()==
                                polyhedron()->nb_space_dimensions() ) ;
   PEL_CHECK_POST( PEL::equal( result->norm(), 1.0 ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_BoundFE:: value( PDE_DiscreteField const* f,
                     size_t level,
                     GE_Point const* ptRefCoord,
                     size_t iComp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BoundFE:: value" ) ;
   PEL_CHECK_PRE( value_PRE( f, level, ptRefCoord, iComp ) ) ;

   double val ;
//   PDE_LocalField const* lf = getLocalField( f ) ;
   if( has_discretization( f ) )
   {
//      PEL_CHECK( adjacentMesh->getLocalField( f ) == 0 ) ;
      PEL_CHECK( !adjacentMesh->has_discretization( f ) ) ;
//??????? à vérifier
      val = PDE_MeshFE::value( f, level, ptRefCoord, iComp ) ;
//???      val = lf->value( level, ptRefCoord, iComp ) ;
   }
   else
   {
      GE_Point* ptRealCoord = GE_Point::create( 0, polyhedron()->nb_space_dimensions() ) ;
      polyhedron()->apply_mapping( ptRefCoord, ptRealCoord ) ;
      GE_Point* meshRefCoord = GE_Point::create( 0, polyhedron()->nb_space_dimensions() ) ;
      adjacentMesh->polyhedron()->apply_inverse_mapping( ptRealCoord, meshRefCoord ) ;
      val = adjacentMesh->value( f, level, meshRefCoord, iComp ) ;
      ptRealCoord->destroy() ;
      meshRefCoord->destroy() ;
   }

   return( val ) ;
}

//----------------------------------------------------------------------
bool
PDE_BoundFE:: isInFlow( size_t level, PDE_DiscreteField const* velocity ) const
//----------------------------------------------------------------------
{
   PEL_CHECK( velocity->nb_components() == polyhedron()->nb_space_dimensions() ) ;

   bool resu = false ;

   GE_Point* vertexRefCoord = GE_Point::create( 0, polyhedron()->nb_space_dimensions()-1 ) ;

   for( size_t iVert=0 ; iVert<polyhedron()->nb_vertices() ; iVert++ )
   {
      polyhedron()->apply_inverse_mapping( polyhedron()->vertex(iVert), vertexRefCoord ) ;
      if( isInFlow( level, velocity, vertexRefCoord ) )
      {
         resu = true ;
         break ;
      }
   }

   vertexRefCoord->destroy() ;
   return( resu ) ;
}

//----------------------------------------------------------------------
bool
PDE_BoundFE:: isInFlow( size_t level,
                        PDE_DiscreteField const* velocity,
                        GE_Point const* ptRefCoord ) const
//----------------------------------------------------------------------
{
   int nbSpDims = polyhedron()->nb_space_dimensions() ;

   PEL_CHECK( velocity->nb_components()== polyhedron()->nb_space_dimensions() ) ;

   GE_Vector* velo = GE_Vector::create( 0, nbSpDims ) ;
   for( int iDim=0 ; iDim<nbSpDims ; iDim++ )
   {
      velo->set_component( iDim,
                     value( velocity, level, ptRefCoord, iDim ) ) ;
   }

   bool resu = ( velo->dot_product( unit_outward_normal() ) < -1.E-8 ) ;

   velo->destroy() ;
   return( resu ) ;
}

//----------------------------------------------------------------------
bool
PDE_BoundFE:: has_adjacent_cell( void ) const
//----------------------------------------------------------------------
{
   bool result = ( adjacentMesh != 0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PDE_BoundFE:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEcell:: print" ) ;

   string space( indent_width, ' ' ) ;

   PDE_MeshFE::print( os, indent_width ) ;

   if( has_adjacent_cell() )
   {
      os << space << "adjacent cell :  " 
                  <<  adjacent_cell()->id_number() << endl ;
   }
   else
   {
      os << space << "no adjacent cell" << endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_BoundFE:: set_parent( PDE_BoundFE* a_parent )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BoundFE:: set_parent" ) ;
   PEL_CHECK_PRE( a_parent->refinement_level() == (refinement_level()-1) ) ;
   PEL_CHECK_PRE( a_parent->color() == color() ) ;
   
   PARENT = a_parent ;
}

//----------------------------------------------------------------------
PDE_BoundFE*
PDE_BoundFE:: parent( void ) const
//----------------------------------------------------------------------
{
   return( PARENT ) ;
}

//----------------------------------------------------------------------
bool
PDE_BoundFE:: is_active( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BoundFE:: is_active" ) ;

   bool result = false ;
   if( adjacentMesh != 0 ) result = adjacentMesh->is_active() ;

   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_BoundFE:: has_adjacent_bound( void ) const
//----------------------------------------------------------------------
{
   bool result = ( ADJ_BOUND != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_BoundFE const*
PDE_BoundFE:: adjacent_bound( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BoundFE:: adjacent_bound" ) ;
   PEL_CHECK_PRE( has_adjacent_bound() ) ;

   return( ADJ_BOUND ) ;
}

//----------------------------------------------------------------------
void
PDE_BoundFE:: insert_adjacent_bound( PDE_BoundFE const* a_bound )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BoundFE:: insert_adjacent_bound" ) ;
   PEL_CHECK_PRE( !has_adjacent_bound() ) ;

   ADJ_BOUND = a_bound ;

   PEL_CHECK_POST( has_adjacent_bound() ) ;
   PEL_CHECK_POST( adjacent_bound() == a_bound ) ;
}

//----------------------------------------------------------------------
PDE_DiscOnMeshFE const*
PDE_BoundFE:: disc( void ) const
//----------------------------------------------------------------------
{
   return( DISC ) ;
}
