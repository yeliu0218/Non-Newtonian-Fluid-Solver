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

#include <PDE_MeshFE.hh>

#include <PEL_Vector.hh>

#include <PEL_Error.hh>
#include <PEL_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PDE_BFvalues.hh>
#include <PDE_BoundFE.hh>
#include <PDE_DiscreteField.hh>
//#include <PDE_DiscOnMeshFE.hh>
#include <PDE_DiscOnMeshFE.hh>
#include <PDE_ReferenceElement.hh>

#include <algorithm>
#include <iostream>
#include <set>

using std::cout ; using std::endl ;
using std::set ;
using std::string ;
using std::vector ;

PDE_RefinementPatternProvider const* PDE_MeshFE::RPP = 0 ;

//----------------------------------------------------------------------
PDE_MeshFE:: PDE_MeshFE( PEL_Object* a_owner,
                         size_t a_number,
                         GE_Mpolyhedron* a_polyhedron,
                         GE_Color const* a_color,
                         size_t a_refinement_level )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , ID( a_number )
   , POLY( a_polyhedron )
   , COLOR( a_color )
   , RLEVEL( a_refinement_level )
{
}

//----------------------------------------------------------------------
PDE_MeshFE:: PDE_MeshFE( PEL_Object* a_owner,
                         size_t a_number,
                         GE_Mpolyhedron* a_polyhedron,
                         GE_Color const* a_color,
                         size_t a_refinement_level,
                         PDE_MeshFE const* a_pattern_mesh )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , ID( a_number )
   , POLY( a_polyhedron )
   , COLOR( a_color )
   , RLEVEL( a_refinement_level )
{
}

//----------------------------------------------------------------------
PDE_MeshFE:: ~PDE_MeshFE( void )
//----------------------------------------------------------------------
{
}

//------------------------------------------------------------------------
bool
PDE_MeshFE:: is_equal( PEL_Object const* other ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;

   bool result = three_way_comparison( other )==0 ;

   PEL_CHECK_POST( is_equal_POST( result, other ) ) ;   
   return( result ) ;
}

//------------------------------------------------------------------------
int
PDE_MeshFE:: three_way_comparison( PEL_Object const* other ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   // Mesh comparison with another mesh :
   PDE_MeshFE const* m = static_cast<PDE_MeshFE const*>( other ) ;
   
   int result = polyhedron()->nb_vertices() - 
      m->polyhedron()->nb_vertices() ;
   if( result == 0 )
   {
      result = polyhedron()->center()->three_way_comparison(
         m->polyhedron()->center() ) ;
         
   }   
   PEL_CHECK( ( ( result == 0) && 
                ( polyhedron()->is_equal( m->polyhedron() ) ) ) || 
              ( ( result != 0) && 
                !( polyhedron()->is_equal( m->polyhedron() ) )
                 ) ) ;
                  
   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t
PDE_MeshFE:: hash_code( void ) const
//------------------------------------------------------------------------
{
   return polyhedron()->center()->hash_code() ;
}

//----------------------------------------------------------------------
size_t
PDE_MeshFE:: id_number( void ) const
//----------------------------------------------------------------------
{
   return( ID ) ;
}

//----------------------------------------------------------------------
GE_Mpolyhedron*
PDE_MeshFE:: polyhedron( void ) const
//----------------------------------------------------------------------
{
   return( POLY ) ;
}

//----------------------------------------------------------------------
GE_Color const*
PDE_MeshFE:: color( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK( COLOR != 0 ) ;   

   return( COLOR ) ;
}

//----------------------------------------------------------------------
void
PDE_MeshFE:: duplicate_discretization( PDE_DiscreteField const* model_f,
                                       PDE_DiscreteField* new_f )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: duplicate_discretization" ) ;
   PEL_CHECK_PRE( model_f != 0 ) ;
   PEL_CHECK_PRE( new_f != 0 ) ;
   PEL_CHECK_PRE( has_discretization( model_f ) ) ;

   size_t ee = disc()->index_of_reference_element( new_f ) ;
   PDE_ReferenceElement const* elm = disc()->reference_element( ee ) ;
   for( size_t ln=0 ; ln<elm->nb_nodes() ; ++ln )
   {
      PDE_BasisFunction* bf = basis_function( ee, ln ) ;
      if( bf != 0 )
      {
         if( !bf->is_in_basis_function_set_of( new_f ) )
         {
            bf->append_one_field( new_f, bf->node_of_DOF( model_f ) ) ;
         }
         else
         {
            PEL_ASSERT( bf->node_of_DOF( new_f ) == 
                        bf->node_of_DOF( model_f ) ) ;
         }
      }
   }
   
   PEL_CHECK_POST( has_discretization( new_f ) ) ;
   PEL_CHECK_POST( index_of_reference_element( new_f ) ==
                              index_of_reference_element( model_f ) ) ;
}

//??? not correct during the building stage
//----------------------------------------------------------------------
bool
PDE_MeshFE:: has_discretization( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: has_discretization" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   bool result = ( disc() != 0 ) ;
   if( result )
   {
      result = disc()->has_discretization( ff ) ;
   }

   return( result ) ;
}

//------------------------------------------------------------------------
size_t
PDE_MeshFE:: index_of_reference_element( PDE_DiscreteField const* ff ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: index_of_reference_element" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   
   size_t result = ( disc() == 0 ? 
                        PEL::bad_index() :
                     disc()->index_of_reference_element( ff ) ) ;
 
   PEL_CHECK_POST( EQUIVALENT( has_discretization( ff ),
                               result != PEL::bad_index() ) ) ;
   PEL_CHECK_POST( EQUIVALENT( has_discretization( ff ),
                               result < nb_reference_elements() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double 
PDE_MeshFE:: value( PDE_DiscreteField const* ff,
		    size_t level,
		    GE_Point const* pt_ref,
		    size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: value" ) ;
   PEL_CHECK_PRE( value_PRE( ff, level, pt_ref, ic ) ) ;
   PEL_CHECK_PRE( has_discretization( ff ) ) ;
   PEL_CHECK_PRE(
      FORALL(
         ( size_t i=0 ;
           i<reference_element( index_of_reference_element( ff ) )->nb_nodes() ;
           ++i ),
         basis_function( index_of_reference_element( ff ), i ) != 0 ) ) ;

   double result = 0.0 ;

   size_t ee = disc()->index_of_reference_element( ff ) ;
   PDE_ReferenceElement const* elm = disc()->reference_element( ee ) ;

   for( size_t ln=0 ; ln<elm->nb_nodes() ; ln++ )
   {
      PDE_BasisFunction* bf = basis_function( ee, ln ) ;
      size_t n = bf->node_of_DOF( ff ) ;
      result += ff->DOF_value( level, n, ic )
                                     * elm->N_local( ln, pt_ref ) ;

//       result += ff->DOF_value( level, lf.global_nodes(i), ic ) 
//               * elm->N_local( i, pt_ref ) ;

   }
   //????????? faux ????? prendre en compte le raffinement
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_MeshFE:: nb_discrete_fields( size_t ee ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: nb_discrete_fields" ) ;
   PEL_CHECK_PRE( ee < nb_reference_elements() ) ;
   
   size_t result = disc()->nb_discrete_fields( ee ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_MeshFE:: discrete_field( size_t ee, size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: discrete_field" ) ;
   PEL_CHECK_PRE( ee<nb_reference_elements() ) ;
   PEL_CHECK_PRE( i<nb_discrete_fields(ee) ) ;
   
   PDE_DiscreteField const* result = disc()->discrete_field( ee, i ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_MeshFE:: index_of_element( PDE_ReferenceElement const* elm ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: index_of_element" ) ;

   size_t result = disc()->index_of_element( elm ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_MeshFE:: nb_reference_elements( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: nb_reference_elements" ) ;

   size_t result = ( disc()==0 ? 0 : disc()->nb_reference_elements() ) ;

   return( result ) ;
}

//-----------------------------------------------------------------------
PDE_ReferenceElement const*
PDE_MeshFE:: reference_element( size_t ee ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: reference_element" ) ;
   PEL_CHECK_PRE( ee < nb_reference_elements() ) ;

   PDE_ReferenceElement const* result = disc()->reference_element( ee ) ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_MeshFE:: nb_basis_functions( size_t ee ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: nb_basis_functions" ) ;
   PEL_CHECK_PRE( ee < nb_reference_elements() ) ;

   size_t result = disc()->nb_basis_functions( ee ) ;

   PEL_CHECK_POST( result == reference_element( ee )->nb_nodes() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_MeshFE:: set_refinement_pattern_provider( 
                            PDE_RefinementPatternProvider const* r )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: set_refinement_pattern_provider" ) ;
   PEL_CHECK_PRE( refinement_pattern_provider() == 0 || 
                  refinement_pattern_provider() == r ) ;

   RPP = r ;

   PEL_CHECK_POST( refinement_pattern_provider() == r ) ;
}

//----------------------------------------------------------------------
PDE_RefinementPatternProvider const*
PDE_MeshFE:: refinement_pattern_provider( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: refinement_pattern_provider" ) ;
   PEL_CHECK_PRE( refinement_pattern_provider() != 0 ) ;
   
   return( RPP ) ;
}

//----------------------------------------------------------------------
size_t
PDE_MeshFE:: refinement_level( void ) const
//----------------------------------------------------------------------
{
   return( RLEVEL ) ;
}

//----------------------------------------------------------------------
void
PDE_MeshFE:: set_color( GE_Color const* col )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: set_color" ) ;
   PEL_CHECK_PRE( col != 0 ) ;
   
   COLOR = col ;
   
   PEL_CHECK_POST( color() == col ) ;
}

//-----------------------------------------------------------------------
void
PDE_MeshFE:: print( std::ostream& os, size_t indent_width ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: print" ) ;
   string bl( indent_width, ' ' ) ;
   os << bl << "id_number : " << id_number() 
      << "  color : \"" << color()->name() << "\"";
   if( is_active() ) 
      os << "  active" ;
   else 
      os << "  inactive" ;
   os << std::endl ;
   polyhedron()->print( os, indent_width ) ;
}

//----------------------------------------------------------------------
bool
PDE_MeshFE:: value_PRE( PDE_DiscreteField const* ff,
                        size_t level,
                        GE_Point const* pt_ref,
                        size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( level < ff->storage_depth() ) ;
   PEL_ASSERT( pt_ref != 0 ) ;
   PEL_ASSERT( polyhedron()->reference_polyhedron()->contains( pt_ref ) ) ;
   PEL_ASSERT( ic < ff->nb_components() ) ;

   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_MeshFE:: basis_function_PRE( size_t ee, size_t ln ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( ee < nb_reference_elements() ) ;
   PEL_ASSERT( ln < nb_basis_functions( ee ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_MeshFE:: basis_function_POST( PDE_BasisFunction* result,
                                  size_t ee, size_t ln ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result!=0, 
                        result->refinement_level() == refinement_level() ) ) ;
   return( true ) ;
}
                                  

