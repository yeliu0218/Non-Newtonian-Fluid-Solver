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

#include <PDE_MortarSideFE.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Vector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Color.hh>

#include <PDE_BoundFE.hh>
#include <PDE_DiscOnMeshFE.hh>
#include <PDE_ReferenceElement.hh>

#include <iostream>
using std::string ;
using std::endl ;

//----------------------------------------------------------------------------
PDE_MortarSideFE*
PDE_MortarSideFE:: create( PEL_Object* a_owner,
                           size_t a_number,
                           GE_Mpolyhedron* a_polyhedron,
                           GE_Color const* a_color,
                           size_t a_refinement_level,
                           PDE_DiscOnMeshFE const* a_disc )
//----------------------------------------------------------------------------
{
   PDE_MortarSideFE* result = new PDE_MortarSideFE( a_owner,
                                                    a_number,
                                                    a_polyhedron,
                                                    a_color,
                                                    a_refinement_level,
                                                    a_disc ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PDE_MortarSideFE:: PDE_MortarSideFE( PEL_Object* a_owner,
                                     size_t a_number,
                                     GE_Mpolyhedron* a_polyhedron,
                                     GE_Color const* a_color,
                                     size_t a_refinement_level,
                                     PDE_DiscOnMeshFE const* a_disc )
//----------------------------------------------------------------------------
   : PDE_MeshFE( a_owner, a_number, a_polyhedron, a_color, a_refinement_level )
   , DISC( a_disc )
   , BOUNDS_0( PEL_Vector::create( this, 0 ) )
   , BOUNDS_1( PEL_Vector::create( this, 1 ) )
{
   if( a_disc != 0 )
   {
      for( size_t e=0 ; e<a_disc->nb_reference_elements() ; ++e )
      {
         PDE_ReferenceElement const* elm = a_disc->reference_element( e );
         std::vector< PDE_BasisFunctionMortarSide* > bf( elm->nb_nodes() ,
                                           (PDE_BasisFunctionMortarSide*) 0 ) ;
         BFS.push_back( bf ) ;
      }
   }
}

//----------------------------------------------------------------------------
PDE_MortarSideFE:: ~PDE_MortarSideFE( void )
//----------------------------------------------------------------------------
{
}

//------------------------------------------------------------------------
void
PDE_MortarSideFE:: set_basis_function( size_t ee, size_t ln, 
                                 PDE_BasisFunctionMortarSide* bf )
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
PDE_BasisFunctionMortarSide*
PDE_MortarSideFE:: basis_function( size_t ee, size_t ln ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MeshFE:: basis_function" ) ;
   PEL_CHECK_PRE( basis_function_PRE( ee, ln ) ) ;

   PDE_BasisFunctionMortarSide* result = BFS[ee][ln] ;

   PEL_CHECK_POST( basis_function_POST( result, ee, ln ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_MortarSideFE:: append_domain_bound( size_t domain_id,
                                        PDE_BoundFE* bd )
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( domain_id==0 || domain_id==1 ) ;

   if( domain_id==0 ) 
   {
      BOUNDS_0->append( bd ) ;
   }
   else 
   {
      BOUNDS_1->append( bd ) ;
   }
}

//----------------------------------------------------------------------------
PEL_Vector const*
PDE_MortarSideFE:: domain_bounds( size_t domain_id ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( domain_id==0 || domain_id==1 ) ;

   PEL_Vector const* result = 0 ;
   if( domain_id==0 ) 
   {
      result = BOUNDS_0 ;
   }
   else 
   {
      result = BOUNDS_1 ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_MortarSideFE:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_MortarSideFE:: print" ) ;

   string space( indent_width, ' ' ) ;

   PDE_MeshFE::print( os, indent_width ) ;

   for( size_t domain_id=0 ; domain_id<2 ; ++domain_id )
   {
      os << space << "bounds of domain " << domain_id << " :" ;
      PEL_VectorIterator* it = 
               PEL_VectorIterator::create( 0, domain_bounds( domain_id ) ) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         PDE_BoundFE const* bound = static_cast<PDE_BoundFE*>( it->item() ) ;
         os << "  " << bound->id_number() ;
      }
      os << endl ;
      it->destroy() ; it=0 ;
   }
}

//----------------------------------------------------------------------
PDE_MortarSideFE*
PDE_MortarSideFE:: parent( void ) const
//----------------------------------------------------------------------
{
   return( 0 ) ;
}

//----------------------------------------------------------------------
bool
PDE_MortarSideFE:: is_active( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
PDE_DiscOnMeshFE const*
PDE_MortarSideFE:: disc( void ) const
//----------------------------------------------------------------------
{
   return( DISC ) ;
}
