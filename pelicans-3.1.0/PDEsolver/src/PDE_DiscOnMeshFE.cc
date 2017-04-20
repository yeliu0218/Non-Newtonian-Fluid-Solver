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

#include <PDE_DiscOnMeshFE.hh>

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
#include <PDE_ReferenceElement.hh>

#include <algorithm>
#include <iostream>
#include <set>

using std::cout ; using std::endl ;
using std::set ;
using std::string ;
using std::vector ;

//----------------------------------------------------------------------
PDE_DiscOnMeshFE*
PDE_DiscOnMeshFE:: create( PEL_Object* a_owner,
                           PDE_DiscreteField const* ff,
                           PDE_ReferenceElement const* elm )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: PDE_DiscOnMeshFE" ) ;
   
   PDE_DiscOnMeshFE* result = new PDE_DiscOnMeshFE( a_owner, ff, elm ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DiscOnMeshFE:: PDE_DiscOnMeshFE( PEL_Object* a_owner,
                                     PDE_DiscreteField const* ff,
                                     PDE_ReferenceElement const* elm )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
{
   add_discretization( ff, elm ) ; //??? do something better!
}

//----------------------------------------------------------------------
PDE_DiscOnMeshFE:: ~PDE_DiscOnMeshFE( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_DiscOnMeshFE:: add_discretization( PDE_DiscreteField const* ff,
                                       PDE_ReferenceElement const* elm )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: add_discretization" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( elm != 0 ) ;
      
   size_t idx = ff->id_number() ;
   if( idx >= FIELD_2_ELMS.size() )
   {
      FIELD_2_ELMS.resize( idx+1, PEL::bad_index() ) ;
   }
   size_t ee = index_of_element( elm ) ;
   if( ee != PEL::bad_index() )
   {
      DFS[ee].push_back( ff ) ;
      FIELD_2_ELMS[idx] = ee ;
   }
   else
   {
      FIELD_2_ELMS[idx] = ELMS.size() ;
      ELMS.push_back( elm ) ;
      DFS.push_back( std::vector< PDE_DiscreteField const* >( 1, ff ) ) ;
   }
   
   PEL_CHECK_POST( has_discretization( ff ) ) ;
   PEL_CHECK_POST( index_of_element( elm ) != PEL::bad_index() ) ;
   PEL_CHECK_POST( reference_element( index_of_element( elm ) ) == elm ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscOnMeshFE:: duplicate_discretization( PDE_DiscreteField const* model_f,
                                             PDE_DiscreteField const* new_f )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: duplicate_discretization" ) ;
   PEL_CHECK_PRE( model_f != 0 ) ;
   PEL_CHECK_PRE( new_f != 0 ) ;
   PEL_CHECK_PRE( has_discretization( model_f ) ) ;

   size_t idx = new_f->id_number() ;
   if( idx >= FIELD_2_ELMS.size() )
   {
      FIELD_2_ELMS.resize( idx+1, PEL::bad_index() ) ;
   }
   size_t const ee = index_of_reference_element( model_f ) ;
   FIELD_2_ELMS[idx] = ee ;
   DFS[ee].push_back( new_f ) ;
   
   PEL_CHECK_POST( has_discretization( new_f ) ) ;
   PEL_CHECK_POST( index_of_reference_element( new_f ) ==
                              index_of_reference_element( model_f ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_DiscOnMeshFE:: has_discretization( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: has_discretization" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   return( index_of_reference_element( ff ) != PEL::bad_index() ) ;
}

//------------------------------------------------------------------------
size_t
PDE_DiscOnMeshFE:: index_of_reference_element( 
                                      PDE_DiscreteField const* ff ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: index_of_reference_element" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   size_t const idx = ff->id_number() ;
   size_t result =
      ( idx < FIELD_2_ELMS.size() ? FIELD_2_ELMS[idx] : PEL::bad_index() ) ;   
   
   PEL_CHECK_POST( EQUIVALENT( has_discretization( ff ),
                               result != PEL::bad_index() ) ) ;
   PEL_CHECK_POST( EQUIVALENT( has_discretization( ff ),
                               result < nb_reference_elements() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DiscOnMeshFE:: nb_discrete_fields( size_t ee ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: nb_discrete_fields" ) ;
   PEL_CHECK_PRE( ee < nb_reference_elements() ) ;
   
   return( (DFS[ ee ]).size() ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_DiscOnMeshFE:: discrete_field( size_t ee, size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: discrete_field" ) ;
   PEL_CHECK_PRE( ee < nb_reference_elements() ) ;
   PEL_CHECK_PRE( i < nb_discrete_fields( ee ) ) ;
   
   PDE_DiscreteField const* result = DFS[ee][i] ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DiscOnMeshFE:: index_of_element( PDE_ReferenceElement const* elm ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: index_of_element" ) ;

   size_t result = PEL::bad_index() ;
   for( size_t ee=0 ; ee<ELMS.size() ; ++ee )
   {
      if( ELMS[ee] == elm )
      {
         result = ee ;
         break ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DiscOnMeshFE:: nb_reference_elements( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: nb_reference_elements" ) ;

   return( ELMS.size() ) ;
}

//-----------------------------------------------------------------------
PDE_ReferenceElement const*
PDE_DiscOnMeshFE:: reference_element( size_t ee ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: reference_element" ) ;
   PEL_CHECK_PRE( ee < nb_reference_elements() ) ;

   PDE_ReferenceElement const* result = ELMS[ ee ] ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DiscOnMeshFE:: nb_basis_functions( size_t ee ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscOnMeshFE:: nb_basis_functions" ) ;
   PEL_CHECK_PRE( ee < nb_reference_elements() ) ;

   size_t result = ELMS[ee]->nb_nodes() ;

   PEL_CHECK_POST( result == reference_element( ee )->nb_nodes() ) ;
   return( result ) ;
}
