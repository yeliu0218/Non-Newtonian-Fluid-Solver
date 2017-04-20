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

#include <GE_Customized_QR.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>

#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

//-----------------------------------------------------------------------------
GE_Customized_QR*
GE_Customized_QR:: create( PEL_Object* a_owner,
                           GE_ReferencePolyhedron const* a_poly,
                           size_t a_order )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Customized_QR:: create" ) ;

   GE_Customized_QR* result =new GE_Customized_QR( a_owner, a_poly, a_order ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_points() == 0 ) ;
   PEL_CHECK_POST( result->order() == a_order ) ;
   PEL_CHECK_POST( result->reference_polyhedron() == a_poly ) ;
   PEL_CHECK_POST( !result->has_been_finalized() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Customized_QR:: GE_Customized_QR( PEL_Object* a_owner,
                                     GE_ReferencePolyhedron const* a_poly,
                                     size_t a_order )
//-----------------------------------------------------------------------------
   : GE_QuadratureRule( a_owner, a_poly, a_order )
   , FINALIZED( false )
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Customized_QR:: ~GE_Customized_QR( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
void
GE_Customized_QR:: reset_points( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Customized_QR:: reset_points" ) ;
   PEL_CHECK_PRE( has_been_finalized() ) ;

   clear_points() ;
   FINALIZED = false ;
   
   PEL_CHECK_POST( nb_points() == 0 ) ;
   PEL_CHECK_POST( !has_been_finalized() ) ;
}

//-----------------------------------------------------------------------------
void
GE_Customized_QR:: insert_point( GE_Point const* pt, double pt_weight )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Customized_QR:: insert_point" ) ;
   PEL_CHECK_PRE( !has_been_finalized() ) ;
   PEL_CHECK_PRE( pt != 0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates()==reference_polyhedron()->dimension() ) ;
   PEL_CHECK_PRE( reference_polyhedron()->contains( pt ) ) ;
   PEL_CHECK_PRE( FORALL( ( size_t i=0; i<nb_points(); ++i ),
                          !pt->is_equal( point(i) ) ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, nb_points, nb_points() ) ;

   append_point( pt->create_clone( this ), pt_weight ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( nb_points() == OLD(nb_points) + 1 ) ;
   PEL_CHECK_POST( point( nb_points()-1 )->is_equal( pt ) ) ;
   PEL_CHECK_POST( point( nb_points()-1 )->owner() == this ) ;
   PEL_CHECK_POST( weight( nb_points()-1 ) == pt_weight ) ;
}

//-----------------------------------------------------------------------------
void
GE_Customized_QR:: finalize( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Customized_QR:: finalize" ) ;
   PEL_CHECK_PRE( !has_been_finalized() ) ;

   double sum = 0.0 ;
   for( size_t i=0 ; i<nb_points() ; ++i )
   {
      sum += weight( i ) ;
   }
   set_sum_of_weights( sum ) ;

   FINALIZED = true ;

   PEL_CHECK_POST( has_been_finalized() ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Customized_QR:: has_been_finalized( void ) const
//-----------------------------------------------------------------------------
{
   return( FINALIZED ) ;
}
