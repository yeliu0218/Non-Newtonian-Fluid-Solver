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

#include <GE_QuadratureRule.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <doubleVector.hh>

#include <iostream>

//-----------------------------------------------------------------------------
GE_QuadratureRule const*
GE_QuadratureRule:: object( std::string a_name )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: object" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   GE_QuadratureRule const* result =
      static_cast<GE_QuadratureRule const*>( plugins_map()->item( a_name ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( PEL_Root::object() ) ) ;   
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( result->sum_of_weights() ==
                                  result->reference_polyhedron()->measure() ) ;

   return( result ) ;
}

//----------------------------------------------------------------------------
GE_QuadratureRule:: GE_QuadratureRule( PEL_Object* a_owner,
                                       GE_ReferencePolyhedron const* poly,
                                       size_t a_order )
//----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , RULE_NAME( "GE_Customized_QR" )
   , POLY( poly )
   , ORDER( a_order )
   , POINTS( PEL_Vector::create( this, 0 ) )
   , WEIGHTS( 0 )
   , TOT_WEIGHT( -1.0 )
{
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_QuadratureRule:: GE_QuadratureRule( std::string a_name,
                                       GE_ReferencePolyhedron const* poly,
                                       size_t a_order )
//-----------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , RULE_NAME( a_name )
   , POLY( poly )
   , ORDER( a_order )
   , POINTS( PEL_Vector::create( this, 0 ) )
   , WEIGHTS( 0 )
   , TOT_WEIGHT( -1.0 )
{
   PEL_LABEL( "GE_QuadratureRule:: GE_QuadratureRule" ) ;

   plugins_map()->register_item( a_name, this ) ;
   
   PEL_CHECK_INV( invariant() ) ;   
}

//-----------------------------------------------------------------------------
GE_QuadratureRule:: ~GE_QuadratureRule( void )
//-----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
std::string const&
GE_QuadratureRule:: name(void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: name" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( RULE_NAME ) ;
}

//----------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_QuadratureRule:: reference_polyhedron( void ) const
//----------------------------------------------------------------------------
{
   GE_ReferencePolyhedron const* result = POLY ;
   PEL_CHECK_POST( result!=0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
size_t
GE_QuadratureRule:: order( void ) const
//----------------------------------------------------------------------------
{
   return( ORDER ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_QuadratureRule:: nb_points( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: nb_points" ) ;
   PEL_CHECK_INV( invariant() ) ;   

   return( WEIGHTS.size() ) ;
}

//----------------------------------------------------------------------------
GE_Point const*
GE_QuadratureRule:: point( size_t i ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: point" ) ;
   PEL_CHECK_PRE( i<nb_points() ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_CHECK( dynamic_cast<GE_Point*>( POINTS->at( i ) ) != 0 ) ;
   GE_Point const* result = static_cast<GE_Point*>( POINTS->at( i ) ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( reference_polyhedron()->contains( result ) ) ; 
   return( result ) ;
}

//-----------------------------------------------------------------------------
double
GE_QuadratureRule:: weight( size_t i ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: weight" ) ;
   PEL_CHECK_PRE( i<nb_points() ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( WEIGHTS( i ) ) ;
}

//-----------------------------------------------------------------------------
double
GE_QuadratureRule:: sum_of_weights( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: sum_of_weights" ) ;
   PEL_CHECK_INV( invariant() ) ;

   double result = 0.0 ;
   if( TOT_WEIGHT > 0.0 )
   {
      result = TOT_WEIGHT ;
   }
   else
   {
      for( size_t i=0 ; i<WEIGHTS.size() ; ++i )
      {
         result += WEIGHTS( i ) ;
      }
   }

   PEL_CHECK_POST( sum_of_weights_POST( result ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_QuadratureRule:: print( std::ostream& os, size_t indent_width ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: print" ) ;

   std::string space( indent_width, ' ' ) ;
   os << space << " Quadrature rule \"" << RULE_NAME << "\" :" << std::endl ;
   for( size_t i=0 ; i<nb_points() ; i++ )
   {
      os << space << "   Pt("<<i<<") : " << point( i ) << " weight : " 
         << weight( i ) << std::endl ;
   }
}

//-----------------------------------------------------------------------------
void
GE_QuadratureRule:: append_point( GE_Point* pt, double pt_weight )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: append_point" ) ;
   PEL_CHECK_PRE( pt!=0 ) ;
   PEL_CHECK_PRE( pt->owner()==this ) ;
   PEL_CHECK_PRE( pt->nb_coordinates()==reference_polyhedron()->dimension() ) ;
   PEL_CHECK_PRE( reference_polyhedron()->contains( pt ) ) ;
   PEL_CHECK_PRE( FORALL( ( size_t i=0; i<nb_points(); ++i ),
                          !pt->is_equal( point(i) ) ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( size_t, nb_points, nb_points() ) ;
   
   POINTS->append( pt ) ; 
   WEIGHTS.append( pt_weight ) ;
   TOT_WEIGHT = -1. ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( nb_points() == OLD(nb_points) + 1 ) ;
   PEL_CHECK_POST( point( nb_points()-1 ) == pt ) ;
   PEL_CHECK_POST( weight( nb_points()-1 ) == pt_weight ) ;
}

//-----------------------------------------------------------------------------
void
GE_QuadratureRule:: set_sum_of_weights( double sum )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: set_sum_of_weights" ) ;
   PEL_CHECK_PRE( set_sum_of_weights_PRE( sum ) ) ;
   PEL_CHECK_INV( invariant() ) ;
                                             
   TOT_WEIGHT = sum ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( sum_of_weights() == sum ) ;
}

//-----------------------------------------------------------------------------
void
GE_QuadratureRule:: clear_points( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QuadratureRule:: clear_points" ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<POINTS->index_limit() ; ++i )
   {
      PEL_CHECK( dynamic_cast<GE_Point*>( POINTS->at( i ) ) != 0 ) ;
      GE_Point* pt = static_cast<GE_Point*>( POINTS->at( i ) ) ;
      destroy_possession( pt ) ;
   }
   POINTS->clear() ;
   WEIGHTS.re_initialize(0) ;
   TOT_WEIGHT = -1. ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( nb_points() == 0 ) ;
}
                  
//-----------------------------------------------------------------------------
bool
GE_QuadratureRule:: set_sum_of_weights_PRE( double sum ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( sum > 0. ) ;
   double sum_weight_i = 0.0 ;
   for( size_t i=0 ; i<nb_points() ; ++i ) sum_weight_i+=weight(i) ;
   PEL_ASSERT( PEL::abs(sum_weight_i-sum)<1.E-8 ) ;
   return( true ) ;
}
 
                  
//-----------------------------------------------------------------------------
bool
GE_QuadratureRule:: sum_of_weights_POST( double result ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( result > 0. ) ;
   double sum_weight_i = 0.0 ;
   for( size_t i=0 ; i<nb_points() ; ++i ) sum_weight_i+=weight(i) ;
   PEL_ASSERT( PEL::abs(sum_weight_i-result)<1.E-8 ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------------
bool
GE_QuadratureRule:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( IMPLIES( POINTS != 0, POINTS->count() == WEIGHTS.size() ) ) ;
   double sum_weight_i = 0.0 ;
   for( size_t i=0 ; i<WEIGHTS.size() ; ++i ) sum_weight_i+=WEIGHTS(i) ;
   PEL_ASSERT( IMPLIES( TOT_WEIGHT>0 , PEL::abs(sum_weight_i-TOT_WEIGHT)<1.E-8 ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
GE_QuadratureRule:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "GE_QuadratureRule descendant" ) ;
   return( result ) ;
}
