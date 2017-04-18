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

#include <GE_Vector.hh>

#include <GE_Point.hh>

#include <PEL.hh>
#include <PEL_Error.hh>

#include <iostream>

using std::endl ;

//----------------------------------------------------------------------
GE_Vector*
GE_Vector:: create( PEL_Object* a_owner, size_t a_dimension )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: create( a_owner, a_dimension )" ) ;
   PEL_CHECK_PRE( a_dimension>0 ) ;
   GE_Vector* result = new GE_Vector( a_owner, a_dimension ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_components()==a_dimension ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_components() ; ++i ),
                           result->component(i)==0. ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Vector:: GE_Vector( PEL_Object* a_owner, size_t a_dimension )
//----------------------------------------------------------------------
   : PEL_Object( a_owner ),
     xx( a_dimension )
{
   PEL_LABEL( "GE_Vector:: GE_Vector( a_owner, a_dimension )" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Vector*
GE_Vector:: create( PEL_Object* a_owner,
                    GE_Point const* end,
                    GE_Point const* start )
//----------------------------------------------------------------------
{
   PEL_LABEL( " GE_Vector:: create( a_owner, end, start )" ) ;
   PEL_CHECK_PRE( end!=0 && start!=0 ) ;
   PEL_CHECK_PRE( end->nb_coordinates()==start->nb_coordinates() ) ;
   GE_Vector* result = new GE_Vector( a_owner, end, start ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_components()==start->nb_coordinates() &&
                   result->nb_components()==end->nb_coordinates() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_components() ; ++i ),
                           FORMAL( result->component(i)==
                                      end->coordinate( i )
                                            -start->coordinate( i ) ) ) ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Vector:: GE_Vector( PEL_Object* a_owner,
                       GE_Point const* end,
                       GE_Point const* start )
//----------------------------------------------------------------------
   : PEL_Object( a_owner ),
     xx( end->nb_coordinates() ) 
{
   PEL_LABEL( " GE_Vector:: GE_Vector( a_owner, end, start )" ) ;
   for( size_t i=0 ; i<xx.size() ; i++ )
   {
      xx( i ) = end->coordinate( i ) - start->coordinate( i ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Vector*
GE_Vector:: create( PEL_Object* a_owner, doubleVector const& vector )
//----------------------------------------------------------------------
{
   PEL_LABEL( " GE_Vector:: create( a_owner, vector )" ) ;
   PEL_CHECK_PRE( vector.size()>0 ) ;
   GE_Vector* result = new GE_Vector( a_owner, vector ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_components()==vector.size() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_components() ; ++i ),
                           result->component(i)==vector(i) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Vector:: GE_Vector( PEL_Object* a_owner, doubleVector const& vector )
//----------------------------------------------------------------------
   : PEL_Object( a_owner ),
     xx( vector )
{
   PEL_LABEL( " GE_Vector:: GE_Vector( a_owner, vector )" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Vector:: ~GE_Vector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( " GE_Vector:: ~GE_Vector" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Vector*
GE_Vector:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( " GE_Vector:: create_clone" ) ;
   GE_Vector* result = new GE_Vector( a_owner, this ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_components()==nb_components() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_components() ; ++i ),
                           result->component(i)==component(i) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Vector:: GE_Vector( PEL_Object* a_owner, GE_Vector const* other )
//----------------------------------------------------------------------
   : PEL_Object( a_owner ),
     xx( other->xx )
{
   PEL_LABEL( " GE_Vector:: GE_Vector( a_owner, other )" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
GE_Vector:: is_equal( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Vector const* oo = static_cast<GE_Vector const*>( other ) ;

   bool result = ( nb_components() == oo->nb_components() ) ;
   for( size_t i=0 ; result && i<xx.size(); ++i )
   {
      result = result && 
               ( PEL::abs( xx( i ) - oo->component( i ) ) < 1.E-08 ) ;
   }

//?????????? pas de postcondition car le hash_code n'est pas implemente
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
GE_Vector:: hash_code( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: hash_code" ) ;
   PEL_Error::object()->raise_not_implemented( this, "hash_code" ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
size_t
GE_Vector:: nb_components( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: nb_components" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( xx.size() ) ;
}

//----------------------------------------------------------------------
double
GE_Vector:: component( size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: component" ) ;
   PEL_CHECK_PRE( ic<nb_components() ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( xx( ic ) ) ;
}

//----------------------------------------------------------------------
doubleVector const&
GE_Vector:: component_vector( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: component_vector" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( xx ) ;
}

//----------------------------------------------------------------------
double
GE_Vector:: norm( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: norm" ) ;
   PEL_CHECK_INV( invariant() ) ;
   double nrm = 0.0 ;
   for( size_t iDir=0 ; iDir<xx.size() ; iDir++ )
   {
      nrm += xx(iDir)*xx(iDir) ;
   }
   return( PEL::sqrt(nrm) ) ;
}

//----------------------------------------------------------------------
double
GE_Vector:: dot_product( GE_Vector const* vec ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: dot_product" ) ;
   PEL_CHECK_PRE( nb_components() == vec->nb_components() ) ;
   PEL_CHECK_INV( invariant() ) ;

   double scal = 0.0 ;
   for( size_t iDir=0 ; iDir<xx.size() ; iDir++ )
   {
      scal += xx(iDir) * vec->xx(iDir) ;
   }
   return( scal ) ;
}

//----------------------------------------------------------------------
double
GE_Vector:: cosine( GE_Vector const* vec ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: cosine" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( nb_components() == vec->nb_components() ) ;
   
   return( dot_product( vec ) / norm() / vec->norm() ) ;
}

//----------------------------------------------------------------------
void
GE_Vector:: re_initialize( GE_Point const* end, GE_Point const* start )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: re_initialize" ) ;
   PEL_CHECK_PRE( start!=0 ) ;
   PEL_CHECK_PRE( end!=0 ) ;
   PEL_CHECK_PRE( start->nb_coordinates() == end->nb_coordinates() ) ;
   PEL_CHECK_INV( invariant() ) ;

   xx.re_initialize( start->nb_coordinates() ) ;
   for( size_t i=0 ; i<xx.size() ; i++ )
   {
      xx( i ) = end->coordinate( i ) - start->coordinate( i ) ;
   }

   PEL_CHECK_POST( nb_components()==start->nb_coordinates() &&
                   nb_components()==end->nb_coordinates() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_components() ; i++ ),
                           FORMAL( component(i)==
                                     end->coordinate( i )
                                          -start->coordinate( i ) ) ) ) ;   
}

//----------------------------------------------------------------------
void
GE_Vector:: copy( GE_Vector const* vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: copy" ) ;
   PEL_CHECK_PRE( vec!= 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   xx = vec->xx ;

   PEL_CHECK_POST( nb_components()==vec->nb_components() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_components() ; i++ ),
                           component(i)==vec->component(i) ) ) ;   
}

//----------------------------------------------------------------------
void
GE_Vector:: set_component( size_t ic, double x )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: set_component" ) ;
   PEL_CHECK_PRE( ic<nb_components() ) ;
   PEL_CHECK_INV( invariant() ) ;
   xx( ic ) =  x ;
   PEL_CHECK_POST(
      FORMAL(
         FORALL( (size_t i=0 ; i<nb_components() ; ++i ),
                 ( IMPLIES( i!=ic, component(i)==OLD(component(i)) ) ) &&
                 ( IMPLIES( i==ic, component(ic)==x ) ) ) ) ) ;
}

//----------------------------------------------------------------------
void
GE_Vector:: set( GE_Vector const* vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: set" ) ;
   PEL_CHECK_PRE( nb_components()==vec->nb_components() ) ;
   PEL_CHECK_INV( invariant() ) ;
   xx = vec->xx ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_components() ; i++ ),
                           component(i)==vec->component(i) ) ) ;      
}

//----------------------------------------------------------------------
void
GE_Vector:: set_as_cross_product( GE_Vector const* vec1, 
                                  GE_Vector const* vec2 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: set_as_cross_product" ) ;
   PEL_CHECK_PRE( nb_components()==3 ) ;
   PEL_CHECK_PRE( vec1 !=0 ) ;
   PEL_CHECK_PRE( vec2 !=0 ) ;
   PEL_CHECK_PRE( vec1->nb_components()>1 ) ;
   PEL_CHECK_PRE( vec1->nb_components()==vec2->nb_components() ) ;
   PEL_CHECK_INV( invariant() ) ;   

   double x0 = vec1->component( 0 ), y0 = vec1->component( 1 ) ;
   double x1 = vec2->component( 0 ), y1 = vec2->component( 1 ) ;
   set_component( 0,            0. ) ;
   set_component( 1,            0. ) ;
   set_component( 2, x0*y1 - y0*x1 ) ;

   if( vec1->nb_components() == 3 )
   {
      double z0 = vec1->component( 2 ), z1 = vec2->component( 2 ) ;
      set_component( 0, y0*z1 - z0*y1 ) ;
      set_component( 1, z0*x1 - x0*z1 ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
GE_Vector:: scale( double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: scale" ) ;
   PEL_CHECK_INV( invariant() ) ;   
   PEL_SAVEOLD( doubleVector,
                component,
                doubleVector( nb_components() ) ;
                for( size_t i=0 ; i<nb_components() ; ++i )
                {
                   OLD(component(i)) = component(i) ;
                } ) ;
   for( size_t i=0 ; i<xx.size() ; i++ )
   {
      xx(i) *= alpha ;
   }
   PEL_CHECK_POST(
      FORALL( (size_t i=0 ; i<nb_components() ; ++i ),
              FORMAL( component(i)==alpha*OLD(component(i)) ) ) ) ;
   
}

//----------------------------------------------------------------------
void
GE_Vector:: sum( double alpha, GE_Vector const* vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: sum" ) ;
   PEL_CHECK_PRE( vec!=0 ) ;
   PEL_CHECK_PRE( nb_components()==vec->nb_components() ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_SAVEOLD( doubleVector,
                component,
                doubleVector( nb_components() ) ;
                for( size_t i=0 ; i<nb_components() ; ++i )
                {
                   OLD(component(i)) = component(i) ;
                } ) ;
   for( size_t i=0 ; i<xx.size() ; i++ )
   {
      xx(i) += alpha*vec->component( i )  ;
   }
   PEL_CHECK_POST(
      FORALL( (size_t i=0 ; i<nb_components() ; ++i ),
              FORMAL( component(i)==OLD(component(i))
                                        +alpha*vec->component(i) ) ) ) ;
}

//----------------------------------------------------------------------
void
GE_Vector:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Vector:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string space( indent_width, ' ' ) ;
   os << space << "GE_Vector: ( " << component( 0 ) ;
   size_t i = 1 ;
   while( i<nb_components() )   
   {
      os << " , " << component( i ) ;
      i++ ;
   }
   os << " ) " ;

   PEL_CHECK_INV( invariant() ) ;
}
