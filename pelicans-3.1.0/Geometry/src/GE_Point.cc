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

#include <GE_Point.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Vector.hh>

#include <PEL_Error.hh>
#include <PEL_DoubleComparatorFloat.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <PEL.hh>

#include <iomanip>
#include <iostream>

using std::ostream ;

#ifdef OUTLINE
#include <GE_Point.icc>
#endif
   
//----------------------------------------------------------------------
GE_Point*
GE_Point:: create( PEL_Object* a_owner,
                   size_t a_dimension )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: create(a_owner,a_dimension)" ) ;
   GE_Point* result = new GE_Point( a_owner, a_dimension ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_coordinates()==a_dimension ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<a_dimension ; ++i ),
              result->coordinate(i)==0. ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point:: GE_Point( PEL_Object* a_owner,
                     size_t a_dimension )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , xx( a_dimension )
{
   PEL_LABEL( "GE_Point:: GE_Point(a_owner,dimension) ") ;
   xx.set( 0.0 ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Point*
GE_Point:: create( PEL_Object* a_owner,
                   doubleVector const& a_coordinates_table )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: create(a_owner,a_coordinates_table)" ) ;
   GE_Point* result = new GE_Point( a_owner, a_coordinates_table ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_coordinates()==a_coordinates_table.size() ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<a_coordinates_table.size() ; ++i ),
              result->coordinate(i)==a_coordinates_table(i) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point:: GE_Point( PEL_Object* a_owner, 
                     doubleVector const& a_coordinates_table )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , xx( a_coordinates_table )
{
   PEL_LABEL( "GE_Point:: GE_Point(a_owner,a_coordinates_table ) ") ;
   PEL_CHECK_INV( invariant() ) ;   
}

//----------------------------------------------------------------------
GE_Point*
GE_Point:: create( PEL_Object* a_owner,
                   double x0 ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: create(a_owner,x0)" ) ;
   GE_Point* result = new GE_Point( a_owner, x0 ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_coordinates()==1 ) ;
   PEL_CHECK_POST( result->coordinate(0)==x0 ) ;
   return( result ) ;   
}

//----------------------------------------------------------------------
GE_Point:: GE_Point( PEL_Object* a_owner,
                     double x0 )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , xx( 1 )
{
   PEL_LABEL( "GE_Point:: GE_Point(a_owner,x0)" ) ;
   xx( 0 ) = x0 ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Point*
GE_Point:: create( PEL_Object* a_owner, double x0, double x1 ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: create(a_owner,x0,x1)" ) ;
   GE_Point* result = new GE_Point( a_owner, x0, x1 ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_coordinates()==2 ) ;
   PEL_CHECK_POST( result->coordinate(0)==x0 &&
                   result->coordinate(1)==x1 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point:: GE_Point( PEL_Object* a_owner, double x0, double x1 )
//----------------------------------------------------------------------
   : PEL_Object( a_owner ),
     xx( 2 )
{
   PEL_LABEL( "GE_Point:: GE_Point(a_owner,x0,x1)" ) ;
   xx( 0 ) = x0 ;
   xx( 1 ) = x1 ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Point*
GE_Point:: create( PEL_Object* a_owner,
                   double x0, double x1, double x2 ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: create(a_owner,x0,x1,x2)" ) ;
   GE_Point* result = new GE_Point( a_owner, x0, x1, x2 ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->nb_coordinates()==3 ) ;
   PEL_CHECK_POST( result->coordinate(0)==x0 &&
                   result->coordinate(1)==x1 &&
                   result->coordinate(2)==x2 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point:: GE_Point( PEL_Object* a_owner,
                     double x0, double x1, double x2 )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , xx( 3 )
{
   PEL_LABEL( "GE_Point:: GE_Point(a_owner,x0,x1,x2)" ) ;
   xx( 0 ) = x0 ;
   xx( 1 ) = x1 ;
   xx( 2 ) = x2 ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Point:: ~GE_Point( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: ~GE_Point" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Point const*
GE_Point:: origin( size_t dimension )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: origin" ) ;

   static PEL_Vector* vecOrigin =
                           PEL_Vector::create( PEL_Root::object(), 0 ) ;
   PEL_CHECK( vecOrigin!=0 ) ;

   if( vecOrigin->index_limit()<dimension+1 )
   {
      vecOrigin->resize( dimension+1 ) ;
   }
   if( vecOrigin->at( dimension ) == 0 )
   {
      GE_Point* pt0 = GE_Point::create( PEL_Root::object(), dimension ) ;
      vecOrigin->set_at( dimension, pt0 ) ;
   }

   PEL_CHECK( dynamic_cast<GE_Point*>( vecOrigin->at( dimension ) ) != 0 ) ;
   GE_Point* result = static_cast<GE_Point*>( vecOrigin->at( dimension ) ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==PEL_Root::object() ) ;
   PEL_CHECK_POST( result->nb_coordinates()==dimension ) ;
   PEL_CHECK_POST( FORALL( (size_t i=0 ; i<dimension ; ++i ),
                           result->coordinate(i)==0. ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point*
GE_Point:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: create_clone" ) ;
   GE_Point* result = new GE_Point( a_owner, this ) ;
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   PEL_CHECK_POST( result->nb_coordinates()==nb_coordinates() ) ;
   PEL_CHECK_POST( FORALL( (size_t i=0 ; i<nb_coordinates() ; ++i ),
                           result->coordinate(i)==coordinate(i) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point:: GE_Point( PEL_Object* a_owner, GE_Point const* other )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , xx( other->xx )
{
   PEL_LABEL( "GE_Point:: GE_Point(a_owner,other)" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_DoubleComparator const*
GE_Point:: coordinates_comparator( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: coordinates_comparator" ) ;
   
   static PEL_DoubleComparator const* result =
      PEL_DoubleComparatorFloat::create( PEL_Root::object(), 1.E-12 ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
GE_Point:: is_equal( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   static PEL_DoubleComparator const* dbl_cmp = coordinates_comparator() ;
   
   GE_Point const* pt = static_cast<GE_Point const*>( other ) ;
   size_t const nb_coords = nb_coordinates() ;
   bool result = ( nb_coords==pt->nb_coordinates() ) ;
   if( result )
   {
      for( size_t i = 0 ; i<nb_coords ; ++i )
      {
         if( dbl_cmp->three_way_comparison( xx(i), pt->xx(i) ) != 0 )
         {
            result = false ;
            break ;
         }
      }
   }

//??? has_code pas implemente   PEL_CHECK_POST( is_equal_POST( other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
int
GE_Point:: three_way_comparison( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   static PEL_DoubleComparator const* dbl_cmp = coordinates_comparator() ;
   
   GE_Point const* pt = static_cast<GE_Point const*>( other ) ;
   size_t const nb_coords = nb_coordinates() ;
   size_t const nb_coords2 = pt->nb_coordinates() ;
   int result = 0 ;

   if( nb_coords == nb_coords2 )
   {
      for( size_t i = 0 ; i<nb_coords && result==0 ; ++i )
      {
         result = dbl_cmp->three_way_comparison( xx(i), pt->xx(i) ) ;
      }
   }
   else if( nb_coords < nb_coords2 )
   {
      result = -1 ;
   }
   else
   {
      result = 1 ;
   }

   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
GE_Point:: hash_code( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: hash_code" ) ;
   PEL_Error::object()->raise_not_implemented( this, "hash_code" ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
double
GE_Point:: distance( GE_Point const* other ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: distance" ) ;
   PEL_CHECK_PRE( nb_coordinates() == other->nb_coordinates()  ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   double dist = 0. ;
   for( size_t i=0 ; i<nb_coordinates() ; ++i )
   {
      double const x = xx(i) - other->xx(i) ;
      dist += x*x ;
   }
   return( PEL::sqrt(dist) ) ;
}

//----------------------------------------------------------------------
void
GE_Point:: copy( GE_Point const* pt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: copy" ) ;
   PEL_CHECK_PRE( pt!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   xx = pt->xx ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( nb_coordinates()==pt->nb_coordinates() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_coordinates() ; ++i ),
                           coordinate(i)==pt->coordinate(i) ) ) ;
}

//----------------------------------------------------------------------
void
GE_Point:: set( GE_Point const* pt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: set" ) ;
   PEL_CHECK_PRE( pt!=0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates()==nb_coordinates() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   xx = pt->xx ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_coordinates() ; ++i ),
                           coordinate(i)==pt->coordinate(i) ) ) ;
}

//----------------------------------------------------------------------
void
GE_Point:: set_coordinates( doubleVector const& a_coordinates_table )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: set_coordinates" ) ;
   PEL_CHECK_PRE( nb_coordinates()==a_coordinates_table.size() ) ;
   PEL_CHECK_INV( invariant() ) ;
   xx = a_coordinates_table ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_coordinates() ; ++i ),
                           coordinate(i)==a_coordinates_table(i) ) ) ;
}

//----------------------------------------------------------------------
void
GE_Point:: set_as_barycenter( double alpha,
                              GE_Point const* ptA,
                              GE_Point const* ptB )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: set_as_barycenter" ) ;
   PEL_CHECK_PRE( ptA!=0 ) ;
   PEL_CHECK_PRE( nb_coordinates() == ptA->nb_coordinates() ) ;
   PEL_CHECK_PRE( ptB!=0 ) ;
   PEL_CHECK_PRE( nb_coordinates() == ptB->nb_coordinates() ) ;
   PEL_CHECK_INV( invariant() ) ;
   for( size_t i=0 ; i<nb_coordinates() ; ++i )
   {
      xx(i) = (1.0-alpha)*ptA->xx(i) + alpha*ptB->xx(i) ;
   }
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST(
      FORALL( (size_t i=0 ; i<nb_coordinates() ; ++i ),
              FORMAL( coordinate(i)==
                           (1.0-alpha)*ptA->coordinate(i)
                                        +alpha*ptB->coordinate(i) ) ) ) ;
}

//----------------------------------------------------------------------
void
GE_Point:: translate( double alpha, GE_Vector const* vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: translate" ) ;
   PEL_CHECK_PRE( vec!=0 ) ;
   PEL_CHECK_PRE( vec->nb_components() == nb_coordinates() ) ;
   PEL_SAVEOLD( doubleVector,
                coordinate,
                doubleVector( nb_coordinates() ) ;
                for( size_t i=0 ; i<nb_coordinates() ; ++i )
                {
                   OLD(coordinate(i)) = coordinate(i) ;
                } ) ;
   PEL_CHECK_INV( invariant() ) ;
   for( size_t i=0 ; i<nb_coordinates() ; ++i )
   {
      xx(i) += alpha * vec->component(i) ;
   }
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST(
      FORALL( (size_t i=0 ; i<nb_coordinates() ; ++i ),
              FORMAL( coordinate(ic)==
                              OLD(coordinate(i)
                                     +alpha*vec->component(i) ) ) ) ) ;
   
}

//----------------------------------------------------------------------
void
GE_Point:: print( std::ostream& os, size_t indent_width ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Point:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   std::string space( indent_width, ' ' ) ;
   std::ios_base::fmtflags original_flags = os.flags() ;
   os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
   std::streamsize p = os.precision() ;
   os << std::setprecision(9) ;
   os << space << "( " ;
   for( size_t i=0 ; i<nb_coordinates() ; ++i )
   {
      os << xx(i);
      if( i<nb_coordinates()-1 )
      {
         os << " , " ;
      }
   }
   os << " )" ;
   os << std::setprecision(p) ;
   os.flags( original_flags ) ;
}
