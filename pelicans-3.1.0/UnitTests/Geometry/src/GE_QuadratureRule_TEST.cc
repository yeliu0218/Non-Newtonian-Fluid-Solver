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

#include <GE_QuadratureRule_TEST.hh>

#include <GE_Customized_QR.hh>
#include <GE_Point.hh>
#include <GE_ReferenceTriangle.hh>
 
#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::string ;

//----------------------------------------------------------------------------
GE_QuadratureRule_TEST const*
GE_QuadratureRule_TEST::unique_instance = new GE_QuadratureRule_TEST() ;
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
GE_QuadratureRule_TEST:: GE_QuadratureRule_TEST( void )
//----------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_QuadratureRule", "GE_QuadratureRule_TEST" )
{
   PEL_LABEL( "GE_QuadratureRule_TEST:: GE_QuadratureRule_TEST" ) ;

   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
GE_QuadratureRule_TEST:: ~GE_QuadratureRule_TEST( void )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: process_all_tests( void )
//----------------------------------------------------------------------------
{
   PEL_ObjectTest::process_all_tests() ;

   do_extendable_rule_tests() ;
}

//----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   MY_EPS = exp->double_data( "dbl_epsilon" ) ;
   MY_MIN = exp->double_data( "dbl_minimum" ) ;

   string const& type = exp->string_data( "type" ) ;
   stringVector const& names = exp->stringVector_data( "rule_names" ) ;
   for( size_t i=0 ; i<names.size() ; ++i )
   {
      string rule_name = names(i) ;
      if( type == "segment_rule" )
      {
         test_one_segment_rule( rule_name ) ;
      }
      else if( type == "triangle_rule" ) 
      {
         test_one_triangle_rule( rule_name ) ;
      }
      else if( type == "square_rule" )
      {
         test_one_square_rule( rule_name ) ;
      }
      else if( type == "tetrahedron_rule" )
      {
         test_one_tetrahedron_rule( rule_name ) ;
      }
      else if( type == "cube_rule" )
      {
         test_one_cube_rule( rule_name ) ;
      }
      else if( type == "square_product_rule" )
      {
         test_one_square_product_rule( rule_name ) ;
      }
      else if( type == "cube_product_rule" )
      {
         test_one_cube_product_rule( rule_name ) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: test_one_segment_rule( std::string const& rule_name )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   GE_QuadratureRule const* qr = GE_QuadratureRule::object( rule_name )  ;

   // Comparison for the integration of all monomials of degree up to order. 
   size_t order = qr->order() ;

   for( size_t i=0 ; i<order+1 ; ++i )
   {
      double exact = 1./(i+1.) ;
      do_a_1D_comparison( qr, i, exact ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: test_one_triangle_rule( std::string const& rule_name )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   GE_QuadratureRule const* qr = GE_QuadratureRule::object( rule_name )  ;

   // Comparison for the integration of all monomials of degree up to order. 
   size_t order = qr->order() ;

   for( size_t i=0; i<order+1; ++i )
   {
      for( size_t j=0; i+j<order+1; ++j )
      {
         double exact =  factorial( i ) * factorial( j ) /
                         double( factorial( i + j +2 ) ) ;
         do_a_2D_comparison( qr, i, j, exact ) ;
      }
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: test_one_square_rule( std::string const& rule_name )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   GE_QuadratureRule const* qr = GE_QuadratureRule::object( rule_name )  ;

   // Integration of all monomials of degree up to order and comparison.
   size_t order = qr->order() ;

   for( size_t i=0; i<order+1; ++i )
   {
      for( size_t j=0; i+j<order+1; ++j )
      {
         double exact_solution = 1./( i+1. )/( j+1. ) ;
         do_a_2D_comparison( qr, i, j, exact_solution ) ;
      }
   }
}

//----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: test_one_tetrahedron_rule( 
                                                 std::string const& rule_name )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   GE_QuadratureRule const* qr = GE_QuadratureRule::object( rule_name )  ;

   // Integration of all monomials of degree up to order and comparison.
   size_t order = qr->order() ;

   for( size_t i=0; i<order+1; ++i )
   {
      for( size_t j=0; i+j<order+1; ++j )
      {
         for( size_t k=0; i+j+k<order+1; ++k )
	 {
            double exact = factorial(i)*factorial(j)*factorial(k) /
	                   double( factorial( i + j + k + 3 ) ) ;
            do_a_3D_comparison( qr, i, j, k, exact ) ;
	 }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: test_one_cube_rule( std::string const& rule_name )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   GE_QuadratureRule const* qr = GE_QuadratureRule::object( rule_name )  ;

   // Integration of all monomials of degree up to order and comparison.
   size_t order = qr->order() ;

   for( size_t i=0; i<order+1; ++i )
   {
      for( size_t j=0; i+j<order+1; ++j )
      {
         for( size_t k=0; i+j+k<order+1; ++k )
	 {
            double exact =  1./( i + 1. )/( j + 1. )/( k + 1. ) ;
            do_a_3D_comparison( qr, i, j, k, exact ) ;
	 }
      }
   }

   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: test_one_square_product_rule( 
                                               std::string const& rule_name )
//-----------------------------------------------------------------------------
{
   GE_QuadratureRule const* qr = GE_QuadratureRule::object( rule_name )  ;
   
   size_t order = qr->order() ;

   for( size_t i=0; i<order+1; ++i )
   {
      for( size_t j=0; j<order+1; ++j )
      {
         double exact = 1./( i+1. )/( j+1. ) ;
         do_a_2D_comparison( qr, i, j, exact ) ;
      }
   }
}

//-----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: test_one_cube_product_rule( 
                                               std::string const& rule_name )
//-----------------------------------------------------------------------------
{
   GE_QuadratureRule const* qr = GE_QuadratureRule::object( rule_name )  ;
   
   size_t order = qr->order() ;

   for( size_t i=0; i<order+1; ++i )
   {
      for( size_t j=0; j<order+1; ++j )
      {
         for( size_t k=0; k<order+1; ++k )
         {
            double exact = 1./( i+1. )/( j+1. )/( k+1. ) ;
            do_a_3D_comparison( qr, i, j, k, exact ) ;
	 }
      }
   }
}

//----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: do_extendable_rule_tests( void )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   
   // Use a GE_Customized_QR to build a 5 order integration rule on the unit
   // triangle.
   // cf : G. Dhatt and G. Touzot,
   //      une présentation de la méthode des éléments finis, chapter 5.1
   //      Maloine, PARIS, 1989.
   GE_Customized_QR* qr =
      GE_Customized_QR::create( 0, GE_ReferenceTriangle::object(), 5 ) ;

   double const a = ( 6. + PEL::sqrt( 15. ) )/21 ;
   double const b = 4./7.-a ;
   double const c = 1./3. ;
   double const wa = ( 155. + PEL::sqrt( 15. ) )/2400. ;
   double const wb = 31./240. - wa ;
   double const wc = 9./80. ;

   qr->insert_point( GE_Point::create( qr, c,    c    ), wc ) ;
   qr->insert_point( GE_Point::create( qr, a,    a    ), wa ) ;
   qr->insert_point( GE_Point::create( qr, 1.-2.*a, a ), wa ) ;
   qr->insert_point( GE_Point::create( qr, a , 1-2.*a ), wa ) ;
   qr->insert_point( GE_Point::create( qr, b, b       ), wb ) ;
   qr->insert_point( GE_Point::create( qr, 1.-2.*b, b ), wb ) ;
   qr->insert_point( GE_Point::create( qr, b, 1.-2.*b ), wb ) ;

   size_t order = qr->order() ;
   size_t nb_points = qr->nb_points() ;

   // Integration of all monomials of degree up to order and comparison
   bool result = true ;
   for( size_t i=0; i<order+1; i++ )
   {
      for( size_t j=0; i+j<order+1; j++ )
      {
         double exact_solution = factorial( i ) * factorial( j ) /
	                         double( factorial( i + j +2 ) ) ;
         double approximation = 0. ;
         for( size_t l=0 ; l<nb_points; l++ )
         {
	    approximation += 
                    qr->weight( l )*PEL::pow( qr->point( l )->coordinate( 0 ), i )
                                   *PEL::pow( qr->point( l )->coordinate( 1 ), j ) ; 
         }
         result = result && PEL::equal( approximation, exact_solution ) ;
      }
   }

   qr->destroy() ;

   notify_one_test_result( "GE_Customized_QR", result ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: do_a_1D_comparison( GE_QuadratureRule const* qr,
                                              int i, double exact )
//-----------------------------------------------------------------------------
{
   double approx = 0. ;
   for( size_t l=0 ; l<qr->nb_points(); l++ )
   {
      approx += qr->weight( l )
                        * PEL::pow( qr->point( l )->coordinate( 0 ), i ) ;
   }
   std::ostringstream mesg ;
   mesg << qr->name() << "  (x^" << i << ")" ;
   bool eq = PEL::double_equality( approx, exact, MY_EPS, MY_MIN ) ;
   notify_one_test_result( mesg.str(), eq ) ;
   if( !eq ) display_error( approx, exact ) ;
}

//-----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: do_a_2D_comparison( GE_QuadratureRule const* qr,
                                              int i, int j, double exact )
//-----------------------------------------------------------------------------
{
   double approx = 0. ;
   for( size_t l=0 ; l<qr->nb_points(); l++ )
   {
      approx += qr->weight( l )
                        * PEL::pow( qr->point( l )->coordinate( 0 ), i )
        	        * PEL::pow( qr->point( l )->coordinate( 1 ), j ) ;
   }
   std::ostringstream mesg ;
   mesg << qr->name() << "  (x^" << i << ")(y^" << j << ")" ;
   bool eq = PEL::double_equality( approx, exact, MY_EPS, MY_MIN ) ;
   notify_one_test_result( mesg.str(), eq ) ;
   if( !eq ) display_error( approx, exact ) ;
}

//-----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: do_a_3D_comparison( GE_QuadratureRule const* qr,
                                             int i, int j, int k,
                                             double exact )
//-----------------------------------------------------------------------------
{
   double approx = 0. ;
   for( size_t l=0 ; l<qr->nb_points(); l++ )
   {
      approx += qr->weight( l )
                        * PEL::pow( qr->point( l )->coordinate( 0 ), i )
        	        * PEL::pow( qr->point( l )->coordinate( 1 ), j )
        	        * PEL::pow( qr->point( l )->coordinate( 2 ), k ) ;
   }
   std::ostringstream mesg ;
   mesg << qr->name() << "  (x^" << i << ")(y^" << j << ")(z^" << k << ")" ;
   bool eq = PEL::double_equality( approx, exact, MY_EPS, MY_MIN ) ;
   notify_one_test_result( mesg.str(), eq ) ;
   if( !eq ) display_error( approx, exact ) ;
}

//----------------------------------------------------------------------------
size_t
GE_QuadratureRule_TEST:: factorial( size_t i ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( i < 20 ) ;

   size_t f = 1 ;
   for( size_t j=2; j<i+1; j++ )
   {
      f *= j ;
   }

   return( f ) ;
}

//----------------------------------------------------------------------------
void
GE_QuadratureRule_TEST:: display_error( double approx, double exact ) const
//----------------------------------------------------------------------------
{
   std::cout << " approx=" << approx << " exact=" << exact << std::endl ;
   std::cout << " 1-approx/exact=" << approx/exact-1.0 << std::endl ;
}

//----------------------------------------------------------------------------
bool
GE_QuadratureRule_TEST:: invariant( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_ObjectTest::invariant() ) ;

   return true ;
}


