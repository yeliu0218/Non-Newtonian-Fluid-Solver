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

#include <GE_Point_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_List.hh> 
#include <PEL_Vector.hh> 
#include <doubleVector.hh> 

#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PEL_ModuleExplorer.hh>
 
#include <iostream>

using std::cout ;
using std::endl ;
using std::string ;

//-------------------------------------------------------------------------
GE_Point_TEST*
GE_Point_TEST::unique_instance = new GE_Point_TEST() ;
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
GE_Point_TEST:: GE_Point_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_Point", "GE_Point_TEST" )
{
}


//----------------------------------------------------------------------------
GE_Point_TEST:: ~GE_Point_TEST( void )
//----------------------------------------------------------------------------
{
}


//----------------------------------------------------------------------------
void
GE_Point_TEST:: process_all_tests( void )
//----------------------------------------------------------------------------
{
   size_t dim = data_deck_explorer()->int_data( "dimension" ) ;
   PEL_ASSERT( dim>0 && dim<4 ) ;

   doubleVector const& coords =
                    data_deck_explorer()->doubleVector_data( "coordinates" ) ;
   PEL_ASSERT( coords.size() == dim ) ;
   
   // Test GE_Point::create
   GE_Point const* point = GE_Point::create( 0, coords ) ;
   bool ok = true ;
   for( size_t i=0 ; ok && i<dim ; ++i )
   {
      ok = ok && ( point->coordinate(i) == coords(i) ) ;
   }
   notify_one_test_result( "GE_Point::coordinate", ok ) ;

   // Test GE_Point::nb_coordinates
   notify_one_test_result(
      "GE_Point::nb_coordinates",
      ( dim == point->nb_coordinates() ) ) ;

   GE_Point* point2 = 0 ;

   // Test GE_Point::create from coordinate list
   if( dim==1 )
   {
      point2 = GE_Point::create( 0,
                                 point->coordinate( 0 ) ) ;
   }
   else if( dim==2 )
   {
      point2 = GE_Point::create( 0,
                                 point->coordinate( 0 ),
                                 point->coordinate( 1 ) ) ;
   }
   else if( dim==3 )
   {
      point2 = GE_Point::create( 0,
                                 point->coordinate( 0 ),
                                 point->coordinate( 1 ),
                                 point->coordinate( 2 ) ) ;
   }
   ok = true ;
   for( size_t i=0 ; ok && i<dim ; ++i )
   {
      ok = ok && ( point2->coordinate(i) == coords(i) ) ;
   }
   notify_one_test_result( "GE_Point::create", ok ) ;

   // Test GE_Point::is_equal and GE_Point::set_coordinate
   ok = point2->is_equal( point ) ;
   point2->set_coordinate( 0, point->coordinate(0)+1.) ;
   notify_one_test_result(
      "GE_Point::set_coordinate",
      point2->coordinate(0) == point->coordinate(0)+1. ) ;
   ok = ok && ( !point2->is_equal( point ) ) ;
   notify_one_test_result( "GE_Point::is_equal", ok ) ;

   // Test GE_Point::set_coordinates
   point2->set_coordinates( coords ) ;
   notify_one_test_result( "GE_Point::set_coordinates",
                           point2->is_equal( point ) ) ;

   // Test GE_Point::copy :
   point2->copy( point ) ;
   notify_one_test_result( "GE_Point::copy", point2->is_equal( point ) ) ;

   // Test GE_Point::translate :
   GE_Point const* orig = point2->origin( dim ) ;
   GE_Vector* vector = GE_Vector::create( 0, point2, orig ) ;
   GE_Vector* transf = GE_Vector::create( 0, coords ) ;
   point2->translate( 2., transf ) ;
   notify_one_test_result(
      "GE_Point::translate",
      PEL::abs( point2->distance( orig )-3.*transf->norm() )<1.E-08 ) ;

   // Test GE_Point::set_as_barycenter
   point2->set_as_barycenter( 0.5, point, point2 ) ;
   notify_one_test_result(
      "GE_Point::set_as_barycenter",
      PEL::abs( point2->distance( orig )-2.*transf->norm() )<1.E-08 ) ;

   point2->destroy() ; point2 = 0 ;
   vector->destroy() ; vector = 0 ;
   transf->destroy() ; transf = 0 ;
   point->destroy()  ; point  = 0 ;
}
