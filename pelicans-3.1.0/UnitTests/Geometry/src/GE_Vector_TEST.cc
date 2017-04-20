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

#include <GE_Vector_TEST.hh>

#include <PEL.hh>
#include <PEL_List.hh> 
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh> 

#include <doubleVector.hh>

#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <iostream>

using std::cout ;
using std::endl ;
using std::string ;

//-------------------------------------------------------------------------
GE_Vector_TEST*
GE_Vector_TEST::unique_instance = new GE_Vector_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
GE_Vector_TEST:: GE_Vector_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_Vector", "GE_Vector_TEST" )
{
}


//----------------------------------------------------------------------------
GE_Vector_TEST:: ~GE_Vector_TEST( void )
//----------------------------------------------------------------------------
{
}




//----------------------------------------------------------------------------
void
GE_Vector_TEST:: process_all_tests( void )
//----------------------------------------------------------------------------
{
   size_t dim = data_deck_explorer()->int_data( "dimension" ) ;
   PEL_ASSERT( dim>0 && dim<4 ) ;

   doubleVector const& compVec = data_deck_explorer()->doubleVector_data( "components" ) ;

   GE_Vector* vector = GE_Vector::create( 0, compVec ) ;

   bool ok =  ( dim == vector->nb_components() ) ;

   // Test all creation media
   GE_Vector* vector2 = GE_Vector::create( 0, dim ) ;
   for( size_t i = 0; i<dim; i++ )
   {
      vector2->set_component( i, vector->component( i ) ) ;
   }
   ok = ok && vector->is_equal( vector2 ) ;
   GE_Point* endpoint2 = GE_Point::create( 0, compVec ) ;
   GE_Point const* endpoint1 = endpoint2->origin( dim ) ;
   GE_Vector* vector3 = GE_Vector::create( 0, endpoint2,
                                              endpoint1 ) ;
   ok = ok && vector->is_equal( vector3 ) ;
   vector2->copy( vector ) ;
   ok = ok && vector->is_equal( vector2 ) ;
   vector3->re_initialize( endpoint2, endpoint1 ) ;
   ok = ok && vector->is_equal( vector3 ) ;

   // Test various methods
   if( ok )
   {
      ok = ok && ( PEL::abs( vector->norm() 
                   - PEL::sqrt( vector->dot_product( vector ) ) )<1.E-08 ) ;
      vector2->sum( 1., vector2 ) ;   
      vector2->scale( 0.5 ) ;
      ok = ok && vector->is_equal( vector3 ) ;
   }
   if( ok && ( dim>1 ) )
   {
      double sum = 0. ;
      for( size_t i = 1; i<dim; i++ )
      {
         vector2->set_component( i, vector->component( 0 ) ) ;
         sum += vector->component( i ) ;
      }
      vector2->set_component( 0, -1.*sum ) ;

      ok = ok && ( PEL::abs( vector->cosine( vector2 ) ) < 1.E-08 ) ;

      GE_Vector* vector4 = GE_Vector::create( 0, 3 ) ;
      vector4->set_as_cross_product( vector, vector2 ) ;
      if( vector->nb_components()==2 )
      {
         ok = ok && ( PEL::abs ( PEL::abs( vector4->component( 2 ) ) - vector->norm()*vector2->norm() ) < 1.E-08 ) ;
      }
      else if( vector->nb_components()==3 )
      {
         ok = ok && ( PEL::abs( vector->cosine( vector4 ) ) < 1.E-08 ) ;
      }
      vector4->destroy() ;
   }

   endpoint2->destroy() ;
   vector->destroy() ;
   vector2->destroy() ;
   vector3->destroy() ;

   notify_one_test_result( "Geometrical properties test", ok ) ;
}
