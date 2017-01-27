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

#ifndef GE_POLYGON_TEST_HH
#define GE_POLYGON_TEST_HH

#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectTest.hh>

class GE_Polygon ;
class GE_Polygon2D ;
class GE_SimplePolygon2D ;

/** Unit tests for the `GE_Polygon' class.
*/

class PEL_EXPORT GE_Polygon_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      GE_Polygon_TEST( void ) ;
     ~GE_Polygon_TEST( void ) ;
      GE_Polygon_TEST( GE_Polygon_TEST const& other ) ;
      GE_Polygon_TEST const& operator=( GE_Polygon_TEST const& other ) ;

  //-- Unit test management
      void process_all_tests( void )  ;

      void process_polygon2D_tests( PEL_ModuleExplorer const* exp ) ;
      void process_simple_polygon2D_tests( PEL_ModuleExplorer const* exp ) ;

      void set_geometrical_algorithms( PEL_ModuleExplorer const* exp,
                                       GE_Polygon* polygon ) const ;

   //-- tests on Polygon

      bool test_vertices_chain( PEL_ModuleExplorer const* exp,
                                GE_Polygon const* polygon ) const ;

      bool test_area( PEL_ModuleExplorer const* exp,
                      GE_Polygon const* polygon ) const ;
      
   //-- Specific tests on SimplePolygon2D

      bool test_orientation( PEL_ModuleExplorer const* exp,
                             GE_SimplePolygon2D const* polygon ) const ;

      bool test_triangulation( PEL_ModuleExplorer const* exp,
                               GE_SimplePolygon2D const* polygon ) const ;

      bool test_inclusion( PEL_ModuleExplorer const* exp,
                           GE_SimplePolygon2D const* polygon ) const ;

      bool test_inclusion_on_boundary(
                           PEL_ModuleExplorer const* exp,
                           GE_SimplePolygon2D const* polygon ) const ;

      bool test_simplification( GE_SimplePolygon2D const* polygon ) const ;
      
   //-- Specific tests on Polygon2D

      bool test_split_of_simple_polygons(
                           PEL_ModuleExplorer const* exp,
                           GE_Polygon2D const* polygon ) const ;

   //-- Class attributes

      static GE_Polygon_TEST* unique_instance ;      
} ;


#endif
