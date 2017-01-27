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

#ifndef GE_M_POLYHEDRON_TEST_HH
#define GE_M_POLYHEDRON_TEST_HH

#include <PEL_ObjectTest.hh>

#include <GE_Matrix.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <doubleVector.hh>

class GE_Mpolyhedron ;

/* Unit test for the GE_Mpolyhedron class.

It is essentially coumpound on a test that checks the general properties
of the polyhedron ( name, measure, ).*/


class PEL_EXPORT GE_Mpolyhedron_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      GE_Mpolyhedron_TEST( void ) ;
     ~GE_Mpolyhedron_TEST( void ) ;
      GE_Mpolyhedron_TEST( GE_Mpolyhedron_TEST const& other ) ;
      GE_Mpolyhedron_TEST& operator=( GE_Mpolyhedron_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;
      
      GE_Mpolyhedron const* create_polyhedron( PEL_Object* a_owner,
                                               PEL_ModuleExplorer const* exp ) ;

      void adapt_to_polyhedron( GE_Mpolyhedron const* poly ) ;

      void test_generalities( PEL_ModuleExplorer const* exp,
                              GE_Mpolyhedron const* poly ) ;
      
      void test_mapping( GE_Mpolyhedron const* poly ) ;
      
      void test_mapping_derivative( GE_Mpolyhedron const* poly,
                                    double dx, double d_eps, double d_min ) ;

      void test_contains( GE_Mpolyhedron const* poly ) ;

      void test_if_fv_center_is_circumcenter( GE_Mpolyhedron const* poly ) ;

      double estimate_measure( GE_Mpolyhedron const* poly ) const ;

      double compute_inter_vert_max_dist( GE_Mpolyhedron const* poly ) const ;

      double compute_equi_ball_diameter( GE_Mpolyhedron const* poly ) const ;

      double compute_inter_vert_max_dist( GE_Mpolyhedron const* poly,
                                          size_t dim ) const ;

      //-- Class attributes

      static GE_Mpolyhedron_TEST* REGISTRATOR ;

      //-- Attributes

      GE_Point* PT_REF ;
      GE_Point const* ORIGIN ;
      GE_Point* PT ;
      GE_Point* PT2 ;

      GE_Vector* VECTOR ;

      GE_Matrix* JACOB ;
} ;


#endif
