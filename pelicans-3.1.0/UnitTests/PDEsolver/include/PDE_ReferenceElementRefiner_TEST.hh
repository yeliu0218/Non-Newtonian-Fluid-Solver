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

#ifndef PDE_REFERENCE_ELEMENT_REFINER_TEST
#define PDE_REFERENCE_ELEMENT_REFINER_TEST

#include <PEL_ObjectTest.hh>

#include <vector>

class size_t_vector ;

class GE_Mpolyhedron ;
class GE_QuadratureRule ;
class GE_SetOfPoints ;

class PDE_ReferenceElementRefiner ;

class PEL_EXPORT PDE_ReferenceElementRefiner_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_ReferenceElementRefiner_TEST( void ) ;
     ~PDE_ReferenceElementRefiner_TEST( void ) ;
      PDE_ReferenceElementRefiner_TEST( 
                          PDE_ReferenceElementRefiner_TEST const& other ) ;
      PDE_ReferenceElementRefiner_TEST& operator=( 
                          PDE_ReferenceElementRefiner_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      GE_Mpolyhedron const* create_coarse_polyhedron( 
                                             GE_SetOfPoints* a_owner,
                                             PEL_ModuleExplorer const* exp ) ;

      void build_new_geometry( GE_SetOfPoints* verts ) ;

      void check_refi_coefficients( void ) ;

      void check_refi_equation( GE_QuadratureRule const* qr ) ;

      void display_error_coef( size_t ic, 
                               size_t child_n,
                               size_t pn, 
                               double xx_1, double xx_2 ) const ;

      void display_error( size_t pn, double xx_1, double xx_2 ) const ;

   //-- Class attributes

      static PDE_ReferenceElementRefiner_TEST* REGISTRATOR ;

   //-- Attributes

      std::string TEST_NAME ;
      size_t DIM ;
      size_t VERB_LEVEL ;
      PDE_ReferenceElementRefiner const* ELRF ;
      GE_Mpolyhedron const* CPOLY ;
      std::vector< GE_Mpolyhedron const* > RPOLYS ;
      double MY_EPS ;
      double MY_MIN ;
} ;

#endif
