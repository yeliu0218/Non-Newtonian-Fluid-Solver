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

#ifndef PDE_COARSENING_TEST
#define PDE_COARSENING_TEST

#include <PEL_ObjectTest.hh>

#include <iosfwd>

class PEL_Vector ;

class LA_Matrix ;

class GE_QRprovider ;

class PDE_AdapterCHARMS ;
class PDE_BlockAssembledSystem ;
class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalFEcell ;
class PDE_LocalEquation ;
class PDE_SystemNumbering ;

class PEL_EXPORT PDE_Coarsening_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_Coarsening_TEST( void ) ;
     ~PDE_Coarsening_TEST( void ) ;
      PDE_Coarsening_TEST( PDE_Coarsening_TEST const& other ) ;
      PDE_Coarsening_TEST& operator=( 
                           PDE_Coarsening_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      void check_coarsening( PDE_SystemNumbering* nmb,
                             LA_Matrix const* coarsest_mat ) ;

      LA_Matrix* create_assembled_matrix( 
                             PEL_Object* a_owner,
                             PDE_SystemNumbering* nmb ) ;

      void compare_matrices( std::string const& fname,
                             LA_Matrix const* m1,
                             LA_Matrix const* m2 ) ;

      void  display_error( size_t row1, size_t col1, double xx_1, 
                           size_t row2, size_t col2, double xx_2 ) const ;

   //-- Class attributes

      static PDE_Coarsening_TEST* REGISTRATOR ;

   //-- Attributes

      std::string TEST_NAME ;
      double MY_EPS ;
      double MY_MIN ;

      PDE_DomainAndFields const* DOM ;
      PDE_AdapterCHARMS* DA ;

      PDE_LocalEquation* ELEMENT_EQ ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;

      LA_Matrix const* MAT_PROTO ;
} ;

#endif
