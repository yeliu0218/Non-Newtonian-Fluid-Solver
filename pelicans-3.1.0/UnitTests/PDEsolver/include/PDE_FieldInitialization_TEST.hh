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

#ifndef PDE_FIELD_INITIALIZATION_TEST
#define PDE_FIELD_INITIALIZATION_TEST

#include <PEL_ObjectTest.hh>

#include <stringVector.hh>

class PEL_Data ;
class PEL_ModuleExplorer ;
class PEL_Vector ;
class PEL_ContextSimple ;
class PEL_DoubleVector ;
class doubleVector ;

class GE_Point ;

class PDE_DiscreteField ;
class PDE_LocalFE ;
class PDE_LocalFEbound ;

class PEL_EXPORT PDE_FieldInitialization_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_FieldInitialization_TEST( void ) ;
     ~PDE_FieldInitialization_TEST( void ) ;
      PDE_FieldInitialization_TEST( 
                               PDE_FieldInitialization_TEST const& other ) ;
      PDE_FieldInitialization_TEST& operator=( 
                               PDE_FieldInitialization_TEST const& other ) ;

  //-- Unit test management
      
      virtual void process_one_test( PEL_ModuleExplorer const* texp ) ;
      

      void process_one_test( PEL_ModuleExplorer* exp ) ;

      bool check_from_meshes( PDE_LocalFE* fe,
                              PDE_DiscreteField const* ff,
                              PEL_Data const* val  ) ;

      void fill_context_with_coordinates( GE_Point const* pt ) ;
      
      void display_error( GE_Point const* pt,
                          size_t ic,
                          double xx_data_deck,
                          double xx_calculated ) const ;

   //-- Class attributes

      static PDE_FieldInitialization_TEST const* REGISTRATOR ;

   //-- Attributes

      PEL_Vector* FES ;
      PEL_Vector* I_FIELDS ;
      PEL_Vector* I_VALUES ;
      PEL_Vector* JACOBIANS ;
      PEL_Vector* BD_FIELDS ;
      PEL_Vector* BD_VALUES ;
      PEL_ContextSimple* CONTEXT ;
      PEL_DoubleVector* COORDS ;
} ;
#endif
