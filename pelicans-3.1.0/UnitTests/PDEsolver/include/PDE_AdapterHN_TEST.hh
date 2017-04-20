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

#ifndef PDE_ADAPTER_HN_TEST
#define PDE_ADAPTER_HN_TEST

#include <PEL_ObjectTest.hh>

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalFE ;
class PDE_SetOfDiscreteFields ;

class PEL_EXPORT PDE_AdapterHN_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_AdapterHN_TEST( void ) ;
     ~PDE_AdapterHN_TEST( void ) ;
      PDE_AdapterHN_TEST( PDE_AdapterHN_TEST const& other ) ;
      PDE_AdapterHN_TEST& operator=( PDE_AdapterHN_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      void process_domain( size_t ref_step,
                           PDE_DomainAndFields const* dom,
                           std::ostream& os ) ;

      void iterate( size_t ref_step,
                    PDE_SetOfDiscreteFields* sdf,
                    PDE_LocalFE* fe,
                    std::ostream& os ) ;

      void check_node_location( size_t ref_step,
                                PDE_SetOfDiscreteFields* sdf,
                                PDE_LocalFE* fe,
                                bool& ok_node_loc ) const ;

      void iterate_fields( size_t ref_step,
                           PDE_SetOfDiscreteFields* sdf,
                           std::ostream& os ) const ;
      
      void display_error( PDE_LocalFE const* fe,
                          double theo, double xx ) const ;
      
      void print_field( PDE_DiscreteField const* ff,
                        std::ostream& os ) const ;

   //-- Class attributes

      static PDE_AdapterHN_TEST* REGISTRATOR ;

   //-- Attributes

      std::string TEST_NAME ;

      GE_QRprovider const* QRP ;
} ;

#endif
