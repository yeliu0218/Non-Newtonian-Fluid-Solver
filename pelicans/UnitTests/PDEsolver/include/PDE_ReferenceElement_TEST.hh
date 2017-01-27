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

#ifndef PDE_REFERENCE_ELEMENT_TEST
#define PDE_REFERENCE_ELEMENT_TEST

#include <PEL_ObjectTest.hh>

class PDE_ReferenceElement ;

class PEL_ModuleExplorer ;
class PEL_Randomizer ;

class doubleVector ;

class PEL_EXPORT PDE_ReferenceElement_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_ReferenceElement_TEST( void ) ;
     ~PDE_ReferenceElement_TEST( void ) ;
      PDE_ReferenceElement_TEST( PDE_ReferenceElement_TEST const& other ) ;
      PDE_ReferenceElement_TEST& operator=( 
                                 PDE_ReferenceElement_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* texp ) ;

      void process_geometry( std::string const& name,
                             PDE_ReferenceElement const* ref,
                             PEL_ModuleExplorer const* exp ) ;

      void process_bf_value_at_nodes(
                             PDE_ReferenceElement const* ref,
                             PEL_ModuleExplorer const* exp ) ;

      void process_bf_value_at_random_points(
                             PDE_ReferenceElement const* ref,
                             PEL_ModuleExplorer const* exp ) ;
      
      void process_dbf_value_at_random_points(
                             PDE_ReferenceElement const* ref,
                             PEL_ModuleExplorer const* exp ) ;
      
      void process_d2bf_value_at_random_points(
                             PDE_ReferenceElement const* ref,
                             PEL_ModuleExplorer const* exp ) ;

      void process_value_at_random_points(
                             PDE_ReferenceElement const* ref,
                             PEL_ModuleExplorer const* exp ) ;

      void initialize_randomizer( void ) ;
      void set_random_point( PDE_ReferenceElement const* ref,
                             doubleVector& coords ) ;
      
   //-- Class attributes

      static PDE_ReferenceElement_TEST* REGISTRATOR ;

   //-- Attributes

      PEL_Randomizer* RANDOMIZER ;

} ;
#endif
