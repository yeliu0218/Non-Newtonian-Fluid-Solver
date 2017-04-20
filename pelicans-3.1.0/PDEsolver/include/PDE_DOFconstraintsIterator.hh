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

#ifndef PDE_DOF_CONSTRAINTS_ITERATOR_HH
#define PDE_DOF_CONSTRAINTS_ITERATOR_HH

#include <PEL_Object.hh>

#include <vector>

#include <PDE_DiscreteField.hh>
#include <PDE_DOFconstraints.hh>

class PEL_EXPORT PDE_DOFconstraintsIterator : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Characteristics
      
      PDE_DiscreteField const* field( void ) const ;

   //-- Constraints
      
      void start( size_t n, size_t ic ) ;
      
      bool is_valid( void ) const ;
      
      void go_next( void ) ;

      size_t node_of_constraining_DOF( void ) const ;
      
      size_t component_of_constraining_DOF( void ) const ;
      
      double constraint_coefficient( void ) const ;
                  
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_DOFconstraintsIterator( void ) ;
     ~PDE_DOFconstraintsIterator( void ) ;
      PDE_DOFconstraintsIterator( PDE_DOFconstraintsIterator const& other ) ;
      PDE_DOFconstraintsIterator& operator=( 
                                  PDE_DOFconstraintsIterator const& other ) ;

      PDE_DOFconstraintsIterator( PEL_Object* a_owner,
                                  PDE_DiscreteField const* a_ff,
                                  PDE_DOFconstraints const* a_cstr ) ;
      
   //-- Friends
      
      friend PDE_DOFconstraintsIterator* 
         PDE_DiscreteField::create_constraints_iterator( 
                                               PEL_Object* a_owner ) const ;

   //-- Attributes

      PDE_DiscreteField const* FF ;
      PDE_DOFconstraints const* CSTR ;
      std::vector< PDE_DOFconstraints::ConstraintElement >::const_iterator IT ;
      std::vector< PDE_DOFconstraints::ConstraintElement > const* VEC_OF_IT ;
      bool OK_IT ;
} ;


#endif

