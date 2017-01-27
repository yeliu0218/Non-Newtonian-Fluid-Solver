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

#ifndef PEL_CONTEXT_SIMPLE_HH
#define PEL_CONTEXT_SIMPLE_HH

#include <PEL_Context.hh>

class PEL_Vector ;
class PEL_Variable ;

class PEL_EXPORT PEL_ContextSimple : public PEL_Context
{

   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      // Create, initialize and return an instance.
      static PEL_ContextSimple* create( PEL_Object* a_owner ) ;

      // Provide new context pointing on `self' values and upper context.
      virtual PEL_ContextSimple* create_clone( PEL_Object* a_owner ) const ;

   //-- Access

      // variable number
      virtual size_t nb_variables( void ) const ;

      // `i'th variable
      virtual PEL_Variable const* variable( size_t i ) const ;

      // Is `var' included as a variable ?
      virtual bool has_variable( PEL_Variable const* var ) const ;

      // data identified by `var' (possibly null)
      virtual PEL_Data* value( PEL_Variable const* var ) const ;

   //-- Modifier
      
      // Does the evaluation of `a_value' requires `var' ?
      static bool has_circular_definition( PEL_Variable const* var,
                                           PEL_Data const* a_value ) ;
      
      // Ensure that `var' is included as a variable with associated `a_value'.
      void extend( PEL_Variable const* var, PEL_Data const* a_value ) ;

      // Extend `self' with all variable contained in `other'.
      // All value are cloned to `self', upper context is not took in account.
      void extend( PEL_Context const* other ) ;

      // Does `a_value' the value associated to `var'.
      void set_value_of( PEL_Variable const* var,
                         PEL_Data const* a_value ) ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      PEL_ContextSimple( void ) ;
     ~PEL_ContextSimple( void ) ;
      PEL_ContextSimple( PEL_ContextSimple const& other ) ;
      PEL_ContextSimple& operator=( PEL_ContextSimple const& other ) ;

      PEL_ContextSimple( PEL_Object* a_owner ) ;

   //-- Instance delivery and initialization

      virtual void update( void ) ;

      virtual void update_for_destruction_of( PEL_Context const* subject ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      PEL_Vector* VALUES ;
};

#endif
