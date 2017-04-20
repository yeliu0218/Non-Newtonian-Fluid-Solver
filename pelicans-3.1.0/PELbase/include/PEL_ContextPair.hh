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

#ifndef PEL_CONTEXT_PAIR_HH
#define PEL_CONTEXT_PAIR_HH

#include <PEL_Context.hh>

class PEL_Data ;
class PEL_Vector ;
class PEL_Variable ;

class PEL_EXPORT PEL_ContextPair : public PEL_Context
{

   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      // Create, initialize and return an instance (the result is the merge
      // of `first' and `second' context, `second' having priority).
      static PEL_ContextPair* create( PEL_Object* a_owner,
                                      PEL_Context const* first,
                                      PEL_Context const* second) ;
      
      // Reinitialize the internal state, as if the `::create' method was
      // just completed (`self' is the merge of `first' and `second' context,
      // `second' having priority).
      void re_initialize( PEL_Context const* first,
                          PEL_Context const* second ) ;

      // Provide new context pointing on `self' values and upper context.
      virtual PEL_ContextPair* create_clone( PEL_Object* a_owner ) const ;

   //-- Access

      // variable number
      virtual size_t nb_variables( void ) const ;

      // `i'-th variable
      virtual PEL_Variable const* variable( size_t i ) const ;

      // Is `var' included as a variable ?
      virtual bool has_variable( PEL_Variable const* var ) const ;

      // data identified by `var' (possibly null)
      virtual PEL_Data* value( PEL_Variable const* var ) const ;

   protected: //-------------------------------------------------------

   private: //---------------------------------------------------------

      PEL_ContextPair( void ) ;
     ~PEL_ContextPair( void ) ;
      PEL_ContextPair( PEL_ContextPair const& other ) ;
      PEL_ContextPair& operator=( PEL_ContextPair const& other ) ;

      PEL_ContextPair( PEL_Object* a_owner,
                       PEL_Context const* first,
                       PEL_Context const* second ) ;

   //-- Instance delivery and initialization

      virtual void update( void ) ;

      virtual void update_for_destruction_of( PEL_Context const* subject ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes
      
      PEL_Context const* CT1 ;
      PEL_Context const* CT2 ;
      PEL_Vector* VALUES ;
};

#endif
