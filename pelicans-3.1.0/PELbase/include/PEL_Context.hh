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

#ifndef PEL_CONTEXT_HH
#define PEL_CONTEXT_HH

#include <PEL_Object.hh>

class PEL_Data ;
class PEL_List ;
class PEL_Variable ;

/*
Contexts in which to interpret expressions.

Contexts are mappings from variables to values.
   - variables are non NULL `PEL_Variables::' objects.
   - values that are `PEL_Data::' objects that are possibly NULL (in which
     case the value is said to be cleared).
     
One context can be said to inherit (`::inherit') from another one : it
means that it shares variables from father context with its owns.
Moreover, one context can be extended (`::extend') from another one : it
means that all variables from source one will be copied in new one.
*/

class PEL_EXPORT PEL_Context : public PEL_Object
{

   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      virtual PEL_Context* create_clone( PEL_Object* a_owner ) const = 0 ;

   //-- Access

      // variable number
      virtual size_t nb_variables( void ) const = 0 ;

      // `i'th variable
      virtual PEL_Variable const* variable( size_t i ) const = 0 ;

      // Is `var' included as a variable ?
      virtual bool has_variable( PEL_Variable const* var ) const = 0 ;

      // data identified by `var'
      virtual PEL_Data* value( PEL_Variable const* var ) const = 0 ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      void print( std::ostream& os, size_t indent_width,
                  PEL_Context const* ctx ) const ;
      
  //-- Observer pattern : subject
      
      void attach_observer( PEL_Context* observer ) const ;

      void detach_observer( PEL_Context* observer ) const ;

   protected: //-------------------------------------------------------

      virtual ~PEL_Context( void ) ;

      PEL_Context( PEL_Object* a_owner ) ;

   //-- Observer pattern : subject

      void update_observers( void ) const ;

      void notify_observers_of_my_destruction( void ) const ;

   //-- Instance delivery and initialization

      virtual void update( void ) = 0 ;

      virtual void update_for_destruction_of( 
                                            PEL_Context const* subject ) = 0 ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool variable_PRE( size_t i  ) const ;
      virtual bool variable_POST( PEL_Variable const* result ) const ;
      virtual bool has_variable_PRE( PEL_Variable const* var ) const ;
      virtual bool value_PRE( PEL_Variable const* var ) const ;
      virtual bool value_POST( PEL_Data* result,
                               PEL_Variable const* var ) const ;
      
   private: //---------------------------------------------------------

      PEL_Context( void ) ;
      PEL_Context( PEL_Context const& other ) ;
      PEL_Context& operator=( PEL_Context const& other ) ;

   //-- Attributes

      PEL_List* OBSERVERS ;
};

#endif
