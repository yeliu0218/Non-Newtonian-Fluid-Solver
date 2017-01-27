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

#ifndef PEL_DATA_WITH_CONTEXT_EXP_HH
#define PEL_DATA_WITH_CONTEXT_EXP_HH

#include <PEL_TransferExp.hh>

/*

Expressions associated to a set of variables.

---
name     : data_with_context
argument : expression, [ variable name (String), variable value ]
type     : first argument type

The returned value is the evaluation of the expression associated to a
context defined with a set of pairs (variable name, variable value).

Example:
  
  MODULE titi
     $DS_X = 3.
     toto = data_with_context( 3.*$DS_X*$DS_ALPHA, "DS_ALPHA", 2. )
  END MODULE titi

PUBLISHED
*/

class PEL_EXPORT PEL_DataWithContextExp : public PEL_TransferExp
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      static PEL_DataWithContextExp* create(
                      PEL_Object* a_owner,
                      PEL_Data const* data, PEL_Context const* ct ) ;
      
   //-- Context
      
      virtual void declare( PEL_List* lst ) const ;
      
      virtual bool context_has_required_variables(
                                      PEL_Context const* ct ) const ;

      
   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual bool value_can_be_evaluated( PEL_Context const* ct ) const ;

      virtual stringVector const& undefined_variables(
                                           PEL_Context const* ct ) const ;
 
   //-- Formal calculus

      virtual bool is_constant( void ) const ;

   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------
      
      PEL_DataWithContextExp( void ) ;
     ~PEL_DataWithContextExp( void ) ;
      PEL_DataWithContextExp( PEL_DataWithContextExp const& other ) ;
      PEL_DataWithContextExp& operator=(
                              PEL_DataWithContextExp const& other ) ;
      
      PEL_DataWithContextExp( PEL_Object* a_owner,
                              std::string const& a_name,
                              PEL_Sequence const* argument_list ) ;

   //-- Plug in

      PEL_DataWithContextExp( std::string const& a_name ) ;

      virtual PEL_DataWithContextExp* create_replica( 
                      PEL_Object* a_owner,
                      PEL_Sequence const* argument_list ) const ;

      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Transfer implementation
      
      virtual PEL_Data const* data( PEL_Context const* ct ) const ;

   //-- Class attributes

      static PEL_DataWithContextExp const* PROTO ;
      
   //-- Attributes
      
      PEL_Data const* DATA ;
} ;

#endif
