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

#ifndef PEL_ARITHMETIC_EXP_HH
#define PEL_ARITHMETIC_EXP_HH

#include <PEL_Expression.hh>

/*

Common algebraic binary operators : {+,-,*,/,modulo}.
Arguments must be integer of real.

PUBLISHED
*/

class PEL_EXPORT PEL_ArithmeticExp : public PEL_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( PEL_Context const* ct ) const ;

      virtual int to_int( PEL_Context const* ct ) const ;

      virtual std::string const& to_string( PEL_Context const* ct ) const ;

   //-- Formal calculus

      virtual PEL_Data* create_derivative( PEL_Object* a_owner,
                                           PEL_Variable const* var,
                                           PEL_Context const* ct ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
         
   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------

      PEL_ArithmeticExp( void ) ;
     ~PEL_ArithmeticExp( void ) ;
      PEL_ArithmeticExp( PEL_ArithmeticExp const& other ) ;
      PEL_ArithmeticExp& operator=( PEL_ArithmeticExp const& other ) ;

      enum AlgebraicOperator { M,L,T,D,MOD } ;
      
      PEL_ArithmeticExp( PEL_Object* a_owner,
                         std::string const& a_name,
                         PEL_Sequence const* argument_list,
                         AlgebraicOperator a_op ) ;

      PEL_Data const* alternative_result( void ) const ;
      
   //-- Plug in

      PEL_ArithmeticExp( std::string const& a_name,
                         AlgebraicOperator a_op ) ;

      virtual PEL_ArithmeticExp* create_replica( 
                                  PEL_Object * a_owner,
                                  PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments(
                                PEL_Sequence const* some_arguments ) const ;
      
   //-- Formal calculus

      virtual PEL_Data* create_operator_simplification( PEL_Object* a_owner ) ;

   //-- Class attributes
            
      static PEL_ArithmeticExp const* PROTOTYPE_M ;
      static PEL_ArithmeticExp const* PROTOTYPE_L ;
      static PEL_ArithmeticExp const* PROTOTYPE_T ;
      static PEL_ArithmeticExp const* PROTOTYPE_D ;
      static PEL_ArithmeticExp const* PROTOTYPE_MOD ;

   //-- Attributes
      
      AlgebraicOperator const OP ;
      PEL_Data const* const ARG0 ;
      PEL_Data const* const ARG1 ;

      mutable std::string RESULT_STR ;
      
} ;

#endif
