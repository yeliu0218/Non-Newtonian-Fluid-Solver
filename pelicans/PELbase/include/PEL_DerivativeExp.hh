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

#ifndef PEL_DERIVATIVE_EXP_HH
#define PEL_DERIVATIVE_EXP_HH

#include <PEL_Expression.hh>

/* Differential operators d , dnum.

PUBLISHED
*/

class PEL_EXPORT PEL_DerivativeExp : public PEL_Expression
{
   public: //-------------------------------------------------------
      
   //-- Instance delivery and initialization

      virtual PEL_DerivativeExp* create_clone( PEL_Object* a_owner ) const ;

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( PEL_Context const* ct = 0 ) const ;

      virtual doubleVector const& to_double_vector( 
                                PEL_Context const* ct = 0 ) const ;
      
   //-- Formal calculus

      virtual PEL_Data* create_derivative( PEL_Object* a_owner,
                                           PEL_Variable const* var,
                                           PEL_Context const* ct ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      PEL_DerivativeExp( void ) ;
     ~PEL_DerivativeExp( void ) ;
      PEL_DerivativeExp( PEL_DerivativeExp const& other ) ;
      PEL_DerivativeExp& operator=( PEL_DerivativeExp const& other ) ;
      
   //-- Plug in

      enum OP_TYPE { d, dnum } ;
      
      PEL_DerivativeExp( std::string const& a_name, OP_TYPE a_op ) ;

      virtual PEL_DerivativeExp* create_replica(
                                   PEL_Object* a_owner,
                                   PEL_Sequence const* argument_list ) const ;

      PEL_DerivativeExp( PEL_Object* a_owner,
                         std::string const& a_name,
                         OP_TYPE a_op,
                         PEL_Sequence const* argument_list ) ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Formal calculus

      virtual PEL_Data* create_non_const_simplification( 
                                                 PEL_Object* a_owner ) const ;

      PEL_Data const* derivative( PEL_Context const* ct ) const ;
      
   //-- Class attributes
      
      static PEL_DerivativeExp const* PROTOTYPE_d ;
      static PEL_DerivativeExp const* PROTOTYPE_dnum ;
      
   //-- Attribute

      OP_TYPE const OP ;
      PEL_Data const* const EXP ;
      PEL_Variable const* const VAR ;
      mutable PEL_Data const* DERIVATIVE ;
} ;

#endif
