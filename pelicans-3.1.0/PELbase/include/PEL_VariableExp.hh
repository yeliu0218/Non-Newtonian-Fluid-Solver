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

#ifndef PEL_VARIABLE_EXP_HH
#define PEL_VARIABLE_EXP_HH

#include <PEL_Expression.hh>

/*
Operators on `PEL_Variable::' objects.

- check if a given variable in an HDS

  if( is_defined( "DS_value" ) )
  MODULE toto
     ...
  END MODULE toto

  Module toto is built only if a variable of name $DS_value has been
  defined.

- value of a variable with default value

  my_variable = value( "DS_value", 1.E-5 )

  my_variable is set with the value of the variable $DS_value this variable
  is defined, or with the default value 1.E-5 elsewhere.
  
*/

class PEL_EXPORT PEL_VariableExp : public PEL_Expression
{
   public: //---------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual bool to_bool( PEL_Context const* ct ) const ;
      virtual double to_double( PEL_Context const* ct ) const ;
      virtual int to_int( PEL_Context const* ct ) const ;
      virtual std::string const& to_string( PEL_Context const* ct ) const ;
      virtual doubleVector const& to_double_vector( PEL_Context const* ct ) const ;
      virtual intVector const& to_int_vector( PEL_Context const* ct ) const ;
      virtual stringVector const& to_string_vector(
                                              PEL_Context const* ct ) const ;
      virtual boolVector const& to_bool_vector( PEL_Context const* ct ) const ;
      virtual doubleArray2D const& to_double_array2D( PEL_Context const* ct ) const ;
      virtual intArray2D const& to_int_array2D( PEL_Context const* ct ) const ;
      virtual boolArray2D const& to_bool_array2D( 
                                        PEL_Context const* ct = 0 ) const ;      
      virtual stringArray2D const& to_string_array2D( 
                                        PEL_Context const* ct = 0 ) const ;
 
   //-- Formal calculus

      virtual bool is_constant( void ) const ;

   protected: //-------------------------------------------------------

   private: //---------------------------------------------------------

      PEL_VariableExp( void ) ; 
     ~PEL_VariableExp( void ) ;
      PEL_VariableExp( PEL_VariableExp const& other ) ;
      PEL_VariableExp& operator=( PEL_VariableExp const& other ) ;

      enum VarExp{ var_def, var_value } ;

      PEL_VariableExp( PEL_Object* a_owner,
                       VarExp exp_id,
                       std::string const& a_name,
                       PEL_Sequence const* argument_list ) ;
      
   //-- Plug in

      PEL_VariableExp( VarExp exp_id, std::string const& a_name ) ;
      
      virtual PEL_VariableExp* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;
   
      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static PEL_VariableExp const* PROTOTYPE_var_def ;
      static PEL_VariableExp const* PROTOTYPE_var_value ;

   //-- Attributes

      VarExp const OP ;
} ;

#endif
