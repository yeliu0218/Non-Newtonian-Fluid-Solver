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

#ifndef PEL_VARIABLE_HH
#define PEL_VARIABLE_HH

#include <PEL_Data.hh>

class PEL_Context ;
class PEL_ContextSimple ;

/*
Variables on which expressions depend.

The result of a variable evaluation in a context is its value
in that context.
Type of variable depends on its name (after $ symbol indicating that
identifier stands for a variable name):
 First character stands for element type ( D for Double, I for
  Integer, B for Boolean, S for String );
 Second one stands for dimension ( S for Scalar, V for Vector and A for
  2D-Array ).
  
Examples :
$SS_HW = "Hello, world!"   // A scalar string
$DVvar = < 0.0 1.0 2.0 >   // A vector of doubles
$BSistrue = true           // A scalar boolean
$IA_itab = array( < 0 0 > , < 1 1 > ) // An integer array
*/

class PEL_EXPORT PEL_Variable : public PEL_Data
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      virtual PEL_Variable* create_clone( PEL_Object* a_owner ) const ;

      // number of instances
      static size_t nb_objects( void ) ;

      // an instance identified by `a_name'
      static PEL_Variable const* object( std::string const& a_name ) ;
      
      static PEL_Variable const* object( size_t id ) ;
      
   //-- Identification

      // number, uniquely determining `self'      
      size_t id_number( void ) const ;

      std::string const& name( void ) const ;

   //-- Context

      // IMPLEMENTATION : ensure that `lst' contains `self'.
      virtual void declare( PEL_List* lst ) const ;

      virtual bool context_has_required_variables( 
                                             PEL_Context const* ct ) const ;

   //-- Type

      // type of data of name `a_name' (raise an error in case of bad syntaxe)
      static PEL_Data::Type data_type( std::string const& a_name ) ;
      
      virtual PEL_Data::Type data_type( void ) const ;
      
  //-- Value
      
      // IMPLEMENTATION : true if `ct' associates a value to `self'
      virtual bool value_can_be_evaluated( PEL_Context const* ct ) const ;
      
      virtual stringVector const& undefined_variables(
                                           PEL_Context const* ct ) const ;

      virtual bool to_bool( PEL_Context const* ct ) const ;

      virtual double to_double( PEL_Context const* ct ) const ;

      virtual int to_int( PEL_Context const* ct ) const ;

      virtual std::string const& to_string( PEL_Context const* ct ) const ;

      virtual doubleVector const& to_double_vector(
                                            PEL_Context const* ct ) const ;

      virtual intVector const& to_int_vector( PEL_Context const* ct ) const ;

      virtual stringVector const& to_string_vector( 
                                            PEL_Context const* ct ) const ;

      virtual boolVector const& to_bool_vector( 
                                            PEL_Context const* ct ) const ;

      virtual doubleArray2D const& to_double_array2D( 
                                            PEL_Context const* ct ) const ;

      virtual intArray2D const& to_int_array2D( PEL_Context const* ct ) const ;
      
      virtual boolArray2D const& to_bool_array2D( 
                                        PEL_Context const* ct = 0 ) const ;
      
      virtual stringArray2D const& to_string_array2D( 
                                        PEL_Context const* ct = 0 ) const ;
      
      virtual doubleArray3D const& to_double_array3D( 
                                            PEL_Context const* ct ) const ;

      virtual intArray3D const& to_int_array3D( PEL_Context const* ct ) const ;
      
   //-- Formal calculus

      virtual PEL_Data* create_derivative( PEL_Object* a_owner,
                                           PEL_Variable const* var,
                                           PEL_Context const* ct ) const ;
      virtual bool is_raw_data( void ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      PEL_Variable( void ) ;
     ~PEL_Variable( void ) ;
      PEL_Variable( PEL_Variable const& other ) ;
      PEL_Variable& operator=( PEL_Variable const& other ) ;

      PEL_Variable( PEL_Object* a_owner, std::string const& a_name ) ;

      PEL_Variable( PEL_Object* a_owner, PEL_Variable const* other ) ;

      // value of `self' in `ct'.
      PEL_Data const* data( PEL_Context const* ct ) const ;
      
      static PEL_List*  variable_list( void ) ;

   //-- Attributes

      std::string const NAME ;
      PEL_Data::Type const KIND ;
      size_t ID ;

      mutable bool EVALUATING ;
};

#endif
