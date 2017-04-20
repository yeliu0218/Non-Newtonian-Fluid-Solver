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

#ifndef PEL_DATA_HH
#define PEL_DATA_HH

#include <PEL_Object.hh>

class boolVector ;
class boolArray2D ;
class doubleArray2D ;
class doubleArray3D ;
class doubleVector ;
class intArray2D ;
class intArray3D ;
class intVector ;
class stringVector ;
class stringArray2D ;
class PEL_Context ;
class PEL_List ;
class PEL_ContextSimple ;
class PEL_Variable ;

/*
Data of the PELICANS Hierarchical Data System.

A Data has two characteristics : 
   - a type, and
   - a value.

The value of a data might depend of a context (but not its type).

A context is  represented by a `PEL_Context::' object (possibly NULL in
which case the type and the value are context independant). 

To evaluate a data in a given context (resulting in its type and value) :
   1. the required variables should occur in the context (which can be
      achieved by calling `declare' ;
   2. the context should assign relevant values to these required variables.
*/

class PEL_EXPORT PEL_Data : public PEL_Object
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      virtual PEL_Data* create_clone( PEL_Object* a_owner ) const = 0 ;
      
   //-- Context(10.)
                  
      // Ensure that `lst' contains all necessary instances of
      // `PEL_Variable::' for the evaluation of `self'.
      // IMPLEMENTATION : do nothing : `lst' is left unchanged.
      virtual void declare( PEL_List* lst ) const ;

      // Does `ct' contains all necessary variables for the evaluation 
      // of `self' (the variables of `ct' can not necessary be evaluated) ?
      // IMPLEMENTATION : true
      virtual bool context_has_required_variables( 
                                               PEL_Context const* ct ) const ;
      
   //-- Type(20.)

      enum Type { Undefined, Double, Int, Bool, String,
                  DoubleVector, IntVector, BoolVector, StringVector,
                  DoubleArray2D, IntArray2D, BoolArray2D, StringArray2D,
                  DoubleArray3D, IntArray3D } ;

      static std::string type_name( Type kind ) ;
      
      // type of self
      virtual Type data_type( void ) const = 0 ;
      
  //-- Value(30.)
      
      // Is it possible to evaluate the value according to `ct' ?
      // IMPLEMENTATION : true
      virtual bool value_can_be_evaluated( PEL_Context const* ct ) const ;

      // undefined variable names for the evaluation of the value according to `ct'
      // IMPLEMENTATION : empty vector
      virtual stringVector const& undefined_variables(
                                           PEL_Context const* ct ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual bool to_bool( PEL_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual double to_double( PEL_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual int to_int( PEL_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual std::string const& to_string( PEL_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual doubleVector const& to_double_vector(
                                       PEL_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual intVector const& to_int_vector( 
                                       PEL_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual stringVector const& to_string_vector(
                                        PEL_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual boolVector const& to_bool_vector(
                                        PEL_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual doubleArray2D const& to_double_array2D(
                                        PEL_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual intArray2D const& to_int_array2D( 
                                        PEL_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual boolArray2D const& to_bool_array2D( 
                                        PEL_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual stringArray2D const& to_string_array2D( 
                                        PEL_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual doubleArray3D const& to_double_array3D(
                                        PEL_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual intArray3D const& to_int_array3D( 
                                        PEL_Context const* ct = 0 ) const ;

      // value represented as a string evaluated according to `ct'
      virtual std::string value_as_string( 
                                        PEL_Context const* ct = 0 ) const ;

   //-- Formal calculus

      // Return partial derivative of `self' with respect to `var'
      // expanding expressions from `ct'.
      virtual PEL_Data* create_derivative( PEL_Object* a_owner,
                                           PEL_Variable const* var,
                                           PEL_Context const* ct  ) const ;

      // Create a simplified expression of `self'.
      PEL_Data* create_simplification( PEL_Object* a_owner ) const ;
      
      // Is `self' a constant value ?
      virtual bool is_constant( void ) const ;
      
      // Is `self' a raw data ?
      // IMPLEMENTATION : returns true
      virtual bool is_raw_data( void ) const ;
      
   protected: //-------------------------------------------------------
      
      virtual ~PEL_Data( void ) ;

      PEL_Data( PEL_Object* a_owner ) ;

   //-- Formal calculus

      // Called by `::create_simplification'.
      virtual PEL_Data* create_non_const_simplification(
                                               PEL_Object* a_owner ) const ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool declare_PRE( PEL_List const* lst ) const ;
      virtual bool declare_POST( PEL_List const* lst ) const ;
      virtual bool context_has_required_variables_PRE( 
                                             PEL_Context const* ct ) const ;
      virtual bool to_double_PRE( PEL_Context const* ct ) const ;
      virtual bool to_int_PRE( PEL_Context const* ct ) const ;
      virtual bool to_bool_PRE( PEL_Context const* ct ) const ;
      virtual bool to_string_PRE( PEL_Context const* ct ) const ;
      virtual bool to_double_vector_PRE( PEL_Context const* ct ) const ;
      virtual bool to_int_vector_PRE( PEL_Context const* ct ) const ;
      virtual bool to_string_vector_PRE( PEL_Context const* ct ) const ;
      virtual bool to_bool_vector_PRE( PEL_Context const* ct ) const ;
      virtual bool to_double_array2D_PRE( PEL_Context const* ct ) const ;
      virtual bool to_int_array2D_PRE( PEL_Context const* ct ) const ;
      virtual bool to_bool_array2D_PRE( PEL_Context const* ct ) const ;
      virtual bool to_string_array2D_PRE( PEL_Context const* ct ) const ;
      virtual bool to_double_array3D_PRE( PEL_Context const* ct ) const ;
      virtual bool to_int_array3D_PRE( PEL_Context const* ct ) const ;
      virtual bool create_derivative_PRE( PEL_Object* a_owner,
                                          PEL_Variable const* var,
                                          PEL_Context const* ct  ) const ;
      virtual bool create_derivative_POST( PEL_Object* a_owner,
                                           PEL_Variable const* var,
                                           PEL_Data const* result ) const ;
      virtual bool create_non_const_simplification_PRE(
                                           PEL_Object* a_owner ) const ;
      virtual bool create_non_const_simplification_POST(
                                           PEL_Object* a_owner,
                                           PEL_Data const* result ) const ;
     
      virtual bool invariant( void ) const ;
      virtual bool is_raw_data_POST( bool result ) const ;
      
   private: //-------------------------------------------------------

      PEL_Data( void ) ;
      PEL_Data( PEL_Data const& other ) ;
      PEL_Data& operator=( PEL_Data const& other ) ;

      void exitWithError( std::string const& mess ) const ;
};

#endif
