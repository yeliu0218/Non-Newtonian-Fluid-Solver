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

#ifndef PEL_VECTOR_EXP_HH
#define PEL_VECTOR_EXP_HH

#include <PEL_Expression.hh>

#include <doubleVector.hh>
#include <intVector.hh>
#include <boolVector.hh>
#include <stringVector.hh>

class PEL_Iterator ;

/*
Operators to form vector from simple items:
   - vector of doubles:
       vector( 0., 3., 4., -1.E3 )              ->  < 0. 3. 4. -1.E3 >
   - vector of integers:
       vector( 0, 3, 4, -5 )                    ->  < 0 3 4 -5 >
   - vector of strings:
       vector( "PELICANS", "is", "beautiful" )  ->  < "PELICANS" "is" "beautiful" >
   - vector of booleans:
       vector( true, false )                    ->  < true, false >

   - nvector( 2, 1.0 )                          ->  < 1.0 1.0 >
   
   - nvector( 3, true )                         ->  < true true true >
   
Operator to retrieve the size of such vectors:
   size( vector( 0., 3., 4., -1.E3 ) )          ->  4

Operator to retrieve component from vector:
   First arg is a vector and second one an integer.
   component( vector( "PELICANS", "is", "beautiful") ), 1 ) -> "PELICANS"

Operator to test if the elements of a vector are in increasing order
   increasing( < 0.0 1.0 1.0 2.0 > )  -> true
   increasing( < 0 1 2 1 > )          -> false

Operator to test if all the elements of a vector are greater than or equal to
a given value
   greater( < -10 3 >, 1 )      -> false
   greater( < 1.0 3.0 >, 0.5 )  -> true

Operator to form vector from another one:
   apply( < 1. 4. 9. 16. >, $DS_x*$DS_x, "DS_x" ) -> < 1. 4. 9. 16. >
   apply( < 1  4  9  16  >, $IS_x*$IS_x, "IS_x" ) -> < 1  4  9  16  >
   apply( < true true false >, ! $BS_x, "BS_x" )  -> < false false true >
   apply( < "titi" "toto" >, $SS_x+"0", "SS_x" )  -> < "titi0" "toto0" >
   apply( < "titi" "toto" >, $SS_x+to_string($IS_ic), "SS_x", "SS_ic" )
                                                  -> < "titi0" "toto1" >

Operator to sum the elements of a vector :
   sum( < 0. 1.0 2.0 > )  ->  3.0
   
Operator to reverse order of vector :
   reverse( < 0. 1.0 2.0 > )  ->  < 2.0 1.0 0. >
   
PUBLISHED
*/

class PEL_EXPORT PEL_VectorExp : public PEL_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
 
   //-- Context
                  
      virtual void declare( PEL_List* lst ) const ;

      virtual bool context_has_required_variables( 
                                           PEL_Context const* ct ) const ;
      
   //-- Value
      
      virtual bool value_can_be_evaluated( PEL_Context const* ct ) const ;
      
      virtual stringVector const& undefined_variables(
                                           PEL_Context const* ct ) const ;
      
      virtual bool to_bool( PEL_Context const* ct = 0 ) const ;

      virtual double to_double( PEL_Context const* ct = 0 ) const ;

      virtual int to_int( PEL_Context const* ct = 0 ) const ;

      virtual std::string const& to_string( PEL_Context const* ct = 0 ) const ;
      
      virtual doubleVector const& to_double_vector(
                                       PEL_Context const* ct = 0 ) const ;
      
      virtual intVector const& to_int_vector( 
                                       PEL_Context const* ct = 0 ) const ;

      virtual stringVector const& to_string_vector( 
                                       PEL_Context const* ct = 0 ) const ;
      
      virtual boolVector const& to_bool_vector(
                                       PEL_Context const* ct = 0 ) const ;
      
   //-- Formal calculus
      
      virtual PEL_Data* create_derivative( PEL_Object* a_owner,
                                           PEL_Variable const* var,
                                           PEL_Context const* ct ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      PEL_VectorExp( void ) ;
     ~PEL_VectorExp( void ) ;
      PEL_VectorExp( PEL_VectorExp const& other ) ;
      PEL_VectorExp& operator=( PEL_VectorExp const& other ) ;

      enum VectorExp{ vect, size, component,
                      increasing, greater,
                      nvector, cond_vect,
                      apply, reverse, sum } ;

      PEL_VectorExp( PEL_Object* a_owner,
                     VectorExp exp_id,
                     std::string const& a_name,
                     PEL_Sequence const* argument_list ) ;

      PEL_Context const* create_apply_context(
                     PEL_Object* a_owner,
                     PEL_Context const* ctx,
                     PEL_Variable const* v,
                     PEL_Variable const* vic ) const ;
                     
   //-- Plug in

      PEL_VectorExp( VectorExp exp_id, std::string const& a_name ) ;
      
      virtual PEL_VectorExp* create_replica( 
                                    PEL_Object * a_owner,
                                    PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes      

      static PEL_VectorExp const* PROTOTYPE_CONDITIONAL_VECTOR ;
      static PEL_VectorExp const* PROTOTYPE_VECTOR ;
      static PEL_VectorExp const* PROTOTYPE_NVECTOR ;
      static PEL_VectorExp const* PROTOTYPE_SIZE ;
      static PEL_VectorExp const* PROTOTYPE_COMPO ;
      static PEL_VectorExp const* PROTOTYPE_INCREASING ;
      static PEL_VectorExp const* PROTOTYPE_GREATER ;
      static PEL_VectorExp const* PROTOTYPE_APPLY ;
      static PEL_VectorExp const* PROTOTYPE_REVERSE ;
      static PEL_VectorExp const* PROTOTYPE_SUM ;

   //-- Attributes

      VectorExp const OP ;
      
      mutable PEL_Iterator* MY_IT ; // To speed up evaluation
      
      mutable doubleVector RESULT_D ;
      mutable intVector    RESULT_I ;
      mutable boolVector   RESULT_B ;
      mutable stringVector RESULT_S ;
      
} ;

#endif
