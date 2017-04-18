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

#ifndef PEL_EXPRESSION_HH
#define PEL_EXPRESSION_HH

#include <PEL_Data.hh>

class PEL_ContextSimple ;
class PEL_Iterator ;
class PEL_ObjectRegister ;
class PEL_Sequence ;

/*
Context dependant data.

An expression is identified by :
   - its name, and
   - its list of arguments, where each argument is itself a `PEL_Data::'
     object.

IMPLEMENTATION :
   A register is maintained that maps the name of each concrete subclass
   with a prototype of this subclass. That prototype is an instance
   with a NULL argument list (and is thus in an invalid state) that serves
   as a model to create and initialize instances with the `create_replica'
   method.

FRAMEWORK INSTANTIATION :
   1. Derive a concrete subclass, say MyExp (derivation of 
      abstract subclasses is improbable although possible). 
   2. Choose a name for `MyExp', say "my_exp". 
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a contructor (eg the default constructor) that initializes
          the `PEL_Expression::' subobject by calling
                   `PEL_Expression( std::string const& )'
          with "my_exp" as argument.
          Example of pseudo-code :
          | MyExp:: MyExp( void ) : PEL_Expression( "my_exp" ) {}
      5.2 Define and initialize a static instance by calling the default 
          constructor.
             declaration :
             |    static MyExp const* PROTOTYPE ;
             definition :
             |    MyExp const* MyExp::PROTOTYPE = new MyExp() ;
   6. Implement the `create_replica' methods, that calls a constructor.
   7. Implement the constructor called by `create_replica'. The 
      `PEL_Expression' subobject is initialized by calling     
      `PEL_Expression( PEL_Object*, std::string const&, PEL_Sequence const* )'
      with "my_exp" as a second argument.
   8. Implement all necessary method evaluating the type and 
      value of MyExp (eg `to_double', `to_int', `data_type' for an expression
      whose value can be of type int and double).
*/

class PEL_EXPORT PEL_Expression : public PEL_Data
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      // Create, initialize and return an instance.
      static PEL_Expression* create( PEL_Object* a_owner,
                                     std::string const& a_name,
                                     PEL_Sequence const* argument_list,
                                     std::string const& a_comment = "" ) ;
      
      virtual PEL_Expression* create_clone( PEL_Object* a_owner ) const ;

   //-- Identification

      // name
      std::string const& name( void ) const ;
      
   //-- Context
      
      virtual void declare( PEL_List* lst ) const ;

      // IMPLEMENTATION :
      //    `PEL_Data::context_has_required_variables'( `ct' ) for all
      //    the data of `self' ?
      virtual bool context_has_required_variables( 
                                           PEL_Context const* ct ) const ;

   //-- Value

      // IMPLEMENTATION :
      //    `PEL_Data::'value_can_be_evaluated( `ct' ) for all
      //    the data of `self' ?
      virtual bool value_can_be_evaluated( PEL_Context const* ct ) const ;

      // IMPLEMENTATION :
      //    vector of all `PEL_Data::'undefined_variables( `ct' ) for all
      //    the data of `self'
      virtual stringVector const& undefined_variables(
                                           PEL_Context const* ct ) const ;

   //-- Formal calculus

      virtual PEL_Data* create_derivative( PEL_Object* a_owner,
                                           PEL_Variable const* var,
                                           PEL_Context const* ct  ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      static void print_prototypes( std::ostream& os, size_t indent_width ) ;

      bool external_brackets_are_set( void ) const ;
      void set_external_brackets( void ) ;
      void unset_external_brackets( void ) ;
      
   //-- Registered instances (1010.0)
      
      static stringVector const& registered_expressions( void ) ;
      
      static std::string const& usage_of( std::string const& a_name ) ;
      
      static bool valid_arguments_of( std::string const& a_name,
                                      PEL_Sequence const* argument_list ) ;
      
   protected: //-------------------------------------------------------
      
   //-- Plug in

      virtual ~PEL_Expression( void ) ;

      // for prototype registration only
      PEL_Expression( std::string const& a_name ) ;

      // Construction of an instance called `a_name', whose
      // owner is `a_owner' with the items of `argument_list' as arguments.
      PEL_Expression( PEL_Object* a_owner,
                      std::string const& a_name,
                      PEL_Sequence const* argument_list ) ;

      // Is `self' a prototype ?
      bool is_a_prototype( void ) const ;
      
      virtual PEL_Expression* create_replica( 
                              PEL_Object* a_owner,
                              PEL_Sequence const* argument_list ) const  = 0 ;

      
   //-- Characteristics

      // user documentation of `self'
      virtual std::string const& usage( void ) const = 0 ;
      
      // Is the list arguments stored in `some_arguments' valid ? 
      virtual bool valid_arguments(
                              PEL_Sequence const* some_arguments ) const = 0 ;

      // comment defined for `self'
      std::string const& comment( void ) const ;

   //-- Formal calculus

      virtual PEL_Data* create_non_const_simplification(
                                           PEL_Object* a_owner ) const ;
      virtual bool is_raw_data( void ) const ;
     
      // Create operator-specific omptimization if any.
      virtual PEL_Data* create_operator_simplification( PEL_Object* a_owner ) ;

   //-- Arguments

      // number of arguments
      size_t nb_arguments( void ) const ;
      
      // `idx'-th argument
      PEL_Data const* arg( size_t idx ) const ;

      // `idx'-th item of `some_arguments'
      static PEL_Data const* extract_arg( PEL_Sequence const* some_arguments,
                                          size_t idx ) ;

   //-- Error

      void raise_error( std::string const& message ) const ;
      
   //-- Preconditions, Postconditions, Invariant      

      virtual bool create_replica_PRE( 
                              PEL_Object const* a_owner,
                              PEL_Sequence const* argument_list ) const ;

      virtual bool create_replica_POST(
	                      PEL_Expression const* result,
                              PEL_Object const* a_owner,
                              PEL_Sequence const* argument_list ) const ;

      virtual bool valid_arguments_PRE(
                              PEL_Sequence const* some_arguments ) const ;

      virtual bool create_operator_simplification_POST(
                              PEL_Object const* a_owner,
                              PEL_Data const* result ) const ;
      
      virtual bool invariant( void ) const ;
      
   private: //-------------------------------------------------------

      PEL_Expression( void ) ;
      PEL_Expression( PEL_Expression const& other ) ;
      PEL_Expression& operator=( PEL_Expression const& other ) ;

      static stringVector& plugins_name( void ) ;
      static PEL_ObjectRegister* plugins_map( void ) ;
      

   //-- Attributes

      std::string const NAME ;
      bool HAS_BRACKETS ;
      std::string COMMENT ;
      PEL_Sequence const* const ARGUMENTS ; // sequence of PEL_Data*
      mutable PEL_Iterator* IT ;            // iterator on ARGUMENTS
} ;

#endif
