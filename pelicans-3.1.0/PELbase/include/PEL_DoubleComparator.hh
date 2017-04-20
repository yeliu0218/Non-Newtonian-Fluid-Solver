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

#ifndef PEL_DOUBLE_COMPARATOR_HH
#define PEL_DOUBLE_COMPARATOR_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

/*
Server comparing double values.

FRAMEWORK INSTANTIATION

   1. Derive a concrete subclass, say MyComparator.
   2. Choose a name for MyComparator, say "my_comparator".
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a default constructor that initializes the
          `PEL_DoubleComparator::' subobject by calling
               `PEL_DoubleComparator( std::string const& )'
          with "my_comparator" as argument.
          Example of pseudo-code :
          | MyComparator:: MyComparator( void )
          |    : PEL_DoubleComparator( "my_comparator" ) {}
      5.2 Define and initialize a static instance by calling the default
          constructor.
             declaration (in the header file, eg MyComparator.hh) :
             | static MyComparator const* PROTOTYPE ;
             definition (in the implementation file, eg MyComparator.cc) :
             | MyComparator const* MyComparator::PROTOTYPE = new MyComparator() ;'
   6. Implement a private constructor that initializes the 
      `PEL_DoubleComparator::' subobject by calling
                 `PEL_DoubleComparator( PEL_Object* )'
      Example of pseudo-code :
      | MyComparator:: MyComparator( PEL_Object* a_owner,
      |                              PEL_ModuleExplorer const* exp )
      |    : PEL_DoubleComparator( a_owner ), ...
      | { ... }
   7. Implement the `::create_replica' method that allocates an object
      of type `MyComparator' initialized using the private constructor described
      above, and subsequently return a pointer to that object.
      Example of pseudo-code :
      | PEL_DoubleComparator const* MyComparator::create_replica(
      |                                   PEL_Object* a_owner,
      |                                   PEL_ModuleExplorer const* exp ) const
      | {
      |    PEL_LABEL( "MyComparator::create_replica" ) ;
      |    PEL_CHECK_PRE( create_replica_PRE( a_owner, exp ) ) ;
      |    MyComparator const* result = new MyComparator( a_owner, exp ) ;
      |    PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ;
      |    return result ;
      | }
   8. Implement the `::three_way_comparison' method

PUBLISHED
*/

class PEL_EXPORT PEL_DoubleComparator : public PEL_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance of `PEL_DoubleComparator::' according
      // to the data attainable by `exp'.
      static PEL_DoubleComparator const* make(
                                           PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp ) ;

   //-- Comparison

      // if `x' equal to `y', 0 ; if smaller, <0 ; if greater, >0  
      virtual int three_way_comparison( double x, double y ) const = 0 ;
      
   protected: //--------------------------------------------------------------
      
      PEL_DoubleComparator( PEL_Object* a_owner ) ;

   //-- Plug in

      virtual ~PEL_DoubleComparator( void ) ;

      PEL_DoubleComparator( std::string const& name ) ;
    
      PEL_DoubleComparator( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp ) ;

      virtual PEL_DoubleComparator const* create_replica(
                            PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp ) const = 0 ;

      bool is_a_prototype( void ) const ;

   //-- Preconditions, Postconditions, Invariant
      
      bool create_replica_PRE( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const ;
      bool create_replica_POST( PEL_DoubleComparator const* result,
				PEL_Object* a_owner,
				PEL_ModuleExplorer const* exp ) const ;
      
      virtual bool invariant( void ) const ;
      
   private: //----------------------------------------------------------------

      PEL_DoubleComparator( void ) ;
      PEL_DoubleComparator( PEL_DoubleComparator const& other ) ;
      PEL_DoubleComparator& operator=( PEL_DoubleComparator const& other ) ;
      
      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes
      
      bool IS_PROTO ;
} ;

#endif
