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

#ifndef PEL_OBJECT_HH
#define PEL_OBJECT_HH

#include <PEL_export.hh>

#include <iosfwd>
#include <string>

class PEL_ListIdentity ;
class PEL_Module ;
class PEL_ObjectReader ;
class PEL_ObjectWriter ;

/*
Objects of dynamic storage duration,
   - that are referred to, accessed and manipulated exclusively 
     through pointers ;
   - that may be compared according to a total order relation ;
   - that may be hashed into a integer index, for use as keys in hash tables ;
   - whose lifetime is managed with the Ownership Method.
Any developper-written class publicly inherit from PEL_Object 

Each sub-object of type PEL_Object has one owner (possibly NULL).
The owner is set at creation and cannot be modified unless it is NULL. 
If the owner is NULL, the complete object of self is terminated, by
calling explicitely the destroy() method.
Otherwise, the complete object of self is terminated when the owner 
itself is terminated. 
When the complete object of self is terminated, all objects for which
self is the owner are also terminated. 

FRAMEWORK INSTANTIATION
   1. Derive a subclass.
   2. If concrete, implement a private destructor. 
      If abstract, implement a virtual protected destructor.
   3. If concrete, declare all constructors private.
      If abstract, declare the implemented constructors protected and
      declare all other constructors private.
   4. All implemented constructors should initialize the `PEL_Object'
      subobject by calling
         `PEL_Object( PEL_Object*)'
   5. If concrete, implement one or more static methods returning instances.
*/


class PEL_EXPORT PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization(0.1)

      // Create and a clone of `self'.
      // IMPLEMENTATION : raise a fatal error.
      virtual PEL_Object* create_clone( PEL_Object* a_owner ) const ;

      // Perform updating to maintain consistency between `self' and a 
      // set of related objects (observer pattern).
      // IMPLEMENTATION : raise a fatal error.
      virtual void update( void ) ;
      
   //-- Termination(0.2)

      // Terminate `self'. (Pseudo-destructor that has access to the real 
      // destructor : a call to destroy leads to the destruction of `self').
      void destroy( void ) const ;

      // Terminate `a_possession'.
      void destroy_possession( PEL_Object const* a_possession ) ;

   //-- Identification(0.3)

      // address
      size_t address( void ) const ;

   //-- Characteristics(0.5)

      // name of the class of which `self' is an instance
      std::string const& type_name( void ) const ;
            
      // owner (possibly NULL)
      PEL_Object const* owner( void ) const ;

      // Is `other' in the tree of the owners of `self' ?
      bool is_under_ownership_of( PEL_Object const* other ) const ;

      // Is type of `self' identical to type of `other' ?
      bool same_type( PEL_Object const* other ) const ;

   //-- Characteristics setting(0.6)
      
      // Make `a_owner' the owner of self.
      void set_owner( PEL_Object* a_owner ) ;

      // Make `a_owner' the owner of some possession `a_possession'.
      void change_owner( PEL_Object* a_owner, PEL_Object* a_possession ) ;

   //-- Comparison(0.7)

      // Is `other' comparable to `self' ? 
      // IMPLEMENTATION : `::same_type(other)'
      virtual bool comparable( PEL_Object const* other ) const ;

      // Is `other' equal to `self' ? 
      // IMPLEMENTATION : `::has_same_address(other)'
      virtual bool is_equal( PEL_Object const* other ) const ;

      // if `self' equal to `other', O ; if smaller, <0 ; if greater, >0 
      // IMPLEMENTATION : address_comparison(other)
      virtual int three_way_comparison( PEL_Object const* other ) const ;

      // hash code value
      // IMPLEMENTATION : address()
      virtual size_t hash_code( void ) const ;

      // Is `other' identical to `self' (same address) ?
      bool has_same_address( PEL_Object const* other ) const ;

      // if `self''s address equal to the address of `other', O ; 
      // if smaller, -1 ; if greater, +1 */
      int address_comparison( PEL_Object const* other ) const ;

   //-- Persistence(900.0)

      // Use `writer' to store `self' so that ist can be retrieved
      // with `::restore_state'.
      // IMPLEMENTATION : raise a fatal error.
      virtual void save_state( PEL_ObjectWriter* writer ) const ;

      // Retrieve `self' from `reader'.
      // IMPLEMENTATION : raise a fatal error.
      virtual void restore_state( PEL_ObjectReader* reader ) ;

      // Perform updating to restore the consistency between `self' and a 
      // set of related objects (observer pattern). To be used instead
      // of `::update()' in `::restore_state()' because, in the course
      // of a restart, the operations performed in 
      // `::update_for_restore_state()' are closely related to those performed
      // in `::restore_state()'.
      // IMPLEMENTATION : call `::update()'.
      virtual void update_for_restore_state( void ) ;
      
   //-- Input - Output(1001.0)

      // Write text for debugging to `os'. 
      // IMPLEMENTATION : write the type name.
      virtual void display_info( std::ostream& os, size_t indent_width ) const ;

      // Write text to `os' with `indent_width' indentation.
      // IMPLEMENTATION : `os' is left unchanged.
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   //-- Hidden

      // Total number of PEL_Object actually allocated */
      static int GetNumberOf_PEL_objects( void ) ;

      // for debugging purposes only
      static std::ostream& TraceRemainingObjects( std::ostream& out ) ;

      // catch creation and deletion for particular object.
      static void catch_object( PEL_Object const* obj ) ;

      // catch creation for particular object by its creation rank.
      static void catch_object_by_rank( size_t rank ) ;

      static void start_trace_allocating( void ) ;
      static bool trace_allocating( void ) ;
      static void stop_trace_allocating( void ) ;
      static void trace_not_destroyed_object( std::ostream& out ) ;

   protected: //--------------------------------------------------------

   //-- Plug in(1002.0)

      virtual ~PEL_Object( void ) = 0 ;

      // Construction of an instance whose owner is `a_owner'
      PEL_Object( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant(9999.0)

      virtual bool invariant( void ) const ;

      bool create_clone_POST( PEL_Object const* result, 
                              PEL_Object const* a_owner ) const ;

      virtual bool save_state_PRE( PEL_ObjectWriter const* writer ) const ;
      virtual bool save_state_POST( PEL_ObjectWriter const* writer ) const ;

      virtual bool restore_state_PRE( PEL_ObjectReader const* reader ) const ;
      virtual bool restore_state_POST( PEL_ObjectReader const* reader ) const ;
      

      virtual bool is_equal_PRE( PEL_Object const* other ) const ;
      virtual bool is_equal_POST( bool result, 
                                  PEL_Object const* other ) const ;

      virtual bool three_way_comparison_PRE( PEL_Object const* other ) const ;
      virtual bool three_way_comparison_POST( int result, 
                                              PEL_Object const* other ) const ;

   private: //----------------------------------------------------------

      PEL_Object( void ) ;
      PEL_Object( PEL_Object const& other ) ;
      PEL_Object& operator=( PEL_Object const& other ) ;

      // Inserts a aComponent to the list of owned components.
      void insert_possession( PEL_Object* obj ) ;

   //-- Class attributes
      
      // Number of allocated PEL_Object instances remaining on the heap. 
      static size_t ALLOCATED ;
      
   //-- Attributes
      
      PEL_Object* MY_OWNER ;
      PEL_ListIdentity* POSSESSIONS ;

} ;


#endif

