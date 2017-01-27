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

#ifndef BOOL_VECTOR_HH
#define BOOL_VECTOR_HH

#include <PEL_assertions.hh>
#include <PEL_export.hh>

/*
sequences of values, all of type bool, ordered according to one index
in a contiguous interval

The values are stored in a block of contiguous storage. Depending on 
the use case, insertion of elements can cause reallocation; that is, the
sequence of values may be stored in a different area after insertion. The
reason is that if there is no room left in the current block for the new
elements, then a larger block is allocated, all of the old elements and the
new elements are copied to the new block, and the old block is deallocated.

Some control can be exercised over when reallocation occurs.

PUBLISHED */

class PEL_EXPORT boolVector
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      boolVector( size_t dim, bool val=false ) ;

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      boolVector( int dim, bool val=false ) ;

      boolVector( boolVector const& other ) ;

      // Assign `other' to `self' 
      // (causes reallocation iff `other.size()'>`capacity()').
      boolVector const& operator=( boolVector const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed 
      // (causes reallocation iff `size()' and `dim' differ ). 
      void re_initialize( size_t dim, bool val=false ) ;

   //-- Termination

     ~boolVector( void ) ;

   //-- Comparison

     bool operator==( boolVector const& other ) const ;
     
     bool operator!=( boolVector const& other ) const ;
       
   //-- Access

      // size of the currently allocated block
      size_t capacity( void ) const ;
      
      // exclusive upper limit of the index interval
      size_t size( void ) const ;

      // item of index i
      bool const& operator()( size_t i ) const ;
      
      // Does `val' appears in `self' ?
      bool has( bool val ) const ;
      
   //-- Element change

      // Change the exclusive upper limit of the index interval to `dim'
      // (with reallocation if `dim'>`capacity()').
      void resize( size_t dim, bool val=false ) ;
 
      // Assign `val' to all items.
      void set( bool val ) ; 

      // Include `val' at the end of `self'  (with possible reallocation).
      void append( bool val ) ;

      // Remove item at place `idx'.
      void remove_at( size_t idx ) ;

      // item of index i
      bool& operator()( size_t i ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       boolVector const& vec ) ;
      
      friend std::istream& operator>>( std::istream& in, 
                                       boolVector& vec ) ;
                  
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      boolVector( void ) ;
      boolVector( double dim ) ;
      boolVector( bool dim ) ;
      boolVector( char dim ) ;

      void allocate( size_t a_size ) ;
      void deallocate( void ) ;

      bool invariant( void ) const ;

   //-- Attributes

      bool* VECTOR ;
      size_t LENGTH ;
      size_t CAPACITY ;              
} ;


#ifndef OUTLINE
#include <boolVector.icc>
#endif

#endif
