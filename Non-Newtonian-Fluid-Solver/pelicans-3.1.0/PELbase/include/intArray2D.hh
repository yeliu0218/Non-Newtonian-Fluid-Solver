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

#ifndef INT_ARRAY_2D_HH
#define INT_ARRAY_2D_HH

#include <intVector.hh>

class size_t_array2D ;

/*
sequences of values, all of type int, ordered according to two indices
in contiguous intervals
*/

class PEL_EXPORT intArray2D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0' and `dim1' as an 
      // exclusive upper limit.
      intArray2D( size_t dim0, size_t dim1, int val=0 ) ;

      intArray2D( intArray2D const& other ) ;

      intArray2D const& operator=( intArray2D const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, int val=0 ) ;

   //-- Termination

     ~intArray2D( void ) ;

   //-- Comparison
      
      bool operator==( intArray2D const& other ) const ;
      bool operator!=( intArray2D const& other ) const ;
      
   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;

      // Reinitialize `x' by copying all items of self whose `an_index'-th 
      // index is equal to `index-value'.
      void extract_section( size_t an_index,
                            size_t index_value,
                            intVector& x ) const ;

      // item of indices `i0' and `i1'
      int const& operator()( size_t i0, size_t i1 ) const ;

    //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, int val=0 ) ;

      // Assign `val' to all items.
      void set( int val ) ; 

      void set( size_t_array2D const& other ) ;

      // Assign `x' items to the items of self whose `an_index'-th index is 
      // equal to `index-value' and whose other index is given by `x'.
      void set_section( size_t an_index, 
                       size_t index_value, 
                       intVector const& x ) ;

      // item of indices `i0' and `i1'
      int& operator()( size_t i0, size_t i1 ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       intArray2D const& a ) ;

   //-- Hidden

      // pointer to internal data
      // (for time optimization, should be used with extreme care)
      int const* data( void ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      intArray2D( void ) ;

   //-- Attributes
  
      intVector vector ;
      size_t d0 ; 
      size_t d1 ;
             
} ;

#ifndef OUTLINE
#include <intArray2D.icc>
#endif

#endif