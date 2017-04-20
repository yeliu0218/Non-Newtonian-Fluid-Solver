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

#ifndef LA_SORTING_HH
#define LA_SORTING_HH

#include <size_t_vector.hh>

class LA_SeqVector ;
class LA_DenseMatrix ;
class LA_SeqMatrix ;
class PEL_Object ;
class size_t_vector ;
class size_t_array2D ;
class intArray2D ;
class intVector ;

/**
  Sort of the values of a vector
*/
class PEL_EXPORT LA_Sorting
{
   public: //-----------------------------------------------------------

  //-- Sorting
      // Sort `values' in `sortedValues' with insert method.
      // `rank' is the corresponding permutation.
      // This sorting method is stable.
      static void sort( LA_SeqVector const* values,
                        LA_SeqVector* sortedValues,
                        size_t_vector& rank ) ;

      // Sort `values' in `sortedValues' with insert method.
      // `rank' is the corresponding permutation.
      // This sorting method is stable.
      static void sort( intVector const& values,
                        intVector& sortedValues,
                        size_t_vector& rank ) ;

      // Sort `sortedValues' in place with quick sort method.
      // This sorting method is not stable.
      static void quick_sort( LA_SeqVector* sortedValues,
                              size_t start,
                              size_t end ) ;
   //-- Ranking
      
      // Fill in  matrix representing rank of items of `mat' in
      // their respective column using `sort' method.
      static void column_rank( LA_SeqMatrix const* mat,
                               size_t_array2D& rank ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      static void exchange( LA_SeqVector* sortedValues,
                            size_t a,
                            size_t b ) ;
      LA_Sorting( void ) ;
      LA_Sorting( LA_Sorting const& other ) ;
      LA_Sorting const& operator=( LA_Sorting const& other ) ;
     ~LA_Sorting( void ) ;

} ;


#endif
