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

#include <boolArray2D.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <stringArray2D.hh>

/*

Operator to form array from simple vectors.

PUBLISHED
*/

class PEL_EXPORT PEL_ArrayExp : public PEL_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value

      virtual doubleArray2D const& to_double_array2D(
                                        PEL_Context const* ct = 0 ) const ;

      virtual intArray2D const& to_int_array2D( 
                                        PEL_Context const* ct = 0 ) const ;
      
      virtual boolArray2D const& to_bool_array2D( 
                                        PEL_Context const* ct = 0 ) const ;
      
      virtual stringArray2D const& to_string_array2D( 
                                        PEL_Context const* ct = 0 ) const ;
      
      virtual doubleArray3D const& to_double_array3D(
                                        PEL_Context const* ct = 0 ) const ;

      virtual intArray3D const& to_int_array3D( 
                                        PEL_Context const* ct = 0 ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

     ~PEL_ArrayExp( void ) ;
      PEL_ArrayExp( PEL_ArrayExp const& other ) ;
      PEL_ArrayExp& operator=( PEL_ArrayExp const& other ) ;
      
      PEL_ArrayExp( PEL_Object* a_owner,
                    PEL_Sequence const* argument_list ) ;
      
   //-- Plug in

      PEL_ArrayExp( void ) ;

      virtual PEL_ArrayExp* create_replica( 
                                  PEL_Object * a_owner,
                                  PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments(
                                 PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes
            
      static PEL_ArrayExp const* PROTOTYPE ;

   //-- Attributes

      mutable doubleArray2D RESULT_D2D ;
      mutable intArray2D    RESULT_I2D ;
      mutable boolArray2D   RESULT_B2D ;
      mutable stringArray2D RESULT_S2D ;
      mutable doubleArray3D RESULT_D3D ;
      mutable intArray3D    RESULT_I3D ;
} ;

#endif
