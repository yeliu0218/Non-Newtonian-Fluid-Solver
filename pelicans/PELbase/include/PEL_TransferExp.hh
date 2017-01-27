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

#ifndef PEL_TRANSFER_EXP_HH
#define PEL_TRANSFER_EXP_HH

#include <PEL_Expression.hh>

/* Abstract class provide to transfer evaluation to context-dependant
   data.
   
   FRAMEWORK INSTANTIATION :
      Implement `::data' method.

PUBLISHED
*/

class PEL_EXPORT PEL_TransferExp : public PEL_Expression
{
   public: //----------------------------------------------------------
      
   //-- Value
      
      virtual bool to_bool( PEL_Context const* ct ) const ;
      virtual double to_double( PEL_Context const* ct ) const ;
      virtual int to_int( PEL_Context const* ct ) const ;
      virtual std::string const& to_string( PEL_Context const* ct ) const ;
      virtual doubleVector const& to_double_vector( PEL_Context const* ct ) const ;
      virtual intVector const& to_int_vector( PEL_Context const* ct ) const ;
      virtual stringVector const& to_string_vector(
                                              PEL_Context const* ct ) const ;
      virtual boolVector const& to_bool_vector( PEL_Context const* ct ) const ;
      virtual doubleArray2D const& to_double_array2D( PEL_Context const* ct ) const ;
      virtual intArray2D const& to_int_array2D( PEL_Context const* ct ) const ;
      virtual boolArray2D const& to_bool_array2D( 
                                        PEL_Context const* ct = 0 ) const ;      
      virtual stringArray2D const& to_string_array2D( 
                                        PEL_Context const* ct = 0 ) const ;
      virtual doubleArray3D const& to_double_array3D( PEL_Context const* ct ) const ;
      virtual intArray3D const& to_int_array3D( PEL_Context const* ct ) const ;
      
   protected: //-------------------------------------------------------
      
      virtual ~PEL_TransferExp( void ) ;

      PEL_TransferExp( std::string const& a_name ) ;

      PEL_TransferExp( PEL_Object* a_owner,
                      std::string const& a_name,
                      PEL_Sequence const* argument_list ) ;
     
   //-- Transfer implementation
      
      // `PEL_Data::' object depending on context `ct'
      virtual PEL_Data const* data( PEL_Context const* ct ) const = 0 ;
     
   //-- Preconditions, Postconditions, Invariant

      virtual bool data_PRE( PEL_Context const* ct ) const ;
      
      virtual bool data_POST( PEL_Data const* result,
                              PEL_Context const* ct ) const ;
      
   private: //-------------------------------------------------------

      PEL_TransferExp( void ) ;
      PEL_TransferExp( PEL_TransferExp const& other ) ;
      PEL_TransferExp& operator=( PEL_TransferExp const& other ) ;
      
};

#endif
