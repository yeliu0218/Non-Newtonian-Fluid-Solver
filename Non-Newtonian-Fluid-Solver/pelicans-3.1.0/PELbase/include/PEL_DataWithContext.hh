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

#ifndef PEL_DATA_WITH_CONTEXT_HH
#define PEL_DATA_WITH_CONTEXT_HH

#include <PEL_Data.hh>

class PEL_Context ;
class PEL_ContextPair ;

class PEL_EXPORT PEL_DataWithContext : public PEL_Data
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      static PEL_DataWithContext* create( PEL_Object* a_owner,
                                          PEL_Data const* data,
                                          PEL_Context const* ct ) ;

      virtual PEL_DataWithContext* create_clone( PEL_Object* a_owner ) const ;
      
   //-- Context
                  
      virtual void declare( PEL_List* lst ) const ;

   //-- Type

      virtual Type data_type( void ) const ;
      
   //-- Value
      
      PEL_Context const* context( PEL_Context const* ct = 0 ) const ;
      
      virtual bool value_can_be_evaluated( PEL_Context const* ct = 0 ) const ;
      
      virtual stringVector const& undefined_variables(
                                           PEL_Context const* ct = 0 ) const ;

      virtual bool to_bool( PEL_Context const* ct = 0 ) const ;

      virtual double to_double( PEL_Context const* ct = 0 ) const ;

      virtual int to_int( PEL_Context const* ct = 0 ) const ;

      virtual std::string const& to_string(
                                        PEL_Context const* ct = 0 ) const ;

      virtual doubleVector const& to_double_vector(
                                        PEL_Context const* ct = 0 ) const ;

      virtual intVector const& to_int_vector(
                                        PEL_Context const* ct = 0 ) const ;

      virtual stringVector const& to_string_vector(
                                        PEL_Context const* ct = 0 ) const ;

      virtual boolVector const& to_bool_vector(
                                        PEL_Context const* ct = 0 ) const ;

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
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Formal calculus

      virtual bool is_raw_data( void ) const ;

   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------
      
      PEL_DataWithContext( void ) ;
     ~PEL_DataWithContext( void ) ;
      PEL_DataWithContext( PEL_DataWithContext const& other ) ;
      PEL_DataWithContext& operator=( PEL_DataWithContext const& other ) ;

      PEL_DataWithContext( PEL_Object* a_owner,
                           PEL_Data const* data,
                           PEL_Context const* ct ) ;

   //-- Attributes

      PEL_Data const* DATA ;
      PEL_Context* CTX ;
      PEL_ContextPair* TMP_CTX ;
      
      
};

#endif
