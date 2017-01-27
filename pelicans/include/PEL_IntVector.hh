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

#ifndef PEL_INTVECTOR_HH
#define PEL_INTVECTOR_HH

#include <PEL_Data.hh>
#include <intVector.hh>

/* Implements list of int constant.
   The list can be built by adding simple value to it then finalizes it
   or from a intVector.  */

class PEL_Container ;

class PEL_EXPORT PEL_IntVector : public PEL_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_IntVector* create( PEL_Object* a_owner,
                                       intVector const& aIntVector ) ;
      
      static PEL_IntVector* create( PEL_Object* a_owner,
                                       PEL_Container const* aDataList ) ;
      
      virtual PEL_IntVector* create_clone( PEL_Object* a_owner ) const ;
      
   //-- Type

      virtual PEL_Data::Type data_type( void ) const ;

   //-- Value

      virtual intVector const& to_int_vector( 
                                            PEL_Context const* ct = 0 ) const ;  
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Modifier

      void set( intVector const& other ) ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PEL_IntVector( void ) ;
     ~PEL_IntVector( void ) ;
      PEL_IntVector( PEL_IntVector const & other ) ;
      PEL_IntVector& operator=( PEL_IntVector const& other ) ;

      PEL_IntVector( PEL_Object* a_owner,
                      intVector const& aIntVector ) ;

      PEL_IntVector( PEL_Object* a_owner,
                        PEL_Container const* aDataList,
                        bool& error ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      intVector dv ;
} ;


#endif
