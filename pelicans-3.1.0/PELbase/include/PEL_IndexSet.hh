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

#ifndef PEL_INDEX_SET_HH
#define PEL_INDEX_SET_HH

#include <PEL_Object.hh>

#include <size_t_vector.hh>

class PEL_EXPORT PEL_IndexSet : public PEL_Object
{
      
   public : //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_IndexSet* create( PEL_Object* a_owner ) ;
      
      static PEL_IndexSet* create( PEL_Object* a_owner,
                                   size_t_vector const& vec,
                                   size_t a_id ) ;
      
      void re_initialize( size_t_vector const& vec, 
                          size_t a_id ) ;
      
   //-- Identifier

      size_t id( void ) const ;
      
   //-- Element access

      // elements in increasingly order
      size_t_vector const& elements( void ) const ;

   //-- Comparison

      virtual bool is_equal( PEL_Object const* other ) const ;
      
      virtual int three_way_comparison( PEL_Object const* other ) const ;
      
      virtual size_t hash_code( void ) const ;    

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected ://------------------------------------------------------------

   private :  //------------------------------------------------------------

      PEL_IndexSet( void ) ;
     ~PEL_IndexSet( void ) ;
      PEL_IndexSet( PEL_IndexSet const& other ) ;
      PEL_IndexSet operator=( PEL_IndexSet const& other ) ;

      PEL_IndexSet( PEL_Object* a_owner ) ;
      
      PEL_IndexSet( PEL_Object* a_owner,
                    size_t_vector const& vec,
                    size_t a_id ) ;

   //-- Attributes

      size_t ID ;
      size_t_vector SET ;
} ;

#endif
