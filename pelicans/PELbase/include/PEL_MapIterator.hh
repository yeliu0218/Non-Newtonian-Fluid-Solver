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

#ifndef PEL_MAP_ITERATOR_HH
#define PEL_MAP_ITERATOR_HH

#include <PEL_Iterator.hh>

class PEL_Map ;
class PEL_Iterator ;

//--------------------------------------------------------------------------
// Iterators on items of PEL_Map containers
//--------------------------------------------------------------------------
// The sequence of access is the internal ordering of items within the map :
// the first bucket is traversed first, then the second and so on.
//--------------------------------------------------------------------------

class PEL_EXPORT PEL_MapIterator : public PEL_Iterator
{
   public: //---------------------------------------------------------------

   //-- Initialization

      // Create an return an iterator over items of `map'
      static PEL_MapIterator* create( PEL_Object* a_owner,
                                      PEL_Map const* map ) ;

   //-- Status report

      virtual bool is_valid( void ) const ;

   //-- Access

      virtual PEL_Object* item( void ) const ;

      // key of the item at current position
      PEL_Object* key( void ) const ;

   //-- Cursor movement

      virtual void start( void ) ;

      virtual void go_next( void ) ;
       
   protected: //------------------------------------------------------------

      virtual ~PEL_MapIterator( void ) ;

      PEL_MapIterator( PEL_Object* a_owner, 
                       PEL_Map const* hTable ) ;
      
   private: //-------------------------------------------------------------- 

      PEL_MapIterator( void  ) ;
      PEL_MapIterator( PEL_MapIterator const& other ) ;
      PEL_MapIterator const& operator=( PEL_MapIterator const& other ) ;

      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      PEL_Map const* table ;
      PEL_Iterator* where ;
} ;

#endif
