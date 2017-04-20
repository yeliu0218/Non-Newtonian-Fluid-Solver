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

#ifndef PEL_OBJECT_READER_HH
#define PEL_OBJECT_READER_HH

#include <PEL_Object.hh>

#include <stack>
#include <string>

class PEL_Data ;
class PEL_Module ;
class PEL_ModuleExplorer ;
class PEL_ModuleIterator ;

/*
Servers used to retrieve objects stored with
associated `PEL_ObjectWriter::' instances.

PUBLISHED
*/

class PEL_EXPORT PEL_ObjectReader : public PEL_Object
{
   public: //----------------------------------------------------------------

   //-- Instance delivery and initialization
      
      static PEL_ObjectReader* create( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) ;

   //-- Cycles

      PEL_Module* header_module( void ) const ;

      size_t nb_cycles( void ) const ;

      void seek_cycle( size_t cycle_number ) ;

      void close_cycle( void ) ;

      bool positioned_in_a_valid_cycle( void ) const ;

   //-- Object retrieval

      void start_object_retrieval( std::string const& class_name ) ;

      size_t current_object_number( void ) const ;

      bool has_entry( std::string const& keyword ) const ;
      PEL_Data const* data_of_entry( std::string const& keyword ) const ;

      void end_object_retrieval( void ) ; 
    
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

      PEL_ObjectReader( void ) ;
     ~PEL_ObjectReader( void ) ;
      PEL_ObjectReader( PEL_ObjectReader const& other ) ;
      PEL_ObjectReader& operator=( PEL_ObjectReader const& other ) ;

      PEL_ObjectReader( PEL_Object* a_owner,
                        PEL_ModuleExplorer const* exp ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      PEL_Module const* ROOT_MOD ;
      size_t NB_CYCLES ;
      size_t LAST_CYCLE ;
      int iOBJECT ;
      std::stack< PEL_Module const* > MODS ;
      std::stack< PEL_ModuleIterator* > MOD_ITS ;
} ;

#endif
