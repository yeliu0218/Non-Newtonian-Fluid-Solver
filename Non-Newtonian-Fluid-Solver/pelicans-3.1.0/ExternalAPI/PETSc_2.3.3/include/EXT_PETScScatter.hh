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

#ifndef EXT_PETScScatter_HH
#define EXT_PETScScatter_HH

#include <LA_Scatter.hh>

#include <EXT_PETScAPI.hh>
#include <intVector.hh>
#include <size_t_vector.hh>

class EXT_PETScVector ;
class LA_SeqVector ;

/*
  `LA_Scatter::' server for `EXT_PETScImplementation::' objects.
  
  PUBLISHED
*/

class EXT_PETScScatter : public LA_Scatter
{
   public: //-----------------------------------------------------------

      static EXT_PETScScatter* create(
                     PEL_Object* a_owner,
                     EXT_PETScVector const* global_vector,
                     size_t_vector const& a_repatriated_items_table,
                     size_t_vector const& a_local_indices_table ) ;
      
   //-- Characteristics
      
      virtual LA_Implementation const* implementation( void ) const ;
      
      virtual size_t size( void ) const ;

      virtual size_t_vector const& repatriated_items( void ) const ;
      
      virtual size_t_vector const& local_indices( void ) const ;
      
      virtual PEL_DistributedPartition const* distribution( void ) const ;
      
   //-- Repatriate items
      
      virtual void get( LA_Vector const* source,
                        LA_SeqVector* dest ) const ;
      
      virtual void set( LA_SeqVector const* source,
                        LA_Vector* dest ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      EXT_PETScScatter( void ) ;
     ~EXT_PETScScatter( void ) ;
      EXT_PETScScatter( EXT_PETScScatter const& other ) ;
      EXT_PETScScatter& operator=( EXT_PETScScatter const& other ) ;

      EXT_PETScScatter( PEL_Object* a_owner,
                        EXT_PETScVector const* global_vector,
                        size_t_vector const& a_repatriated_items_table,
                        size_t_vector const& a_local_indices_table ) ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool implementation_POST(
                            LA_Implementation const* result ) const ;
      
   //-- Attributes
      
      size_t const SIZE ;
      PEL_DistributedPartition* const DIST ;
      size_t_vector const NEEDED ;
      size_t_vector const LOCAL ;
      
      VecScatter SCATTER ;
      Vec SEQ ;
} ;

#endif
