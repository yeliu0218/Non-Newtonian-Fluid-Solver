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

#ifndef LA_SCATTER_HH
#define LA_SCATTER_HH

#include <PEL_Object.hh>

class PEL_DistributedPartition ;
class size_t_vector ;
class LA_Implementation ;
class LA_SeqVector ;
class LA_Vector ;

/*
  Objects used to repatriate distributed items in local one.

  PUBLISHED
*/

class PEL_EXPORT LA_Scatter : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Characteristics

      // implementation indicator
      virtual LA_Implementation const* implementation( void ) const = 0 ;
      
      virtual size_t size( void ) const = 0 ;

      // table of items which needed to be repatriated
      virtual size_t_vector const& repatriated_items( void ) const = 0 ;
      
      // indices of `::repatriated_items' elements in the sequential vector
      virtual size_t_vector const& local_indices( void ) const = 0 ;
      
      // distribution of vector rows over the processes
      virtual PEL_DistributedPartition const* distribution( void ) const = 0 ;

   //-- Repatriate items
      
      // Repatriate distributed items in distributed `source' to local `dest'.
      virtual void get( LA_Vector const* source,
                        LA_SeqVector* dest ) const = 0 ;
      
      // Repatriate distributed items in local `source' to distributed `dest'.
      virtual void set( LA_SeqVector const* source,
                        LA_Vector* dest ) const = 0 ;
      
   protected: //--------------------------------------------------------

      virtual ~LA_Scatter( void ) ;
      LA_Scatter( PEL_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool implementation_POST(
                            LA_Implementation const* result ) const ;

      virtual bool repatriated_items_POST(
                                size_t_vector const& result ) const ;

      virtual bool local_indices_POST(
                                size_t_vector const& result ) const ;
      
      virtual bool distribution_POST(
                     PEL_DistributedPartition const* result ) const ;
      
      virtual bool get_PRE( LA_Vector const* source,
                            LA_SeqVector const* dest ) const ;
      virtual bool get_POST( LA_Vector const* source,
                            LA_SeqVector const* dest ) const ;
      
      virtual bool set_PRE( LA_SeqVector const* source,
                            LA_Vector const* dest ) const ;
      virtual bool set_POST( LA_SeqVector const* source,
                            LA_Vector const* dest ) const ;
      
   private: //----------------------------------------------------------

      LA_Scatter( void ) ;
      LA_Scatter( LA_Scatter const& other ) ;
      LA_Scatter& operator=( LA_Scatter const& other ) ;

} ;

#endif
