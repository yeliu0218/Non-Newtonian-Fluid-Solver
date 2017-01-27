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

#ifndef PDE_CROSS_PROCESS_NODE_NUMBERING_HH
#define PDE_CROSS_PROCESS_NODE_NUMBERING_HH

#include <PEL_Object.hh>

#include <intVector.hh>
#include <size_t_vector.hh>

#include <PDE_DomainAndFields.hh>

class PEL_VectorIterator ;
class boolVector ;
class doubleArray2D ;

class GE_Point ;

class PDE_DiscreteField ;

/*
Mappings between local (on-process) and the global (cross-process)
node numbering

Each node is handled by exactly one process, but can be accessed by
more than one process.
*/

class PEL_EXPORT PDE_CrossProcessNodeNumbering : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      void attach_field( PDE_DiscreteField* a_field ) ;

   //-- Characteristics
      
      // discrete field linked to `self'
      PDE_DiscreteField const* field( void ) const ;

      // communicator used by `self'
      PEL_Communicator const* communicator( void ) const ;

   //-- On-process to cross-process mapping
      
      // number of different global (cross-process) nodes
      size_t nb_global_nodes( void ) const ;
      
      // global (cross-process) index of the local (on-process) node `n'
      size_t global_node_index( size_t n ) const ;

      // local (on-process) node whose global (cross-process) index
      // is `global_index' (`PEL::bad_index()' if not in the current process)
      size_t local_node( size_t global_index ) const ;

      // Is `n' handled by the current process (each node is
      // handled by exactly one process, but may be viewed by other
      // processes) ? 
      bool current_process_handles_node( size_t n ) const ;

      // rank of process handling `n' (each node is
      // handled by exactly one process, but may be viewed by other
      // processes) 
      size_t rank_of_process_handling( size_t n ) const ;

   //-- Synchronization

      void synchronize_valid_nodes( void ) const ;

      void synchronize_imposed_values( void ) const ;

      // Synchronize nodes added to `::field()' by 
      // `PDE_DiscreteField::add_nodes()' since the last synchronization.
      void synchronize_new_nodes( void ) ;
      
   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      PDE_CrossProcessNodeNumbering( void ) ;
     ~PDE_CrossProcessNodeNumbering( void ) ;
      PDE_CrossProcessNodeNumbering( 
                            PDE_CrossProcessNodeNumbering const& other ) ;
      PDE_CrossProcessNodeNumbering& operator=(
                            PDE_CrossProcessNodeNumbering const& other ) ;

      PDE_CrossProcessNodeNumbering( PEL_Object* a_owner,
                                     PDE_DomainBuilder const* a_dom,
                                     PEL_Communicator const* a_com ) ;
      
      friend  PDE_CrossProcessNodeNumbering* 
      PDE_DomainAndFields::create_CrossProcessNodeNumbering( 
                                              PEL_Object* a_owner) const ;
      
      void globalize( void ) ;
      
      void recover_coordinates( doubleArray2D& coord,
                                boolVector& is_halo ) const  ;
      
      void recover_coordinates( PEL_VectorIterator* mesh_it,
                                doubleArray2D& coord,
                                boolVector& is_halo ) const ;

      void recover_global_numbering(
                                doubleArray2D& coord,
                                boolVector const& is_halo,
                                size_t_vector& global_node,
                                size_t_vector& node_owner,
                                bool verbose ) const ;

      void set_global_to_local( size_t_vector const& global_nodes,
                                size_t_vector const& nodes_owner,
                                size_t_vector& local_nodes ) const ;

  //-- Attributes

      PDE_DomainBuilder const* DOM ;
      size_t DIM ;
      GE_Point* PT ;
      size_t_vector LOCAL_NODE ;
      size_t_vector GLOBAL_NODE ;
      size_t_vector NODE_OWNER ;
      PDE_DiscreteField* FIELD ;
      PEL_Communicator const* COMM ;
} ; 

#endif 
