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

#ifndef GE_SPLITTING_STRATEGY_HH
#define GE_SPLITTING_STRATEGY_HH

#include <PEL_Object.hh>

class GE_Meshing ;

class PEL_Communicator ;
class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

/*
Strategies used by `GE_SplitMeshing' to split a given meshing and distribute
the pieces on several processes.

FRAMEWORK INSTANTIATION:

  1. Derive a concrete subclass, say `MyStrategy'.
  2. Choose a name for `MyStrategy', say "MyStrategy".
  3. Implement a private destructor.
  4. Declare all constructors private.
  5. Define the prototype to be registered:
      5.1 Implement a default constructor that initializes the
          `GE_SplittingStrategy::' subobject by calling
          `::GE_SplittingStrategy( std::string const& )'
          with "MyStrategy" as argument.
      5.2 Define and initialize a static instance by calling the default
          constructor.
          declaration:
               `static MyStrategy const* prototype ;'
          definition :
              `MyStrategy const* MyStrategy::prototype = new MyStrategy() ;'
  6. Implement the `create_replica' methods.
  7. Implement the constructor called by `create_replica' methods.
     The `MyStrategy'  subobject is initialized by calling
        `GE_SplittingStrategy'( PEL_Object*, ...)
  8. Implement `::cell_rank' pure virtual method.

PUBLISHED
*/

class PEL_EXPORT GE_SplittingStrategy : public PEL_Object
{
   public: //------------------------------------------------------------
      
   //-- Instance delivery and initialization
      
      // Create and return an instance devoted to the distribution of `meshing'
      // on the set of processes defined by `com'.
      static GE_SplittingStrategy* create( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp,
                                           GE_Meshing* meshing,
                                           PEL_Communicator const* com ) ;

      // Create and return an instance for testing purposes: the distribution
      // of `meshing' is simulated using only one process.
      //   - the number of processes is simulated to be `nb_rks'
      //   - the current process is simulated to be of rank `rk'.
      static GE_SplittingStrategy* create( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp,
                                           GE_Meshing* meshing,
                                           size_t nb_rks, size_t rk ) ;
      
   //-- Processes(200.)

      size_t nb_ranks( void ) const ;

      size_t rank( void ) const ;
      
   //-- Cell balance(300.)

      size_t nb_cells( void ) const ;
      
      // rank of the process which owns the cell of index `mesh_id', this
      // index being related to the implicit numbering defined by the
      // traversal order of the cell-iterator in `GE_Meshing::'
      virtual size_t cell_rank( size_t mesh_id ) const = 0 ;

   protected: //---------------------------------------------------------

   //-- Plug in
      
      virtual ~GE_SplittingStrategy( void ) ;

      // Registration of an instance as `name'.
      GE_SplittingStrategy( std::string const& name ) ;
      
      // In the constructor called by `::create_replica(com)'.
      GE_SplittingStrategy( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp,
                            GE_Meshing* meshing,
                            PEL_Communicator const* com ) ;
      
      // In the constructor called by `::create_replica(nb_rks,rk)'.
      GE_SplittingStrategy( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp,
                            GE_Meshing* meshing,
                            size_t nb_rks, size_t rk ) ;
      
      bool is_a_prototype( void ) const ;

      virtual GE_SplittingStrategy* create_replica(
                            PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp,
                            GE_Meshing* meshing,
                            size_t nb_rks, size_t rk ) const = 0 ;

      virtual GE_SplittingStrategy* create_replica(
                            PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp,
                            GE_Meshing* meshing,
                            PEL_Communicator const* com ) const = 0 ;
      
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
      
      bool cell_rank_PRE( size_t mesh_id ) const ;
      bool cell_rank_POST( size_t mesh_id, size_t result ) const ;
      
      bool create_replica_PRE(  PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp,
                                GE_Meshing* meshing,
                                PEL_Communicator const* com ) const ;
      bool create_replica_POST( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp,
                                GE_Meshing* meshing,
                                PEL_Communicator const* com,
                                GE_SplittingStrategy const* result ) const ;
      
      bool create_replica_PRE(  PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp,
                                GE_Meshing* meshing,
                                size_t nb_rks, size_t rk ) const ;
      bool create_replica_POST( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp,
                                GE_Meshing* meshing,
                                size_t nb_rks, size_t rk,
                                GE_SplittingStrategy const* result ) const ;
    
   private: //-----------------------------------------------------------

      GE_SplittingStrategy( void ) ;
      GE_SplittingStrategy( GE_SplittingStrategy const& other ) ;
      GE_SplittingStrategy& operator=( GE_SplittingStrategy const& other ) ;
      
      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes
      
      bool const IS_PROTO ;

      size_t const NB_CELLS ;
      size_t const NB_RANKS ;
      size_t const RANK ;
} ;

#endif


