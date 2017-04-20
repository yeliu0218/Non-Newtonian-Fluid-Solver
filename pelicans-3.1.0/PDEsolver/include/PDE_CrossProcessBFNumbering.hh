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

#ifndef PDE_CROSS_PROCESS_BF_NUMBERING_HH
#define PDE_CROSS_PROCESS_BF_NUMBERING_HH

#include <PEL_Object.hh>

#include <intVector.hh>
#include <size_t_vector.hh>
#include <vector>

#include <PDE_DomainAndFields.hh>

class PEL_VectorIterator ;
class boolVector ;
class doubleArray2D ;

class GE_Point ;
class PEL_Communicator ;
class PDE_SetOfBasisFunctions ;

/*
Mappings between local (on-process) and the global (cross-process)
basis functions numbering.

Each basis function is handled by exactly one process,
but can be accessed by more than one process.
*/

class PEL_EXPORT PDE_CrossProcessBFNumbering : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

   //-- Characteristics

      //attached set of reference elements groups `PDE_SetOfBasisFunctions'
      PDE_SetOfBasisFunctions const* set_of_basis_functions( void ) const ;

      // communicator used by `self'
      PEL_Communicator const* communicator( void ) const ;

   // Initialization

      // resize internal data structures
      // may be called when attached_set_of_basis_functions()->
      // attached_ref_elts_grps()->nb_ref_elts_grps() change
      void resize_internals( void ) ;

   //-- On-process to cross-process mapping

      // number of different global (cross-process) basis functions
      size_t nb_global_basis_functions( void ) const ;

      // number of different global (cross-process) basis functions
      // whose reference elements group index is `e' in the attached
      // `PDE_SetOfReferenceElements'
      size_t nb_global_basis_functions( size_t e ) const ;

      // global (cross-process) index
      // of basis function whose local (on-process) index is `local_ind'
      size_t global_basis_function_index( size_t local_ind, size_t e ) const ;

      // local (on-process) index
      // of basis function whose global (cross-process) index
      // is `global_index' (`PEL::bad_index()' if not in the current process)
      size_t local_basis_function_index( size_t global_ind, 
                                         size_t& local_e ) const ;

   //-- Synchronization

      void synchronize_new_basis_functions( void ) ;

      void globalize( void ) ;

      bool is_synchronized( void ) const ;

      void set_unsynchronized_state( void ) ;

   //-- Preconditions, Postconditions, Invariant
      virtual bool invariant( void ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_CrossProcessBFNumbering( void ) ;
     ~PDE_CrossProcessBFNumbering( void ) ;
      PDE_CrossProcessBFNumbering(
                            PDE_CrossProcessBFNumbering const& other ) ;
      PDE_CrossProcessBFNumbering& operator=(
                            PDE_CrossProcessBFNumbering const& other ) ;

      PDE_CrossProcessBFNumbering( PEL_Object* a_owner,
                                   PDE_DomainBuilder const* a_dom,
                                   PEL_Communicator const* a_com ) ;

      friend PDE_CrossProcessBFNumbering*
      PDE_DomainAndFields:: create_CrossProcessBFNumbering(
                                       PEL_Object* a_owner ) const ;


      void recover_coordinates( PDE_SetOfBasisFunctions const* bf_set,
                                size_t e,
                                doubleArray2D& coord ) const ;

      void recover_global_numbering( doubleArray2D& coord,
                                     size_t_vector& global_bf_id,
                                     size_t& nb_glob_bf,
                                     bool verbose ) const ;

      void set_global_to_local( size_t_vector const& global_bf_id,
                                size_t e,
                                size_t_vector& local_bf_id ) const ;

  //-- Attributes

      size_t DIM ;
      PDE_DomainBuilder const* DOMB ;
      GE_Point* PT ;
      std::vector<size_t_vector> LOCAL_BF ;
      std::vector<size_t_vector> GLOBAL_BF ;

      PEL_Communicator const* COMM ;
      size_t_vector NB_GLOBAL_BF ;

      PDE_SetOfBasisFunctions* BFS ;

      bool IS_SYNC ;
} ;

#endif
