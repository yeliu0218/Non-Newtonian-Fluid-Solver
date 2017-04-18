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

#ifndef PDE_SET_OF_BASIS_FUNCTIONS_HH
#define PDE_SET_OF_BASIS_FUNCTIONS_HH

#include <PEL_Object.hh>

#include <size_t_vector.hh>
#include <vector>

class PDE_BasisFunction ;
class PDE_BasisFunctionCell ;
class PDE_CrossProcessBFNumbering ;

class PEL_IndexSet ;
class PEL_Vector ;

class size_t_vector ;
class stringVector ;

class PEL_EXPORT PDE_SetOfBasisFunctions : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      static PDE_SetOfBasisFunctions* create( PEL_Object* a_owner ) ;

   //-- Reference elements groups
      
      void extend_ref_elts_grp( PEL_Vector const* elms, size_t& e ) ;

      size_t nb_ref_elts_grps( void ) const ;
       
   //-- Basis functions
      
      // total number of instances of `PDE_BasisFunction' included
      size_t nb_basis_functions( void  ) const ;

      // number of instances of `PDE_BasisFunction' with associated
      // reference elements group index `e'
      size_t nb_basis_functions( size_t e ) const ;

      // Does self include an instance of PDE_BasisFunction whose
      // id_number is  equal to `id_number'
      // and associated reference elements group index is equal to `e' ?
      bool has( size_t id_number, size_t e ) const ;

      // the included instance of PDE_BasisFunction whose
      // id_number is  equal to `id_number'
      // and associated reference elements group index is equal to `e'
      PDE_BasisFunction* item( size_t id_number, size_t e ) const ;

      // Add `bf' and attribute an id_number to `bf'.
      void add( PDE_BasisFunctionCell* bf, size_t e ) ;

   //-- Iteration over the included instances of PDE_BasisFunction whose
   //-- associated reference element index is equal to `e'

      // Move iterator to the first position.
      void start( size_t e ) const ;

      // Is iterator position valid ?
      bool is_valid( size_t e ) const ;

      // Move iterator one position.
      void go_next( size_t e ) const ;

      // Instance of PDE_BasisFunction at current iterator position
      PDE_BasisFunction* item( size_t e ) const ;

   //-- Distributed processing

      // Does `self' represent a part of a logical global discrete field
      // involved in a processing distributed over several processes ?
      bool is_distributed( void ) const ;

      // cross-process numbering
      PDE_CrossProcessBFNumbering* cross_process_numbering( void ) const ;

      // Perform all necessary tasks to reinterpret `self' as a part of
      // a logical global discrete field in a processing distributed over
      // `numbering'->communicator().
      void build_cross_process_globalization(
                                PDE_CrossProcessBFNumbering* numbering ) ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_SetOfBasisFunctions( void ) ;
     ~PDE_SetOfBasisFunctions( void ) ;
      PDE_SetOfBasisFunctions( PDE_SetOfBasisFunctions const& other ) ;
      PDE_SetOfBasisFunctions& operator=(
                               PDE_SetOfBasisFunctions const& other ) ;

      PDE_SetOfBasisFunctions( PEL_Object* a_owner ) ;
      
      // Return the index of `ref_elts' (`PEL::bad_index()' if not included)
      size_t index_of_ref_elts_grp( PEL_IndexSet const* ref_elts ) const ;

   //-- Attributes
      
      std::vector< PEL_IndexSet const* > REF_ELTS  ;
      
      PDE_CrossProcessBFNumbering* NUMBERING ;

      std::vector< PEL_Vector* > BFS ;
      mutable size_t_vector BFS_I ;
} ;

#endif
