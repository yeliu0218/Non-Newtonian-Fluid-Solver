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

#ifndef PDE_BASIS_FUNCTION_BOUND_HH
#define PDE_BASIS_FUNCTION_BOUND_HH

#include <PDE_BasisFunction.hh>

class PEL_EXPORT PDE_BasisFunctionBound : public PDE_BasisFunction
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_BasisFunctionBound* create( PEL_Object* a_owner ) ;

   //-- Activation and refinement

      virtual size_t refinement_level( void ) const ;

   //-- Parents

      size_t nb_parents( void ) const ;

      PDE_BasisFunctionBound* parent( size_t i ) const ;

      // notify that `self' is a child of `a_parent' (`a_parent' is modified)
      // and that `a_parent' is a parent of `self'
      void set_child_parent_relationship( PDE_BasisFunctionBound* a_parent,
                                          double refinement_coef ) ;

   //-- Childs

      size_t nb_childs( void ) const ;

      PDE_BasisFunctionBound* child( size_t i ) const ;

      double refinement_coefficient( size_t i ) const ;

   //-- Pieces
 
      void extend_pieces( PDE_BoundFE* a_bound,
                          size_t elm_index,
                          size_t node_in_elm ) ;

      size_t nb_bounds( void ) const ;

      // support of the `i'-th piece
      PDE_BoundFE* bound( size_t i ) const ;

      // reference element of the `i'-th piece
      size_t element_index_of_bound( size_t i ) const ;

      // local node index in the `i'-th piece
      size_t local_node_of_bound( size_t i ) const ;

   //-- Location

      virtual bool is_located_in_cell( PDE_CellFE const* cell ) const ;

      virtual bool is_located_on_bound( PDE_BoundFE const* bound ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_BasisFunctionBound( void ) ;
     ~PDE_BasisFunctionBound( void ) ;
      PDE_BasisFunctionBound( PDE_BasisFunctionBound const& other ) ;
      PDE_BasisFunctionBound& operator=( 
                              PDE_BasisFunctionBound const& other ) ;

      PDE_BasisFunctionBound( PEL_Object* a_owner ) ;

      bool check_location_on_bound( PDE_BoundFE const* a_bound ) const ;

   //-- Attributes

      std::vector< PDE_BoundFE* > BOUNDS ;
      std::vector< size_t > ELM_IDXS ;
      std::vector< size_t > ELM_NODES ;

      std::vector< PDE_BasisFunctionBound* > PARENTS ;
      std::vector< double > REFI_COEFS ;
      std::vector< PDE_BasisFunctionBound* > CHILDS ;
} ;

#endif
