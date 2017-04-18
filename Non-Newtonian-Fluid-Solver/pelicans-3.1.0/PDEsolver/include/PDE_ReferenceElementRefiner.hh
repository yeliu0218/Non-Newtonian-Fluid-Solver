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

#ifndef PDE_REFERENCE_ELEMENT_REFINER_HH
#define PDE_REFERENCE_ELEMENT_REFINER_HH

#include <PEL_Object.hh>

#include <doubleArray3D.hh>
#include <size_t_array2D.hh>
#include <size_t_array3D.hh>

class PEL_Vector ;

class GE_Point ;
class GE_ReferencePolyhedron ;
class GE_ReferencePolyhedronRefiner ;

class PDE_ReferenceElement ;

class PEL_EXPORT PDE_ReferenceElementRefiner : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_ReferenceElementRefiner const* object(
                                         std::string const& a_name ) ;

   //-- Characteristics

      GE_ReferencePolyhedronRefiner const* cell_refiner( void ) const ;

      PDE_ReferenceElement const* reference_element( void ) const ;

      enum OneLevelDiffRule{ Supports, Parents, No, Invalid } ;

   //-- Refinement equation

      // number of childs of the coarse node `a_parent_node'
      // in the `ic'-th subcell
      size_t nb_childs( size_t a_parent_node, size_t ic ) const ;

      size_t child_node( size_t a_parent_node, size_t ic, size_t icn ) const ;

      // number of parent nodes for the fine node `a_child_node'
      // of the `ic'-th subcell
      size_t nb_parents( size_t a_child_node, size_t ic ) const ;

      size_t parent_node( size_t a_child_node, size_t ic, size_t ipn ) const ;

      //One level difference rule
      static void set_one_level_difference_rule( OneLevelDiffRule a_rule ) ;

      static OneLevelDiffRule one_level_difference_rule( void ) ;

      //number of refinement dep node for the fine node `a_fine_node'
      //of the `ic'-th subcell
      size_t nb_ascendants( size_t a_fine_node, size_t ic ) const ;

      size_t ascendant_node( size_t a_fine_node,
                             size_t ic, size_t iun ) const ;

      size_t nb_descendants( size_t a_coarse_node, size_t ic ) const ;

      size_t descendant_node( size_t a_coarse_node,
                              size_t ic, size_t iun ) const ;

      // index (according to the numbering of `::reference_element')
      // of the leading child node (defined, for Lagrange elements,
      // as the child node with the same geometric location than its parent,
      // if any)
      size_t leading_child_node( size_t a_parent_node, size_t ic ) const ;

      double refinement_coef( size_t a_parent_node,
                              size_t a_child_node, size_t ic ) const ;

   protected: //--------------------------------------------------------

      virtual ~PDE_ReferenceElementRefiner( void ) ;

      PDE_ReferenceElementRefiner(
                   std::string const& name_of_reference_element,
                   GE_ReferencePolyhedronRefiner const* rcell_refiner ) ;

      void append_child( size_t a_parent_node,
                         size_t a_subcell,
                         size_t a_child_node,
                         double a_child_coef,
                         bool is_leading_child ) ;

   private: //----------------------------------------------------------

      PDE_ReferenceElementRefiner( void ) ;
      PDE_ReferenceElementRefiner( PDE_ReferenceElementRefiner const& other ) ;
      PDE_ReferenceElementRefiner& operator=(
                                   PDE_ReferenceElementRefiner const& other ) ;

   //-- Class attributes

      static OneLevelDiffRule OLRULE ;

   //-- Attributes

      GE_ReferencePolyhedronRefiner const* CELL_REFINER ;
      PDE_ReferenceElement const* ELM ;

      size_t_array2D NB_PARENTS ;
      size_t_array2D NB_CHILDS ;

      size_t_array3D PARENT ;
      size_t_array3D CHILD ;

      size_t_array2D LEADING_CHILD ;
      doubleArray3D  REFI_COEF ;
} ;

#endif
