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

#ifndef PDE_REFERENCE_ELEMENT_HH
#define PDE_REFERENCE_ELEMENT_HH

#include <PEL_Object.hh>

#include <string>

class GE_Mpolyhedron ;
class GE_Point ;
class GE_ReferencePolyhedron ;

class PEL_ObjectRegister ;
class PEL_Vector ;

/*
Lagrange reference finite elements.

A Lagrange finite element of R(n) is a triple (K,P,S) given by :
   1. A closed subset K of R(n) with a nonempty interior and a
      Lipschitz-continuous boundary.
   2. A finite dimensional space P of real valued functions defined over K
      (let N denote the dimension of P).
   3. A finite set S of N distinct points of K (called the set of nodes)
      that is assumed, by definition, to be P-unisolvent.

S being P-unisolvent, the exist N functions pi (i=1,...,N) in P
such that for any aj (j=1,...,N) in P,
   pi(aj) = 1 if i=j
   pi(aj) = 0 otherwise
These N functions pi are called basis functions or shape functions.

For all descendants of PDE_ReferenceElement, the basis functions are
polynomials.

Two finite elements (K1,P1,S1) and (K2,P2,S2) are said to be equivalent
provided  that there exist a one-to-one mapping F such that :
   1. F(K1) = K2
   2. any function of P2 is the composition of a function of P1 with the
      inverse of F ;
   3. any node of S2 is the image of a node of S1.
Note : K1 and K2 might be respectively subsets of R(n) and R(m),
with n and m non identical.

From one reference finite element (Kref,Pref,Sref), a whole family of finite
elements that are equivalent to (Kref, Pref, Sref) can be constructed
(through mappings that are implemented in `GE_Mpolyhderon' instances).

FRAMEWORK INSTANTIATION
*/

class PEL_EXPORT PDE_ReferenceElement : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_ReferenceElement const* object( std::string a_name ) ;

      static PDE_ReferenceElement const* object_with_nodes_at_vertices(
	                                         GE_Mpolyhedron const* poly ) ;

      static size_t nb_objects( void ) ;

   //-- Identification

      size_t id_number( void ) const ;

      virtual std::string const& name( void ) const ;

   //-- Geometrical characteristics

      // number of space dimensions (note that equivalent elements built
      // from `self' through mappings might belong to spaces of dimension
      // bigger than the returned result)
      size_t dimension( void ) const ;

      // geometrical support
      GE_ReferencePolyhedron const* reference_polyhedron( void ) const ;

   //-- Nodes

      // number of nodes
      size_t nb_nodes( void ) const ;

      // Is there a node located at point `pt_ref' ?
      bool has_node( GE_Point const* pt_ref ) const ;

      // index of the node located at `pt_ref'
      size_t local_node( GE_Point const* pt_ref ) const ;

      // location of the node of index `node'
      GE_Point const* node_location( size_t node ) const ;


   //-- Basis functions

      // value the `node'-th basis function at point `pt_ref'
      virtual double N_local( size_t node, GE_Point const* pt_ref ) const = 0 ;

      // value at point `pt_ref' of the derivative with respect to the `a'-th
      // direction of the `node'-th basis functions
      virtual double dN_local( size_t node,
                               size_t a,
                               GE_Point const* pt_ref ) const = 0 ;

      // value at point `pt_ref' of the double derivative with respect to
      // the `a'-th and `b'-th directions of the `node'-th basis functions
      virtual double d2N_local( size_t node,
                                size_t a,
                                size_t b,
                                GE_Point const* pt_ref ) const = 0 ;

   //-- Comparison

      virtual bool comparable( PEL_Object const* other ) const ;

      virtual bool is_equal( PEL_Object const* other ) const ;

      virtual int three_way_comparison( PEL_Object const* other ) const ;

      virtual size_t hash_code( void ) const ;

   protected: //--------------------------------------------------------

      virtual ~PDE_ReferenceElement( void ) ;

      PDE_ReferenceElement( std::string a_name,
                            GE_ReferencePolyhedron const* a_ref_poly ) ;

      void append_node( GE_Point* node_point ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool N_local_PRE( size_t node, GE_Point const* pt_ref ) const ;

      virtual bool dN_local_PRE( size_t node,
                                 size_t a,
                                 GE_Point const* pt_ref ) const ;

      virtual bool d2N_local_PRE( size_t node,
                                  size_t a,
                                  size_t b,
                                  GE_Point const* pt_ref ) const ;

      virtual bool reference_polyhedron_POST(
                               GE_ReferencePolyhedron const* result ) const ;

      virtual bool invariant( void ) const ;

   private: //----------------------------------------------------------

      PDE_ReferenceElement( void ) ;
      PDE_ReferenceElement( PDE_ReferenceElement const& other ) ;
      PDE_ReferenceElement const& operator=(
                            PDE_ReferenceElement const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Class attributes

      static size_t NB_INSTANCES ;

   //-- Attributes

      std::string ONAME ;
      size_t ID ;
      
      size_t NB_BFS ;
      PEL_Vector* NODE_PTS ;
      GE_ReferencePolyhedron const* REF_POLY ;
} ;

#endif
