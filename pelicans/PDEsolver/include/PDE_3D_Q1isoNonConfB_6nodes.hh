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

#ifndef PDE_3D_Q1_ISO_NONCONF_B_6_NODES_HH
#define PDE_3D_Q1_ISO_NONCONF_B_6_NODES_HH

#include <PDE_ReferenceElement.hh>

/*
Isoparametric nonconforming reference finite element with the following 
characteristics :
   Dimension            : 3
   Reference polyhedron : `GE_ReferenceCube::'
   Polynomial degree    : 1
   number of nodes      : 6
   Nodes are located at the side centers.

The nodal values refer to the value of the field at the center of the side
containing the node.

Ref. Rannacher and Turek (1992) "A simple nonconforming quadrilateral
Stokes element", Numer. Meth. Part. Diff. Equ., n�8, pp 97-111.

PUBLISHED
*/

class PEL_EXPORT PDE_3D_Q1isoNonConfB_6nodes : public PDE_ReferenceElement
{

   public: //-----------------------------------------------------------

   //-- Basis functions

      virtual double N_local( size_t node,
                              GE_Point const* pt_ref ) const ;

      virtual double dN_local( size_t node,
                               size_t a,
                               GE_Point const* pt_ref ) const ;

      virtual double d2N_local( size_t node,
                                size_t a, 
                                size_t b,
                                GE_Point const* pt_ref ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_3D_Q1isoNonConfB_6nodes( void ) ;
     ~PDE_3D_Q1isoNonConfB_6nodes( void ) ;
      PDE_3D_Q1isoNonConfB_6nodes( PDE_3D_Q1isoNonConfB_6nodes const& other ) ;
      PDE_3D_Q1isoNonConfB_6nodes& operator=( 
                                   PDE_3D_Q1isoNonConfB_6nodes const& other ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool reference_polyhedron_POST(
                          GE_ReferencePolyhedron const* result ) const ;

   //-- Class attributes

      static PDE_3D_Q1isoNonConfB_6nodes const* REGISTRATOR ;

} ;

#endif