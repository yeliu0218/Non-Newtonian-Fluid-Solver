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

#ifndef GE_REFINED_TETRAHEDRON_QR_HH
#define GE_REFINED_TETRAHEDRON_QR_HH

#include <GE_QuadratureRule.hh>

class GE_ReferencePolyhedron ;

/*
Quadrature rules constructed by gathering quadrature rules on 
the 8 subtetrahedra of the P1isoP2 reference element.

These 8 subtetrahedra are obtained via a subdivision of the unit tetrahedron
obtained as follows: 
   * each face is divided into four new triangles by connecting the midpoint
     of the edges
   * the midpoints (1/2,0,0) and (0,1/2,1/2) (which belong to a pair 
     of non-intersecting edges) are connected.

PUBLISHED
*/

class PEL_EXPORT GE_RefinedTetrahedron_QR : public GE_QuadratureRule
{
   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance defined by setting a quadrature rule 
      // based on `tria_rule' on each of the 8 subtetrahedra 
      // subdividing the reference tetrahedron.
      static GE_RefinedTetrahedron_QR const* create( 
                                      std::string a_name,
                                      GE_QuadratureRule const* tria_rule ) ;
      
   protected: //---------------------------------------------------------------

   private: //-----------------------------------------------------------------
     
      GE_RefinedTetrahedron_QR( void ) ;
     ~GE_RefinedTetrahedron_QR( void ) ;
      GE_RefinedTetrahedron_QR( GE_RefinedTetrahedron_QR const& other ) ;
      GE_RefinedTetrahedron_QR& operator=( GE_RefinedTetrahedron_QR const& other ) ;

      GE_RefinedTetrahedron_QR( std::string a_name,
                                GE_QuadratureRule const* tria_rule ) ;
} ;

#endif
