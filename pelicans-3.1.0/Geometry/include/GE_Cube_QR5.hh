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

#ifndef GE_CUBE_QR5_HH
#define GE_CUBE_QR5_HH

#include <GE_QuadratureRule.hh>

/*
Quadrature rule (in the Gaussian sense) with the following characteristics :
   registration name : "GE_Cube_QR5"
   polyhedron        : `GE_ReferenceCube::'
   number of points  : 14
   precision         : exact for polynomials of degree up to 5

Reference : G. Dhatt and G. Touzot, 
"Une presentation de la methode des elements finis", chapter 5.1
Maloine, PARIS, 1989.

This class has a unique instance, whose creation is immediately followed
by the creation of an additional quadrature rule of type `GE_RefinedCube_QR::'
with the following characteristics: 
   registration name : "GE_Cube8R_QR5"
   polyhedron        : `GE_ReferenceCube::'
   number of points  : 8*14
   precision         : exact for functions whose restrictions to each of the
                       8 subcubes are polynomials of degree up to 5

PUBLISHED 
*/

class PEL_EXPORT GE_Cube_QR5 : public GE_QuadratureRule
{
   public: //------------------------------------------------------------------

   protected: //---------------------------------------------------------------

   private: //-----------------------------------------------------------------

      GE_Cube_QR5( void ) ;
     ~GE_Cube_QR5( void ) ;
      GE_Cube_QR5( GE_Cube_QR5 const& other ) ;
      GE_Cube_QR5& operator=( GE_Cube_QR5 const& other ) ;

      static GE_Cube_QR5 const* unique_instance( void ) ;
      
    //--  Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

    //-- Class Attributes

      static GE_Cube_QR5 const* REGISTRATOR ;
} ;

#endif
