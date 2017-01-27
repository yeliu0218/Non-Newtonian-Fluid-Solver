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

#ifndef GE_Segment_QR9_HH
#define GE_Segment_QR9_HH

#include <GE_QuadratureRule.hh>

/*
Quadrature rule (in the Gaussian sense) with the following characteristics:
   registration name : "GE_Segment_QR9"
   polyhedron        : `GE_ReferenceSegment::'
   number of points  : 5
   precision         : exact for polynomials of degree up to 9

Reference : G. Dhatt and G. Touzot, 
"Une presentation de la methode des elements finis", chapter 5.1,
Maloine, PARIS, 1989.

This class has a unique instance, whose creation is immediately followed
by the creation of two product quadrature rules of type `GE_Product_QR::'
with the characteristics given below :

product rule on the unit square:
   registration name : "GE_Segment_QR9_to_Square"
   polyhedron        : `GE_ReferenceSquare::'
   number of points  : 5*5
   precision         : exact for polynomials of degree up to 9
                       component by component

product rule on the unit cube:
   registration name : "GE_Segment_QR9_to_Cube"
   polyhedron        : `GE_ReferenceCube::'
   number of points  : 5*5*5
   precision         : exact for polynomials of degree up to 9
                       component by component

The creation of the product quadrature rule on the unit square is followed by
the creation an additional quadrature rule of type `GE_RefinedSquare_QR::'
with the following characteristics:
   registration name : "GE_Segment_QR9_to_Square4R"
   polyhedron        : `GE_ReferenceSquare::'
   number of points  : 4*(5*5)
   precision         : exact for functions whose restrictions to each of the
                       4 subsquares are polynomials of degree up to 9
                       component by component

The creation of the product quadrature rule on the unit cube is followed by
the creation an additional quadrature rule of type `GE_RefinedCube_QR::'
with the following characteristics:
   registration name : "GE_Segment_QR9_to_Cube8R"
   polyhedron        : `GE_ReferenceCube::'
   number of points  : 8*(5*5*5)
   precision         : exact for functions whose restrictions to each of the
                       8 subcubes are polynomials of degree up to 9
                       component by component
   
PUBLISHED */

class PEL_EXPORT GE_Segment_QR9 : public GE_QuadratureRule
{
   public: //-----------------------------------------------------------

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      GE_Segment_QR9( void ) ;
     ~GE_Segment_QR9( void ) ;
      GE_Segment_QR9( GE_Segment_QR9 const& other ) ;
      GE_Segment_QR9& operator=( GE_Segment_QR9 const& other ) ;

      static GE_Segment_QR9 const* unique_instance( void ) ;

    //--  Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

    //-- Class Attributes

      static GE_Segment_QR9 const* REGISTRATOR ;

} ;

#endif
