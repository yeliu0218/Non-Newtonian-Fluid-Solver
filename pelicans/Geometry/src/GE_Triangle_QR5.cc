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

#include <GE_Triangle_QR5.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceTriangle.hh>
#include <GE_RefinedTriangle_QR.hh>

GE_Triangle_QR5 const* 
GE_Triangle_QR5:: REGISTRATOR = unique_instance() ;

//----------------------------------------------------------------------
GE_Triangle_QR5:: GE_Triangle_QR5( void )
//----------------------------------------------------------------------
   : GE_QuadratureRule( "GE_Triangle_QR5", GE_ReferenceTriangle::object(), 5 )
{
   double const a = ( 6. + PEL::sqrt( 15. ) )/21 ;
   double const b = 4./7.-a ;
   double const c = 1./3. ;
   double const coef1 = ( 155. + PEL::sqrt( 15. ) )/2400. ;
   double const coef2 = 31./240. - coef1 ;
   double const coef3 = 9./80. ;

   append_point( GE_Point::create( this, c,    c    ), coef3 ) ;
   append_point( GE_Point::create( this, a,    a    ), coef1 ) ;
   append_point( GE_Point::create( this, 1.-2.*a, a ), coef1 ) ;
   append_point( GE_Point::create( this, a , 1-2.*a ), coef1 ) ;
   append_point( GE_Point::create( this, b, b       ), coef2 ) ;
   append_point( GE_Point::create( this, 1.-2.*b, b ), coef2 ) ;
   append_point( GE_Point::create( this, b, 1.-2.*b ), coef2 ) ;

   set_sum_of_weights( 0.5 ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Triangle_QR5:: ~GE_Triangle_QR5( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
GE_Triangle_QR5 const*
GE_Triangle_QR5:: unique_instance( void )
//----------------------------------------------------------------------
{
   static GE_Triangle_QR5 const* rule = 0 ;
   if( rule == 0 )
   {
      rule = new GE_Triangle_QR5() ;

      GE_RefinedTriangle_QR const* rrule = 
                 GE_RefinedTriangle_QR::create( "GE_Triangle4R_QR5", rule ) ;
      
      rrule = 0 ;
   }
   return( rule ) ;
}

//----------------------------------------------------------------------
bool
GE_Triangle_QR5:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle_QR5:: invariant" ) ;
   PEL_ASSERT( GE_QuadratureRule::invariant() ) ;
   PEL_ASSERT( IMPLIES(REGISTRATOR!=0,GE_QuadratureRule::nb_points()==7) ) ;

   return( true ) ;
}
