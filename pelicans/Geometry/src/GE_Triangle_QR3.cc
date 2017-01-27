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

#include <GE_Triangle_QR3.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceTriangle.hh>
#include <GE_RefinedTriangle_QR.hh>

GE_Triangle_QR3 const*
GE_Triangle_QR3:: REGISTRATOR = unique_instance() ;

//----------------------------------------------------------------------
GE_Triangle_QR3:: GE_Triangle_QR3( void )
//----------------------------------------------------------------------
   : GE_QuadratureRule( "GE_Triangle_QR3", GE_ReferenceTriangle::object(), 3 )
{
   append_point( GE_Point::create( this, 1./3., 1./3. ), -27./96. ) ;
   append_point( GE_Point::create( this, 1./5., 1./5. ),  25./96. ) ;
   append_point( GE_Point::create( this, 3./5., 1./5. ),  25./96. ) ;
   append_point( GE_Point::create( this, 1./5., 3./5. ),  25./96. ) ;

   set_sum_of_weights( 0.5 ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_Triangle_QR3:: ~GE_Triangle_QR3( void )
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;

   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
GE_Triangle_QR3 const*
GE_Triangle_QR3:: unique_instance( void )
//----------------------------------------------------------------------
{
   static GE_Triangle_QR3 const* rule = 0 ;
   if( rule == 0 )
   {
      rule = new GE_Triangle_QR3() ;

      GE_RefinedTriangle_QR const* rrule = 
                 GE_RefinedTriangle_QR::create( "GE_Triangle4R_QR3", rule ) ;
      
      rrule = 0 ;
   }
   return( rule ) ;
}

//-----------------------------------------------------------------------------
bool
GE_Triangle_QR3:: invariant( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Triangle_QR3:: invariant" ) ;
   PEL_ASSERT( GE_QuadratureRule::invariant() ) ;
   PEL_ASSERT( IMPLIES(REGISTRATOR!=0,GE_QuadratureRule::nb_points()==4) ) ;

   return( true ) ;
}
