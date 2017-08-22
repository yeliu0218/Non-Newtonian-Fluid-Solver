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

#include <PDE_3D_Q2_27nodes.hh>

#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceCube.hh>

#include <iostream>

PDE_3D_Q2_27nodes const* 
PDE_3D_Q2_27nodes:: REGISTRATOR = new PDE_3D_Q2_27nodes() ;

//----------------------------------------------------------------------
PDE_3D_Q2_27nodes:: ~PDE_3D_Q2_27nodes( void )
//----------------------------------------------------------------------
{
   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
PDE_3D_Q2_27nodes:: PDE_3D_Q2_27nodes ( void )
//----------------------------------------------------------------------
   : PDE_ReferenceElement( "PDE_3D_Q2_27nodes",
                           GE_ReferenceCube::object() )
{
   append_node( GE_Point::create( this, 0.0, 0.0, 0.0 ) ) ; // 0
   append_node( GE_Point::create( this, 0.5, 0.0, 0.0 ) ) ; // 1
   append_node( GE_Point::create( this, 1.0, 0.0, 0.0 ) ) ; // 2
   append_node( GE_Point::create( this, 1.0, 0.5, 0.0 ) ) ; // 3
   append_node( GE_Point::create( this, 0.5, 0.5, 0.0 ) ) ; // 4
   append_node( GE_Point::create( this, 0.0, 0.5, 0.0 ) ) ; // 5
   append_node( GE_Point::create( this, 0.0, 1.0, 0.0 ) ) ; // 6
   append_node( GE_Point::create( this, 0.5, 1.0, 0.0 ) ) ; // 7
   append_node( GE_Point::create( this, 1.0, 1.0, 0.0 ) ) ; // 8

   append_node( GE_Point::create( this, 0.0, 0.0, 0.5 ) ) ; // 9
   append_node( GE_Point::create( this, 0.5, 0.0, 0.5 ) ) ; // 10
   append_node( GE_Point::create( this, 1.0, 0.0, 0.5 ) ) ; // 11
   append_node( GE_Point::create( this, 1.0, 0.5, 0.5 ) ) ; // 12
   append_node( GE_Point::create( this, 0.5, 0.5, 0.5 ) ) ; // 13
   append_node( GE_Point::create( this, 0.0, 0.5, 0.5 ) ) ; // 14
   append_node( GE_Point::create( this, 0.0, 1.0, 0.5 ) ) ; // 15
   append_node( GE_Point::create( this, 0.5, 1.0, 0.5 ) ) ; // 16
   append_node( GE_Point::create( this, 1.0, 1.0, 0.5 ) ) ; // 17

   append_node( GE_Point::create( this, 0.0, 0.0, 1.0 ) ) ; // 18 
   append_node( GE_Point::create( this, 0.5, 0.0, 1.0 ) ) ; // 19
   append_node( GE_Point::create( this, 1.0, 0.0, 1.0 ) ) ; // 20
   append_node( GE_Point::create( this, 1.0, 0.5, 1.0 ) ) ; // 21
   append_node( GE_Point::create( this, 0.5, 0.5, 1.0 ) ) ; // 22
   append_node( GE_Point::create( this, 0.0, 0.5, 1.0 ) ) ; // 23
   append_node( GE_Point::create( this, 0.0, 1.0, 1.0 ) ) ; // 24
   append_node( GE_Point::create( this, 0.5, 1.0, 1.0 ) ) ; // 25
   append_node( GE_Point::create( this, 1.0, 1.0, 1.0 ) ) ; // 26
}

//----------------------------------------------------------------------
double
PDE_3D_Q2_27nodes:: N_local( size_t node,
                             GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_Q2_27nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;
   
   double result = PEL::max_double() ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   switch( node )
   {
      case 0 :
         result = 8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 1 :
         result = 16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 2 :
         result = -8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 3 :
         result = -16.0 * x*(0.5-x) * y*(1.0-y) * (1.0-z)*(0.5-z);
         break ;
      case 4 :
         result = 32.0 * x*(1.0-x) * y*(1.0-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 5 :
         result = 16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 6 :
         result =-8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 7 :
         result = -16.0 * x*(1.0-x) * y*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 8 :
         result = 8.0 * x*(0.5-x) * y*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 9 :
         result = 16.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * z*(1.0-z) ;
         break ;
      case 10 :
         result = 32.0 * x*(1.0-x) * (1.0-y)*(0.5-y) *  z*(1.0-z);
         break ;
      case 11 :
         result = -16.0 * x*(0.5-x) * (1.0-y)*(0.5-y) *  z*(1.0-z) ;
         break ;
      case 12 :
         result = -32.0 * x*(0.5-x) * y*(1.0-y) * z*(1.0-z) ;
         break ;
      case 13 :
         result = 64.0 * x*(1.0-x) * y*(1.0-y) *  z*(1.0-z) ;
         break ;
      case 14 :
         result = 32.0 * (1.0-x)*(0.5-x) * y*(1.0-y) *  z*(1.0-z) ;
         break ;
      case 15 :
         result =-16.0 * (1.0-x)*(0.5-x) * y*(0.5-y) *  z*(1.0-z) ;
         break ;
      case 16 :
         result = -32.0 * x*(1.0-x) * y*(0.5-y) *  z*(1.0-z) ;
         break ;
      case 17 :
         result = 16.0 * x*(0.5-x) * y*(0.5-y) * z*(1.0-z);
         break ;
      case 18 :
         result = -8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * z*(0.5-z) ;
         break ;
      case 19 :
         result = -16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) *  z*(0.5-z);
         break ;
      case 20 :
         result = 8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) *  z*(0.5-z) ;
         break ;
      case 21 :
         result = 16.0 * x*(0.5-x) * y*(1.0-y) * z*(0.5-z) ;
         break ;
      case 22 :
         result = -32.0 * x*(1.0-x) * y*(1.0-y) *  z*(0.5-z) ;
         break ;
      case 23 :
         result = -16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) *  z*(0.5-z) ;
         break ;
      case 24 :
         result = 8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) *  z*(0.5-z) ;
         break ;
      case 25 :
         result = 16.0 * x*(1.0-x) * y*(0.5-y) *  z*(0.5-z) ;
         break ;
      case 26 :
         result = -8.0 * x*(0.5-x) * y*(0.5-y) * z*(0.5-z);
         break ; 
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_3D_Q2_27nodes:: dN_local( size_t node, size_t a,
                              GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_Q2_27nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double result = PEL::max_double() ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   switch( node )
   {
      case 0 :
         switch( a )
         {
            case 0 :
               result = 8.0 * (2.0*x-1.5) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z);
               break ;
            case 1 :
               result = 8.0 * (1.0-x)*(0.5-x) * (2.0*y-1.5)  * (1.0-z)*(0.5-z);
               break ;
            case 2 :
               result = 8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * (2.0*z-1.5) ;
               break ;
         }
         break ;
      case 1 :
         switch( a )
         {
            case 0 :
               result = 16.0 * (-2.0*x+1.0) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z);
               break ;
            case 1 :
               result = 16.0 * x*(1.0-x) * (2.0*y-1.5) * (1.0-z)*(0.5-z);
               break ;
            case 2 :
               result = 16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) * (2.0*z-1.5) ;
               break ;
         }
         break ;
      case 2 :
         switch( a )
         {
            case 0 :
               result = -8.0 * (-2.0*x+0.5) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z);
               break ;
            case 1 :
               result = -8.0 * x*(0.5-x) * (2.0*y-1.5) * (1.0-z)*(0.5-z) ;
               break ;
            case 2 :
               result = -8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) * (2.0*z-1.5) ;
               break ;
         }
         break ;
      case 3 :
         switch( a )
         {
            case 0 :
               result = -16.0 * (-2.0*x+0.5) * y*(1.0-y) * (1.0-z)*(0.5-z) ;
               break ;
            case 1 :
               result = -16.0 * x*(0.5-x) * (-2.0*y+1.0) * (1.0-z)*(0.5-z);
               break ;
            case 2 :
               result = -16.0 * x*(0.5-x) * y*(1.0-y) * (2.0*z-1.5);
               break ;
         }
         break ;
      case 4 :
         switch( a )
         {
            case 0 :
               result = 32.0 * (-2.0*x+1.0) * y*(1.0-y) * (1.0-z)*(0.5-z);
               break ;
            case 1 :
               result = 32.0 * x*(1.0-x) * (-2.0*y+1.0) * (1.0-z)*(0.5-z);
               break ;
            case 2 :
               result = 32.0 * x*(1.0-x) * y*(1.0-y) * (2.0*z-1.5);
               break ;
         }
         break ;
      case 5 :
         switch( a )
         {
            case 0 :
               result = 16.0 * (2.0*x-1.5) * y*(1.0-y) * (1.0-z)*(0.5-z);
               break ;
            case 1 :
               result = 16.0 * (1.0-x)*(0.5-x) * (-2.0*y+1.0) * (1.0-z)*(0.5-z);
               break ;
            case 2 :
               result = 16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) * (2.0*z-1.5);
               break ;
         }
         break ;
      case 6 :
         switch( a )
         {
            case 0 :
               result = -8.0 * (2.0*x-1.5) * y*(0.5-y) * (1.0-z)*(0.5-z);
               break ;
            case 1 :
               result = -8.0 * (1.0-x)*(0.5-x) * (-2.0*y+0.5) * (1.0-z)*(0.5-z);
               break ;
            case 2 :
               result = -8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) * (2.0*z-1.5);
               break ;
         }
         break ;
      case 7 :
         switch( a )
         {
            case 0 :
               result = -16.0 * (-2.0*x+1.0) * y*(0.5-y) * (1.0-z)*(0.5-z);
               break ;
            case 1 :
               result = -16.0 * x*(1.0-x) * (-2.0*y+0.5) * (1.0-z)*(0.5-z);
               break ;
            case 2 :
               result = -16.0 * x*(1.0-x) * y*(0.5-y) * (2.0*z-1.5);
               break ;
         }
         break ;
      case 8 :
         switch( a )
         {
            case 0 :
               result = 8.0 * (-2.0*x+0.5) * y*(0.5-y) * (1.0-z)*(0.5-z);
               break ;
            case 1 :
               result = 8.0 * x*(0.5-x) * (-2.0*y+0.5) * (1.0-z)*(0.5-z);
               break ;
            case 2 :
               result = 8.0 * x*(0.5-x) * y*(0.5-y) * (2.0*z-1.5);
               break ;
         }
         break ;
      case 9 :
         switch( a )
         {
            case 0 :
               result = 16.0 * (2.0*x-1.5) * (1.0-y)*(0.5-y) * z*(1.0-z);
               break ;
            case 1 :
               result = 16.0 * (1.0-x)*(0.5-x) * (2.0*y-1.5)  * z*(1.0-z);
               break ;
            case 2 :
               result = 16.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * (-2.0*z+1.0) ;
               break ;
         }
         break ;
      case 10 :
         switch( a )
         {
            case 0 :
               result = 32.0 * (-2.0*x+1.0) * (1.0-y)*(0.5-y) * z*(1.0-z);
               break ;
            case 1 :
               result = 32.0 * x*(1.0-x) * (2.0*y-1.5) * z*(1.0-z);
               break ;
            case 2 :
               result = 32.0 * x*(1.0-x) * (1.0-y)*(0.5-y) * (-2.0*z+1.0)  ;
               break ;
         }
         break ;
      case 11 :
         switch( a )
         {
            case 0 :
               result = -16.0 * (-2.0*x+0.5) * (1.0-y)*(0.5-y) * z*(1.0-z) ;
               break ;
            case 1 :
               result = -16.0 * x*(0.5-x) * (2.0*y-1.5) * z*(1.0-z) ;
               break ;
            case 2 :
               result = -16.0 * x*(0.5-x) * (1.0-y)*(0.5-y) *  (-2.0*z+1.0) ;
               break ;
         }
         break ;
      case 12 :
         switch( a )
         {
            case 0 :
               result = -32.0 * (-2.0*x+0.5) * y*(1.0-y) * z*(1.0-z) ;
               break ;
            case 1 :
               result = -32.0 * x*(0.5-x) * (-2.0*y+1.0) * z*(1.0-z);
               break ;
            case 2 :
               result = -32.0 * x*(0.5-x) * y*(1.0-y) * (-2.0*z+1.0) ;
               break ;
         }
         break ;
      case 13:
         switch( a )
         {
            case 0 :
               result = 64.0 * (-2.0*x+1.0) * y*(1.0-y) * z*(1.0-z) ;
               break ;
            case 1 :
               result = 64.0 * x*(1.0-x) * (-2.0*y+1.0) * z*(1.0-z);
               break ;
            case 2 :
               result = 64.0 * x*(1.0-x) * y*(1.0-y) *  (-2.0*z+1.0);
               break ;
         }
         break ;
      case 14 :
         switch( a )
         {
            case 0 :
               result = 32.0 * (2.0*x-1.5) * y*(1.0-y) * z*(1.0-z);
               break ;
            case 1 :
               result = 32.0 * (1.0-x)*(0.5-x) * (-2.0*y+1.0) * z*(1.0-z);
               break ;
            case 2 :
               result = 32.0 * (1.0-x)*(0.5-x) * y*(1.0-y) * (-2.0*z+1.0) ;
               break ;
         }
         break ;
      case 15 :
         switch( a )
         {
            case 0 :
               result = -16.0 * (2.0*x-1.5) * y*(0.5-y) * z*(1.0-z);
               break ;
            case 1 :
               result = -16.0 * (1.0-x)*(0.5-x) * (-2.0*y+0.5) * z*(1.0-z);
               break ;
            case 2 :
               result = -16.0 * (1.0-x)*(0.5-x) * y*(0.5-y) * (-2.0*z+1.0) ;
               break ;
         }
         break ;
      case 16 :
         switch( a )
         {
            case 0 :
               result = -32.0 * (-2.0*x+1.0) * y*(0.5-y) * z*(1.0-z);
               break ;
            case 1 :
               result = -32.0 * x*(1.0-x) * (-2.0*y+0.5) * z*(1.0-z);
               break ;
            case 2 :
               result = -32.0 * x*(1.0-x) * y*(0.5-y) *  (-2.0*z+1.0);
               break ;
         }
         break ;
      case 17 :
         switch( a )
         {
            case 0 :
               result = 16.0 * (-2.0*x+0.5) * y*(0.5-y) * z*(1.0-z);
               break ;
            case 1 :
               result = 16.0 * x*(0.5-x) * (-2.0*y+0.5) * z*(1.0-z);
               break ;
            case 2 :
               result = 16.0 * x*(0.5-x) * y*(0.5-y) *  (-2.0*z+1.0) ;
               break ;
         }
         break ;
      case 18 :
         switch( a )
         {
            case 0 :
               result = -8.0 * (2.0*x-1.5) * (1.0-y)*(0.5-y) * z*(0.5-z);
               break ;
            case 1 :
               result = -8.0 * (1.0-x)*(0.5-x) * (2.0*y-1.5)  * z*(0.5-z);
               break ;
            case 2 :
               result = -8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * (-2.0*z+0.5) ;
               break ;
         }
         break ;
      case 19 :
         switch( a )
         {
            case 0 :
               result = -16.0 * (-2.0*x+1.0) * (1.0-y)*(0.5-y) * z*(0.5-z);
               break ;
            case 1 :
               result = -16.0 * x*(1.0-x) * (2.0*y-1.5) * z*(0.5-z);
               break ;
            case 2 :
               result = -16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) * (-2.0*z+0.5) ;
               break ;
         }
         break ;
      case 20 :
         switch( a )
         {
            case 0 :
               result = 8.0 * (-2.0*x+0.5) * (1.0-y)*(0.5-y) * z*(0.5-z);
               break ;
            case 1 :
               result = 8.0 * x*(0.5-x) * (2.0*y-1.5) * z*(0.5-z) ;
               break ;
            case 2 :
               result = 8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) *  (-2.0*z+0.5);
               break ;
         }
         break ;
      case 21 :
         switch( a )
         {
            case 0 :
               result = 16.0 * (-2.0*x+0.5) * y*(1.0-y) * z*(0.5-z) ;
               break ;
            case 1 :
               result = 16.0 * x*(0.5-x) * (-2.0*y+1.0) * z*(0.5-z) ;
               break ;
            case 2 :
               result = 16.0 * x*(0.5-x) * y*(1.0-y) * (-2.0*z+0.5);
               break ;
         }
         break ;
      case 22 :
         switch( a )
         {
            case 0 :
               result = -32.0 * (-2.0*x+1.0) * y*(1.0-y) * z*(0.5-z);
               break ;
            case 1 :
               result = -32.0 * x*(1.0-x) * (-2.0*y+1.0) * z*(0.5-z) ;
               break ;
            case 2 :
               result = -32.0 * x*(1.0-x) * y*(1.0-y) * (-2.0*z+0.5);
               break ;
         }
         break ;
      case 23 :
         switch( a )
         {
            case 0 :
               result = -16.0 * (2.0*x-1.5) * y*(1.0-y) * z*(0.5-z);
               break ;
            case 1 :
               result = -16.0 * (1.0-x)*(0.5-x) * (-2.0*y+1.0) * z*(0.5-z);
               break ;
            case 2 :
               result = -16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) * (-2.0*z+0.5) ;
               break ;
         }
         break ;
      case 24 :
         switch( a )
         {
            case 0 :
               result = 8.0 * (2.0*x-1.5) * y*(0.5-y) * z*(0.5-z);
               break ;
            case 1 :
               result = 8.0 * (1.0-x)*(0.5-x) * (-2.0*y+0.5) * z*(0.5-z);
               break ;
            case 2 :
               result = 8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) * (-2.0*z+0.5) ;
               break ;
         }
         break ;
      case 25 :
         switch( a )
         {
            case 0 :
               result = 16.0 * (-2.0*x+1.0) * y*(0.5-y) * z*(0.5-z);
               break ;
            case 1 :
               result = 16.0 * x*(1.0-x) * (-2.0*y+0.5) * z*(0.5-z);
               break ;
            case 2 :
               result = 16.0 * x*(1.0-x) * y*(0.5-y) * (-2.0*z+0.5)  ;
               break ;
         }
         break ;
      case 26 :
         switch( a )
         {
            case 0 :
               result = -8.0 * (-2.0*x+0.5) * y*(0.5-y) * z*(0.5-z) ;
               break ;
            case 1 :
               result = -8.0 * x*(0.5-x) * (-2.0*y+0.5) * z*(0.5-z);
               break ;
            case 2 :
               result = -8.0 * x*(0.5-x) * y*(0.5-y) * (-2.0*z+0.5);
               break ;
         }
         break ;
   }
   return( result ) ;
}


//----------------------------------------------------------------------
double
PDE_3D_Q2_27nodes:: d2N_local( size_t node,
                               size_t a, size_t b,
                               GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_Q2_27nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double result = PEL::max_double() ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   switch( node )
   {
      case 0 : // <--------
         switch( a )
         {
            case 0 :
               switch( b )
               {
                  case 0 :
                     result = 8.0 * 2.0 * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z) ;
                     break ;
                  case 1 :
                     result = 8.0 * (2.0*x-1.5) * (2.0*y-1.5) * (1.0-z)*(0.5-z) ;
                     break ;
                  case 2 :
                     result = 8.0 * (2.0*x-1.5) * (1.0-y)*(0.5-y) * (2.0*z-1.5) ;
                     break ;
               }
               break ;
            case 1 :
               switch( b )
               {
                  case 0 :
                     result = 8.0 * (2.0*x-1.5)  * (2.0*y-1.5)  * (1.0-z)*(0.5-z)  ;
                     break ;
                  case 1 :
                     result = 8.0 * (1.0-x)*(0.5-x) * 2.0 * (1.0-z)*(0.5-z) ;
                     break ;
                  case 2 :
                     result = 8.0 * (1.0-x)*(0.5-x) * (2.0*y-1.5)  * (2.0*z-1.5) ;
                     break ;
               }
               break ;
            case 2 :
               switch( b )
               {
                  case 0 :
                     result =  8.0 * (2.0*x-1.5) * (1.0-y)*(0.5-y) * (2.0*z-1.5)  ;
                     break ;
                  case 1 :
                     result = 8.0 * (1.0-x)*(0.5-x) * (2.0*y-1.5) * (2.0*z-1.5) ;
                     break ;
                  case 2 :
                     result = 8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * 2.0 ;
                     break ;
               }
               break ;
         }
         break ;
      case 1 : // <--------
         switch( a )
         {
            case 0 :
               switch( b )
               {
                  case 0 :
                     result = 16.0 * -2.0 * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z)  ;
                     break ;
                  case 1 :
                     result = 16.0 * (-2.0*x+1.0) * (2.0*y-1.5) * (1.0-z)*(0.5-z) ;
                     break ;
                  case 2 :
                     result = 16.0 * (-2.0*x+1.0) * (1.0-y)*(0.5-y) * (2.0*z-1.5)  ;
                     break ;
               }
               break ;
            case 1 :
               switch( b )
               {
                  case 0 :
                     result =  16.0 * (-2.0*x+1.0) * (2.0*y-1.5) * (1.0-z)*(0.5-z) ;
                     break ;
                  case 1 :
                     result = 16.0 * x*(1.0-x) * 2.0 * (1.0-z)*(0.5-z) ;
                     break ;
                  case 2 :
                     result = 16.0 * x*(1.0-x) * (2.0*y-1.5) * (2.0*z-1.5) ;
                     break ;
               }
               break ;
            case 2 :
               switch( b )
               {
                  case 0 :
                     result = 16.0 * (-2.0*x+1.0)  * (1.0-y)*(0.5-y) * (2.0*z-1.5)  ;
                     break ;
                  case 1 :
                     result = 16.0 * x*(1.0-x) * (2.0*y-1.5)  * (2.0*z-1.5) ;
                     break ;
                  case 2 :
                     result = 16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) * 2.0 ;
                     break ;
               }
	            break ;
         }
         break ;

      case 2 : // <--------
         switch( a )
         {
            case 0 :
               switch( b )
               {
                  case 0 :
                     result = -8.0 * -2.0 * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z)  ;
                     break ;
                  case 1 :
                     result = -8.0 * (-2.0*x+0.5) * (2.0*y-1.5) * (1.0-z)*(0.5-z)  ;
                     break ;
                  case 2 :
                     result = -8.0 * (-2.0*x+0.5) * (1.0-y)*(0.5-y) * (2.0*z-1.5)  ;
                     break ;
               }
               break ;
	         case 1 :
	            switch( b )
	            {
	               case 0 :
	                  result = -8.0 * (-2.0*x+0.5) * (2.0*y-1.5) * (1.0-z)*(0.5-z);
	                  break ;
	               case 1 :
	                  result = -8.0 * x*(0.5-x) * 2.0 * (1.0-z)*(0.5-z) ;
	                  break ;
	               case 2 :
	                  result = -8.0 * x*(0.5-x) * (2.0*y-1.5) * (2.0*z-1.5) ;
	                  break ;
	            }
	            break ;
	         case 2 :
	            switch( b )
	            {
	               case 0 :
	                  result = -8.0 * (-2.0*x+0.5) * (1.0-y)*(0.5-y) * (2.0*z-1.5) ;
	                  break ;
	               case 1 :
	                  result = -8.0 * x*(0.5-x) * (2.0*y-1.5) * (2.0*z-1.5) ;
	                  break ;
	               case 2 :
	                  result = -8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) * 2.0 ;
	                  break ;		    
	            }
	            break ;
         }
         break ;
      case 3 : // <--------
         switch( a )
         {
            case 0 :
               switch( b )
               {
                  case 0 :
                     result = -16.0 * -2.0 * y*(1.0-y) * (1.0-z)*(0.5-z) ;
                     break ;
                  case 1 :
                     result = -16.0 * (-2.0*x+0.5) * (-2.0*y+1.0) * (1.0-z)*(0.5-z);
                     break ;
                  case 2 :
                     result = -16.0 * (-2.0*x+0.5) * y*(1.0-y) * (2.0*z-1.5);
                     break ;
               }
               break ;
	         case 1 :
	            switch( b )
	            {
	               case 0 :
	                  result = -16.0 * (-2.0*x+0.5) * (-2.0*y+1.0) * (1.0-z)*(0.5-z) ;
	                  break ;
	               case 1 :
	                  result = -16.0 * x*(0.5-x) * -2.0 * (1.0-z)*(0.5-z);
	                  break ;
	               case 2 :
	                  result = -16.0 * x*(0.5-x) * (-2.0*y+1.0) * (2.0*z-1.5);
	                  break ;
	            }
	            break ;
	         case 2 :
	            switch( b )
	            {
	               case 0 :
	                  result = -16.0 * (-2.0*x+0.5) * y*(1.0-y) * (2.0*z-1.5) ;
	                  break ;
	               case 1 :
	                  result = -16.0 * x*(0.5-x) * (-2.0*y+1.0) * (2.0*z-1.5) ;
	                  break ;
	               case 2 :
	                  result = -16.0 * x*(0.5-x) * y*(1.0-y) * 2.0;
	                  break ;	    
	            }
	            break ;
         }
         break ;
      case 4 : // <--------
         switch( a )
         {
            case 0 :
               switch( b )
               {
                  case 0 :
                     result = 32.0 * -2.0 * y*(1.0-y) * (1.0-z)*(0.5-z);
                     break ;
                  case 1 :
                     result = 32.0 * (-2.0*x+1.0) * (-2.0*y+1.0) * (1.0-z)*(0.5-z);
                     break ;
                  case 2 :
                     result = 32.0 * (-2.0*x+1.0) * y*(1.0-y) * (2.0*z-1.5);
                     break ;
               }
               break ;
            case 1 :
               switch( b )
               {
                  case 0 :
                     result = 32.0 * (-2.0*x+1.0) * (-2.0*y+1.0) * (1.0-z)*(0.5-z);
                     break ;
                  case 1 :
                     result = 32.0 * x*(1.0-x) * -2.0 * (1.0-z)*(0.5-z);
                     break ;
                  case 2 :
                     result = 32.0 * x*(1.0-x) * (-2.0*y+1.0) * (2.0*z-1.5);
                     break ;
               }
               break ;
            case 2 :
               switch( b )
               {
                  case 0 :
                     result = 32.0 * (-2.0*x+1.0) * y*(1.0-y) * (2.0*z-1.5) ;
                     break ;
                  case 1 :
                     result = 32.0 * x*(1.0-x) * (-2.0*y+1.0) * (2.0*z-1.5) ;
                     break ;
                  case 2 :
                     result = 32.0 * x*(1.0-x) * y*(1.0-y) * 2.0;
                     break ;	    
               }
               break ;
         }
         break ;
      case 5 : // <--------
         switch( a )
         {
            case 0 :
               switch( b )
               {
                  case 0 :
                     result = 16.0 * 2.0 * y*(1.0-y) * (1.0-z)*(0.5-z);
                     break ;
                  case 1 :
                     result = 16.0 * (2.0*x-1.5) * (-2.0*y+1.0) * (1.0-z)*(0.5-z);
                     break ;
                  case 2 :
                     result = 16.0 * (2.0*x-1.5) * y*(1.0-y) * (2.0*z-1.5);
                     break ;
               }
               break ;
	         case 1 :
	            switch( b )
	            {
	               case 0 :
	                  result = 16.0 * (2.0*x-1.5) * (-2.0*y+1.0) * (1.0-z)*(0.5-z);
	                  break ;
	               case 1 :
	                  result = 16.0 * (1.0-x)*(0.5-x) * -2.0 * (1.0-z)*(0.5-z);
	                  break ;
	               case 2 :
	                  result = 16.0 * (1.0-x)*(0.5-x) * (-2.0*y+1.0)  * (2.0*z-1.5);
	                  break ;
	            }
	            break ;
	         case 2 :
	            switch( b )
	            {
	               case 0 :
	                  result = 16.0 * (2.0*x-1.5) * y*(1.0-y) * (2.0*z-1.5);
	                  break ;
	               case 1 :
	                  result = 16.0 * (1.0-x)*(0.5-x) * (-2.0*y+1.0) * (2.0*z-1.5) ;
	                  break ;
	               case 2 :
	                  result = 16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) * 2.0;
	                  break ;		     
	            }
	            break ;
         }
         break ;
      case 6 : // <--------
         switch( a )
         {
            case 0 :
               switch( b )
               {
                  case 0 :
                     result = -8.0 * 2.0 * y*(0.5-y) * (1.0-z)*(0.5-z);
                     break ;
                  case 1 :
                     result = -8.0 * (2.0*x-1.5) * (-2.0*y+0.5) * (1.0-z)*(0.5-z);
                     break ;
                  case 2 :
                     result = -8.0 * (2.0*x-1.5) * y*(0.5-y) * (2.0*z-1.5);
                     break ;
               }
               break ;
            case 1 :
               switch( b )
               {
                  case 0 :
                     result = -8.0 * (2.0*x-1.5) * (-2.0*y+0.5) * (1.0-z)*(0.5-z);
                     break ;
                  case 1 :
                     result = -8.0 * (1.0-x)*(0.5-x) * -2.0 * (1.0-z)*(0.5-z);
                     break ;
                  case 2 :
                     result = -8.0 * (1.0-x)*(0.5-x) * (-2.0*y+0.5)  * (2.0*z-1.5);
                     break ;   
               }
               break ;
            case 2 :
               switch( b )
               {
                  case 0 :
                     result = -8.0 * (2.0*x-1.5) * y*(0.5-y) * (2.0*z-1.5);
                     break ;
                  case 1 :
                     result = -8.0 * (1.0-x)*(0.5-x) * (-2.0*y+0.5) * (2.0*z-1.5) ;
                     break ;
                  case 2 :
                     result = -8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) * 2.0;
                     break ; 	    
               }
               break ;
         }
         break ;    
      case 7 : // <--------
         switch( a )
         {
            case 0 :
               switch( b )
               {
                  case 0 :
                     result = -16.0 * -2.0 * y*(0.5-y) * (1.0-z)*(0.5-z);
                     break ;
                  case 1 :
                     result = -16.0 * (-2.0*x+1.0) * (-2.0*y+0.5) * (1.0-z)*(0.5-z);
                     break ;
                  case 2 :
                     result = -16.0 * (-2.0*x+1.0) * y*(0.5-y) * (2.0*z-1.5);
                     break ;
               }
               break ;
	         case 1 :
	            switch( b )
	            {
	               case 0 :
	                  result = -16.0 * (-2.0*x+1.0) * (-2.0*y+0.5) * (1.0-z)*(0.5-z);
	                  break ;
	               case 1 :
	                  result = -16.0 * x*(1.0-x) * -2.0 * (1.0-z)*(0.5-z);
	                  break ;
	               case 2 :
	                  result = -16.0 * x*(1.0-x) * (-2.0*y+0.5) * (2.0*z-1.5);
	                  break ;
	            }
	            break ;
	         case 2 :
	            switch( b )
	            {
	               case 0 :
	                  result = -16.0 * (-2.0*x+1.0) * y*(0.5-y) * (2.0*z-1.5) ;
	                  break ;
	               case 1 :
	                  result = -16.0 * x*(1.0-x) * (-2.0*y+0.5) * (2.0*z-1.5) ;
	                  break ;
	               case 2 :
	                  result = -16.0 * x*(1.0-x) * y*(0.5-y) * 2.0;
	                  break ;	    
	            }
	            break ;
         }
         break ;
      case 8 : // <--------
         switch( a )
         {
            case 0 :
               switch( b )
               {
                  case 0 :
                     result = 8.0 * -2.0 * y*(0.5-y) * (1.0-z)*(0.5-z);
                     break ;
                  case 1 :
                     result = 8.0 * (-2.0*x+0.5)  * (-2.0*y+0.5) * (1.0-z)*(0.5-z);
                     break ;
                  case 2 :
                     result = 8.0 * (-2.0*x+0.5)  * y*(0.5-y) * (2.0*z-1.5);
                     break ; 
               }
               break ;
            case 1 :
               switch( b )
               {
                  case 0 :
                     result = 8.0 * (-2.0*x+0.5) * (-2.0*y+0.5) * (1.0-z)*(0.5-z);
                     break ;
                  case 1 :
                     result = 8.0 * x*(0.5-x) *  -2.0 * (1.0-z)*(0.5-z);
                     break ;
                  case 2 :
                     result = 8.0 * x*(0.5-x) * (-2.0*y+0.5) * (2.0*z-1.5);
                     break ;
               }
               break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = 8.0 * (-2.0*x+0.5) * y*(0.5-y) * (2.0*z-1.5);
		     break ;
		  case 1 :
		     result = 8.0 * x*(0.5-x) * (-2.0*y+0.5) * (2.0*z-1.5);
		     break ;
		  case 2 :
		     result = 8.0 * x*(0.5-x) * y*(0.5-y) *  2.0;
		     break ;	    
               }
	       break ;
	 }
	 break ;

      case 9 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
                  case 0 :
                     result = 16.0 * 2.0 * (1.0-y)*(0.5-y) *  z*(1.0-z) ;
                     break ;
                  case 1 :
                     result = 16.0 * (2.0*x-1.5) * (2.0*y-1.5) *   z*(1.0-z);
                     break ;
                  case 2 :
                     result = 16.0 * (2.0*x-1.5) * (1.0-y)*(0.5-y) * (-2.0*z+1.0) ;
                     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
                     result = 16.0 * (2.0*x-1.5)  * (2.0*y-1.5)  *  z*(1.0-z)  ;
                     break ;
                  case 1 :
                     result = 16.0 * (1.0-x)*(0.5-x) * 2.0 *  z*(1.0-z) ;
                     break ;
                  case 2 :
                     result = 16.0 * (1.0-x)*(0.5-x) * (2.0*y-1.5)  * (-2.0*z+1.0) ;
                     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
                  case 0 :
                     result =  16.0 * (2.0*x-1.5) * (1.0-y)*(0.5-y) *  (-2.0*z+1.0) ;
                     break ;
                  case 1 :
                     result = 16.0 * (1.0-x)*(0.5-x) * (2.0*y-1.5) * (-2.0*z+1.0) ;
                     break ;
                  case 2 :
                     result = 16.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * -2.0 ;
                     break ;
               }
	       break ;
	 }
	 break ;

      case 10 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
                  case 0 :
                     result = 32.0 * -2.0 * (1.0-y)*(0.5-y) * z*(1.0-z)  ;
                     break ;
                  case 1 :
                     result = 32.0 * (-2.0*x+1.0) * (2.0*y-1.5) * z*(1.0-z) ;
                     break ;
                  case 2 :
                     result = 32.0 * (-2.0*x+1.0) * (1.0-y)*(0.5-y) * (-2.0*z+1.0) ;
                     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
                     result =  32.0 * (-2.0*x+1.0) * (2.0*y-1.5) * z*(1.0-z) ;
                     break ;
                  case 1 :
                     result = 32.0 * x*(1.0-x) * 2.0 * z*(1.0-z)  ;
                     break ;
                  case 2 :
                     result = 32.0 * x*(1.0-x) * (2.0*y-1.5) * (-2.0*z+1.0) ;
                     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
                  case 0 :
                     result = 32.0 * (-2.0*x+1.0)  * (1.0-y)*(0.5-y) * (-2.0*z+1.0);
                     break ;
                  case 1 :
                     result = 32.0 * x*(1.0-x) * (2.0*y-1.5)  * (-2.0*z+1.0);
                     break ;
                  case 2 :
                     result = 32.0 * x*(1.0-x) * (1.0-y)*(0.5-y) * -2.0 ;
                     break ;
               }
	       break ;
	 }
	 break ;

      case 11 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
                  case 0 :
                     result = -16.0 * -2.0 * (1.0-y)*(0.5-y) * z*(1.0-z)  ;
                     break ;
                  case 1 :
                     result = -16.0 * (-2.0*x+0.5) * (2.0*y-1.5) *  z*(1.0-z) ;
                     break ;
                  case 2 :
                     result = -16.0 * (-2.0*x+0.5) * (1.0-y)*(0.5-y) * (-2.0*z+1.0)  ;
                     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
		     result = -16.0 * (-2.0*x+0.5) * (2.0*y-1.5) * z*(1.0-z);
		     break ;
		  case 1 :
		     result = -16.0 * x*(0.5-x) * 2.0 * z*(1.0-z) ;
		     break ;
		  case 2 :
		     result = -16.0 * x*(0.5-x) * (2.0*y-1.5) * (-2.0*z+1.0) ;
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = -16.0 * (-2.0*x+0.5) * (1.0-y)*(0.5-y) * (-2.0*z+1.0) ;
		     break ;
		  case 1 :
		     result = -16.0 * x*(0.5-x) * (2.0*y-1.5) * (-2.0*z+1.0) ;
		     break ;
		  case 2 :
		     result = -16.0 * x*(0.5-x) * (1.0-y)*(0.5-y) * -2.0 ;
		     break ;		    
               }
	       break ;
	 }
	 break ;

      case 12 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
                  case 0 :
		     result = -32.0 * -2.0 * y*(1.0-y) * z*(1.0-z) ;
		     break ;
		  case 1 :
		     result = -32.0 * (-2.0*x+0.5) * (-2.0*y+1.0) * z*(1.0-z);
	       break ;
		  case 2 :
		     result = -32.0 * (-2.0*x+0.5) * y*(1.0-y) * (-2.0*z+1.0);
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
		     result = -32.0 * (-2.0*x+0.5) * (-2.0*y+1.0) * z*(1.0-z) ;
		     break ;
		  case 1 :
		     result = -32.0 * x*(0.5-x) * -2.0 * z*(1.0-z);
		     break ;
		  case 2 :
		     result = -32.0 * x*(0.5-x) * (-2.0*y+1.0) * (-2.0*z+1.0);
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = -32.0 * (-2.0*x+0.5) * y*(1.0-y) * (-2.0*z+1.0) ;
		     break ;
		  case 1 :
		     result = -32.0 * x*(0.5-x) * (-2.0*y+1.0) * (-2.0*z+1.0) ;
		     break ;
		  case 2 :
		     result = -32.0 * x*(0.5-x) * y*(1.0-y) * -2.0;
		     break ;	    
               }
	       break ;
	 }
	 break ;

      case 13 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
                  case 0 :
		     result = 64.0 * -2.0 * y*(1.0-y) * z*(1.0-z);
		     break ;
		  case 1 :
		     result = 64.0 * (-2.0*x+1.0) * (-2.0*y+1.0) * z*(1.0-z);
		     break ;
		  case 2 :
		     result = 64.0 * (-2.0*x+1.0) * y*(1.0-y) * (-2.0*z+1.0);
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
		     result = 64.0 * (-2.0*x+1.0) * (-2.0*y+1.0) * z*(1.0-z);
		     break ;
		  case 1 :
		     result = 64.0 * x*(1.0-x) * -2.0 * z*(1.0-z) ;
		     break ;
		  case 2 :
		     result = 64.0 * x*(1.0-x) * (-2.0*y+1.0) * (-2.0*z+1.0);
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = 64.0 * (-2.0*x+1.0) * y*(1.0-y) * (-2.0*z+1.0) ;
		     break ;
		  case 1 :
		     result = 64.0 * x*(1.0-x) * (-2.0*y+1.0) * (-2.0*z+1.0) ;
		     break ;
		  case 2 :
		     result = 64.0 * x*(1.0-x) * y*(1.0-y) * -2.0;
		     break ;	    
               }
	       break ;
	 }
	 break ;
      case 14 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = 32.0 * 2.0 * y*(1.0-y) * z*(1.0-z);
		     break ;
		  case 1 :
		     result = 32.0 * (2.0*x-1.5) * (-2.0*y+1.0) * z*(1.0-z);
		     break ;
		  case 2 :
		     result = 32.0 * (2.0*x-1.5) * y*(1.0-y) * (-2.0*z+1.0);
		     break ;
	       }
	       break ;
	    case 1 :
	       switch( b )
               {
		  case 0 :
		     result = 32.0 * (2.0*x-1.5) * (-2.0*y+1.0) * z*(1.0-z);
		     break ;
		  case 1 :
		     result = 32.0 * (1.0-x)*(0.5-x) * -2.0 * z*(1.0-z);
		     break ;
		  case 2 :
		     result = 32.0 * (1.0-x)*(0.5-x) * (-2.0*y+1.0)  * (-2.0*z+1.0);
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = 32.0 * (2.0*x-1.5) * y*(1.0-y) * (-2.0*z+1.0);
		     break ;
		  case 1 :
		     result = 32.0 * (1.0-x)*(0.5-x) * (-2.0*y+1.0) * (-2.0*z+1.0) ;
		     break ;
		  case 2 :
		     result = 32.0 * (1.0-x)*(0.5-x) * y*(1.0-y) * -2.0;
		     break ;		     
               }
	       break ;
	 }
	 break ;
      case 15 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = -16.0 * 2.0 * y*(0.5-y) * z*(1.0-z);
		     break ;
		  case 1 :
		     result = -16.0 * (2.0*x-1.5) * (-2.0*y+0.5) * z*(1.0-z);
		     break ;
		  case 2 :
		     result = -16.0 * (2.0*x-1.5) * y*(0.5-y) * (-2.0*z+1.0);
	       break ;
	       }
	       break ;
	    case 1 :
	       switch( b )
               {
		  case 0 :
		     result = -16.0 * (2.0*x-1.5) * (-2.0*y+0.5) * z*(1.0-z);
		     break ;
		  case 1 :
		     result = -16.0 * (1.0-x)*(0.5-x) * -2.0 * z*(1.0-z);
		     break ;
		  case 2 :
		     result = -16.0 * (1.0-x)*(0.5-x) * (-2.0*y+0.5)  * (-2.0*z+1.0);
		     break ;   
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = -16.0 * (2.0*x-1.5) * y*(0.5-y) * (-2.0*z+1.0);
		     break ;
		  case 1 :
		     result = -16.0 * (1.0-x)*(0.5-x) * (-2.0*y+0.5) * (-2.0*z+1.0) ;
		     break ;
		  case 2 :
		     result = -16.0 * (1.0-x)*(0.5-x) * y*(0.5-y) * -2.0;
		     break ; 	    
               }
	       break ;
	 }
	 break ;    
      case 16 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
                  case 0 :
		     result = -32.0 * -2.0 * y*(0.5-y) * z*(1.0-z);
		     break ;
		  case 1 :
		     result = -32.0 * (-2.0*x+1.0) * (-2.0*y+0.5) * z*(1.0-z);
		     break ;
		  case 2 :
		     result = -32.0 * (-2.0*x+1.0) * y*(0.5-y) * (-2.0*z+1.0);
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
		     result = -32.0 * (-2.0*x+1.0) * (-2.0*y+0.5) * z*(1.0-z);
		     break ;
		  case 1 :
		     result = -32.0 * x*(1.0-x) * -2.0 * z*(1.0-z);
		     break ;
		  case 2 :
		     result = -32.0 * x*(1.0-x) * (-2.0*y+0.5) * (-2.0*z+1.0);
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		 case 0 :
		    result = -32.0 * (-2.0*x+1.0) * y*(0.5-y) * (-2.0*z+1.0) ;
		    break ;
		  case 1 :
		     result = -32.0 * x*(1.0-x) * (-2.0*y+0.5) * (-2.0*z+1.0) ;
		     break ;
		  case 2 :
		     result = -32.0 * x*(1.0-x) * y*(0.5-y) * -2.0;
		     break ;	    
               }
	       break ;
	 }
	 break ;
      case 17 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = 16.0 * -2.0 * y*(0.5-y) * z*(1.0-z) ;
		     break ;
		  case 1 :
		     result = 16.0 * (-2.0*x+0.5)  * (-2.0*y+0.5) * z*(1.0-z);
		     break ;
		  case 2 :
		     result = 16.0 * (-2.0*x+0.5)  * y*(0.5-y) * (-2.0*z+1.0);
		     break ; 
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
		     result = 16.0 * (-2.0*x+0.5) * (-2.0*y+0.5) * z*(1.0-z) ;
		     break ;
		  case 1 :
		     result = 16.0 * x*(0.5-x) *  -2.0 * z*(1.0-z);
		     break ;
		  case 2 :
		     result = 16.0 * x*(0.5-x) * (-2.0*y+0.5) * (-2.0*z+1.0);
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = 16.0 * (-2.0*x+0.5) * y*(0.5-y) * (-2.0*z+1.0);
		     break ;
		  case 1 :
		     result = 16.0 * x*(0.5-x) * (-2.0*y+0.5) * (-2.0*z+1.0);
		     break ;
		  case 2 :
		     result = 16.0 * x*(0.5-x) * y*(0.5-y) *  -2.0;
		     break ;	    
               }
	       break ;
	 }
	 break ;
      case 18 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = -8.0 * 2.0 * (1.0-y)*(0.5-y) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = -8.0 * (2.0*x-1.5) * (2.0*y-1.5)  * z*(0.5-z);
		     break ;
		  case 2 :
		     result = -8.0 * (2.0*x-1.5) * (1.0-y)*(0.5-y) * (-2.0*z+0.5) ;
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
		  case 0 :
		     result = -8.0 * (2.0*x-1.5) * (2.0*y-1.5) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = -8.0 * (1.0-x)*(0.5-x) * 2.0 * z*(0.5-z);
		     break ;
		  case 2 :
		     result = -8.0 * (1.0-x)*(0.5-x) * (2.0*y-1.5) * (-2.0*z+0.5) ;
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		   case 0 :
		      result = -8.0 * (2.0*x-1.5) * (1.0-y)*(0.5-y) * (-2.0*z+0.5) ;
		      break ;
		  case 1 :
		     result = -8.0 * (1.0-x)*(0.5-x) * (2.0*y-1.5)  * (-2.0*z+0.5);
		     break ;
		  case 2 :
		     result = -8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * -2.0 ;
		     break ;
               }
	       break ;
	 }
	 break ;
      case 19 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = -16.0 * -2.0 * (1.0-y)*(0.5-y) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = -16.0 * (-2.0*x+1.0) * (2.0*y-1.5) * z*(0.5-z);
		     break ;
		  case 2 :
		     result = -16.0 * (-2.0*x+1.0) * (1.0-y)*(0.5-y) * (-2.0*z+0.5) ;
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
		  case 0 :
		     result = -16.0 * (-2.0*x+1.0) * (2.0*y-1.5) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = -16.0 * x*(1.0-x) * 2.0 * z*(0.5-z);
		     break ;
		  case 2 :
		     result = -16.0 * x*(1.0-x) * (2.0*y-1.5) * (-2.0*z+0.5) ;
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = -16.0 * (-2.0*x+1.0) * (1.0-y)*(0.5-y) *(-2.0*z+0.5);
		     break ;
		  case 1 :
		     result = -16.0 * x*(1.0-x) * (2.0*y-1.5) * (-2.0*z+0.5);
		     break ;
		  case 2 :
		     result = -16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) * -2.0 ;
		     break ;
               }
	       break ;
	 }
	 break ;
      case 20 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = 8.0 * -2.0 * (1.0-y)*(0.5-y) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = 8.0 * (-2.0*x+0.5) * (2.0*y-1.5) * z*(0.5-z) ;
		     break ;
		  case 2 :
		     result = 8.0 * (-2.0*x+0.5) * (1.0-y)*(0.5-y) *  (-2.0*z+0.5);
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
		  case 0 :
		     result = 8.0 * (-2.0*x+0.5) *  (2.0*y-1.5) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = 8.0 * x*(0.5-x) * 2.0 * z*(0.5-z) ;
		     break ;
		  case 2 :
		     result = 8.0 * x*(0.5-x) *  (2.0*y-1.5) *  (-2.0*z+0.5);
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = 8.0 * (-2.0*x+0.5) * (1.0-y)*(0.5-y) * (-2.0*z+0.5);
		     break ;
		  case 1 :
		     result = 8.0 * x*(0.5-x) * (2.0*y-1.5) * (-2.0*z+0.5) ;
		     break ;
		  case 2 :
		     result = 8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) *  -2.0;
		     break ;
               }
	       break ;
	 }
	 break ;
      case 21 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = 16.0 * -2.0 * y*(1.0-y) * z*(0.5-z) ;
		     break ;
		  case 1 :
		     result = 16.0 *(-2.0*x+0.5) * (-2.0*y+1.0) * z*(0.5-z) ;
		     break ;
		  case 2 :
		     result = 16.0 * (-2.0*x+0.5)  * y*(1.0-y) * (-2.0*z+0.5);
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
		     result = 16.0 * (-2.0*x+0.5) * (-2.0*y+1.0) * z*(0.5-z) ;
		     break ;
		  case 1 :
		     result = 16.0 * x*(0.5-x) * -2.0 * z*(0.5-z) ;
		     break ;
		  case 2 :
		     result = 16.0 * x*(0.5-x) * (-2.0*y+1.0) * (-2.0*z+0.5);
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = 16.0 * (-2.0*x+0.5) * y*(1.0-y) * (-2.0*z+0.5) ;
		     break ;
		  case 1 :
		     result = 16.0 * x*(0.5-x) * (-2.0*y+1.0) * (-2.0*z+0.5) ;
		     break ;
		  case 2 :
		     result = 16.0 * x*(0.5-x) * y*(1.0-y) * -2.0;
		     break ;
               }
	       break ;
	 }
	 break ;
      case 22 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = -32.0 * -2.0 * y*(1.0-y) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = -32.0 * (-2.0*x+1.0) * (-2.0*y+1.0) * z*(0.5-z) ;
		     break ;
		  case 2 :
		     result = -32.0 * (-2.0*x+1.0) * y*(1.0-y) * (-2.0*z+0.5);
		     break ;	 
               }
	       break ;
	    case 1 :
	       switch( b )
               {
		  case 0 :
		     result = -32.0 * (-2.0*x+1.0) * (-2.0*y+1.0) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = -32.0 * x*(1.0-x) * -2.0 * z*(0.5-z) ;
		     break ;
		  case 2 :
		     result = -32.0 * x*(1.0-x) * (-2.0*y+1.0) * (-2.0*z+0.5);
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = -32.0 * (-2.0*x+1.0) * y*(1.0-y) * (-2.0*z+0.5);
		     break ;
		  case 1 :
		     result = -32.0 * x*(1.0-x) * (-2.0*y+1.0) * (-2.0*z+0.5) ;
		     break ;
		  case 2 :
		     result = -32.0 * x*(1.0-x) * y*(1.0-y) * -2.0;
		     break ;
               }
	       break ;
	 }
	 break ;
      case 23 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = -16.0 * 2.0 * y*(1.0-y) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = -16.0 * (2.0*x-1.5) * (-2.0*y+1.0) * z*(0.5-z);
		     break ;
		  case 2 :
		     result = -16.0 * (2.0*x-1.5) * y*(1.0-y) * (-2.0*z+0.5) ;
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
		     result = -16.0 * (2.0*x-1.5) * (-2.0*y+1.0) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = -16.0 * (1.0-x)*(0.5-x) * -2.0 * z*(0.5-z);
		     break ;
		  case 2 :
		     result = -16.0 * (1.0-x)*(0.5-x) * (-2.0*y+1.0) * (-2.0*z+0.5) ;
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = -16.0 * (2.0*x-1.5) * y*(1.0-y) * (-2.0*z+0.5);
		     break ;
		  case 1 :
		     result = -16.0 * (1.0-x)*(0.5-x) * (-2.0*y+1.0) * (-2.0*z+0.5) ;
		     break ;
		  case 2 :
		     result = -16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) * -2.0 ;
		     break ;
               }
	       break ;
	 }
	 break ;
      case 24 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = 8.0 * 2.0 * y*(0.5-y) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = 8.0 * (2.0*x-1.5) * (-2.0*y+0.5) * z*(0.5-z);
		     break ;
		  case 2 :
		     result = 8.0 * (2.0*x-1.5) * y*(0.5-y) * (-2.0*z+0.5) ;
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
		     result = 8.0 * (2.0*x-1.5) * (-2.0*y+0.5) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = 8.0 * (1.0-x)*(0.5-x) * -2.0 * z*(0.5-z);
		     break ;
		  case 2 :
		     result = 8.0 * (1.0-x)*(0.5-x) *(-2.0*y+0.5)  * (-2.0*z+0.5) ;
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = 8.0 * (2.0*x-1.5) * y*(0.5-y) * (-2.0*z+0.5);
		     break ;
		  case 1 :
		     result = 8.0 * (1.0-x)*(0.5-x) * (-2.0*y+0.5) * (-2.0*z+0.5) ;
		     break ;
		  case 2 :
		     result = 8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) * -2.0 ;
		     break ;
               }
	       break ;
	 }
	 break ;
      case 25 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = 16.0 * -2.0 * y*(0.5-y) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = 16.0 * (-2.0*x+1.0) * (-2.0*y+0.5) * z*(0.5-z);
		     break ;
		  case 2 :
		     result = 16.0 * (-2.0*x+1.0) * y*(0.5-y) * (-2.0*z+0.5)  ;
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
                  case 0 :
		     result = 16.0 * (-2.0*x+1.0) * (-2.0*y+0.5) * z*(0.5-z);
		     break ;
		  case 1 :
		     result = 16.0 * x*(1.0-x) * -2.0 * z*(0.5-z);
		     break ;
		  case 2 :
		     result = 16.0 * x*(1.0-x) * (-2.0*y+0.5) * (-2.0*z+0.5)  ;
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = 16.0 * (-2.0*x+1.0) * y*(0.5-y) * (-2.0*z+0.5);
		     break ;
		  case 1 :
		     result = 16.0 * x*(1.0-x) * (-2.0*y+0.5) * (-2.0*z+0.5);
		     break ;
		  case 2 :
		     result = 16.0 * x*(1.0-x) * y*(0.5-y) * -2.0  ;
		     break ;
               }
	       break ;
	 }
	 break ;
      case 26 : // <--------
	 switch( a )
         {
	    case 0 :
	       switch( b )
               {
		  case 0 :
		     result = -8.0 * -2.0 * y*(0.5-y) * z*(0.5-z) ;
		     break ;
		  case 1 :
		     result = -8.0 * (-2.0*x+0.5) * (-2.0*y+0.5) * z*(0.5-z);
		     break ;
		  case 2 :
		     result = -8.0 * (-2.0*x+0.5) * y*(0.5-y) * (-2.0*z+0.5);
		     break ;
               }
	       break ;
	    case 1 :
	       switch( b )
               {
		  case 0 :
		     result = -8.0 * (-2.0*x+0.5) * (-2.0*y+0.5) * z*(0.5-z) ;
		     break ;
		  case 1 :
		     result = -8.0 * x*(0.5-x) * -2.0 * z*(0.5-z);
		     break ;
		  case 2 :
		     result = -8.0 * x*(0.5-x) * (-2.0*y+0.5) * (-2.0*z+0.5);
		     break ;
               }
	       break ;
	    case 2 :
	       switch( b )
               {
		  case 0 :
		     result = -8.0 * (-2.0*x+0.5) * y*(0.5-y) * (-2.0*z+0.5) ;
		     break ;
		  case 1 :
		     result = -8.0 * x*(0.5-x) * (-2.0*y+0.5) * (-2.0*z+0.5);
		     break ;
		  case 2 :
		     result = -8.0 * x*(0.5-x) * y*(0.5-y) * -2.0;
		     break ;
               }
	       break ;
	 }
	 break ;
    
   }
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
PDE_3D_Q2_27nodes:: reference_polyhedron_POST(
                               GE_ReferencePolyhedron const* result ) const 
//----------------------------------------------------------------------------
{
   PEL_ASSERT( PDE_ReferenceElement::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceCube::object() ) ;
   return( true ) ;
}

