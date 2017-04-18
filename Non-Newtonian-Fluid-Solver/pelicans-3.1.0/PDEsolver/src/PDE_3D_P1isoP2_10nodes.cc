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

#include <PDE_3D_P1isoP2_10nodes.hh>

#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceTetrahedron.hh>

PDE_3D_P1isoP2_10nodes const* 
PDE_3D_P1isoP2_10nodes:: REGISTRATOR = new PDE_3D_P1isoP2_10nodes() ;

//----------------------------------------------------------------------
PDE_3D_P1isoP2_10nodes:: ~PDE_3D_P1isoP2_10nodes( void )
//----------------------------------------------------------------------
{
   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
PDE_3D_P1isoP2_10nodes:: PDE_3D_P1isoP2_10nodes( void )
//----------------------------------------------------------------------
   : PDE_ReferenceElement( "PDE_3D_P1isoP2_10nodes",
                           GE_ReferenceTetrahedron::object() )
{
   append_node( GE_Point::create( this, 0.0, 0.0, 0.0 ) ) ; // 0
   append_node( GE_Point::create( this, 1.0, 0.0, 0.0 ) ) ; // 1
   append_node( GE_Point::create( this, 0.0, 1.0, 0.0 ) ) ; // 2
   append_node( GE_Point::create( this, 0.0, 0.0, 1.0 ) ) ; // 3
   append_node( GE_Point::create( this, 0.5, 0.0, 0.0 ) ) ; // 4
   append_node( GE_Point::create( this, 0.5, 0.5, 0.0 ) ) ; // 5
   append_node( GE_Point::create( this, 0.0, 0.5, 0.0 ) ) ; // 6
   append_node( GE_Point::create( this, 0.0, 0.0, 0.5 ) ) ; // 7
   append_node( GE_Point::create( this, 0.5, 0.0, 0.5 ) ) ; // 8
   append_node( GE_Point::create( this, 0.0, 0.5, 0.5 ) ) ; // 9
}

//----------------------------------------------------------------------
double
PDE_3D_P1isoP2_10nodes:: N_local( size_t node,
                                  GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_P1isoP2_10nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;
   
   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;
   
   double result = 0.0 ;
   switch( node )
   {
      case 0 :
         if( t_0467( x, y, z ) ) result = 1.0 - 2.0*( x + y + z ) ;
         break ;
      case 1 :
         if( t_4158( x, y, z ) ) result = 2.0*( x - 0.5 ) ;
         break ;
      case 2 :
         if( t_6529( x, y, z ) ) result = 2.0*( y - 0.5 ) ;
         break ;
      case 3 :
         if( t_7893( x, y, z ) ) result = 2.0*( z - 0.5 ) ;
         break ;
      case 4 :
         if( t_0467( x, y, z ) || t_9647( x, y, z ) )
         {
            result = 2.0*x ;
         }
         else if( t_4158(x,y,z) || t_8549(x,y,z) )
         {
            result = 2.0*( 1.0 - x - y - z ) ;
         }
         else if( t_8794( x, y, z ) )
         {
            result = 1.0 - 2.0*z ;
         }
         else if( t_6459( x, y, z ) )
         {
            result = 1.0 - 2.0*y  ;
         }
         break ;
      case 5 :
         if( t_4158( x, y, z ) )
         {
            result = 2.0*y ;
         }
         else if( t_6529( x, y, z ) ) 
         {
            result = 2.0*x ;
         }
         else if( t_8549( x, y, z ) || t_6459( x, y, z ) ) 
         {
            result = 2.0*( x + y - 0.5 ) ;
         }
         break ;
      case 6 :
         if( t_0467( x, y, z ) )
         {
            result = 2.0*y ;
         }
         else if( t_6529( x, y, z ) )
         {
            result = 2.0*( 1.0 - x - y - z ) ;
         }
         else if( t_6459( x, y, z ) || t_9647( x, y, z ) ) 
         {
            result = 2.0*( 0.5 - x - z ) ;
         }
         break ;
      case 7 :
         if( t_0467( x, y, z ) )
         {
            result = 2.0*z ;
         }
         else if( t_7893( x, y, z ) )
         {
            result = 2.0*( 1.0 - x - y - z ) ;
         }
         else if( t_8794( x, y, z ) || t_9647( x, y, z ) )
         {
            result = 2.0*( 0.5 - x - y ) ;
         }
         break ;
      case 8 :
         if( t_4158( x, y, z ) )
         {
            result = 2.0*z ;
         }
         else if( t_7893( x, y, z ) )
         {
            result = 2.0*x ;
         }
         else if( t_8794( x, y, z ) || t_8549( x, y, z ) )
         {
            result = 2.0*( x + z - 0.5 ) ;
         }
         break ;
      case 9 :
         if( t_6529( x, y, z ) || t_6459( x, y, z ) )
         {
            result = 2.0*z ;
         }
         else if( t_7893( x, y, z ) || t_8794( x, y, z ) )
         {
            result = 2.0*y ;
         }
         else if( t_8549( x, y, z ) ) 
         {
            result = 1.0 - 2.0*x ;
         }
         else if( t_9647( x, y, z ) )
         {
            result = 2.0*( x + y + z - 0.5 ) ;
         }
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_3D_P1isoP2_10nodes:: dN_local( size_t node, size_t a,
                                   GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_P1isoP2_10nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;
   
   double result = 0.0 ;
   switch( node )
   {
      case 0 :
         if( t_0467( x, y, z ) ) result = -2.0 ;
         break ;
      case 1 :
         if( ( a == 0 ) && t_4158( x, y, z ) ) result = 2.0 ;
         break ;
      case 2 :
         if( ( a == 1 ) && t_6529( x, y, z ) ) result = 2.0 ;
         break ;
      case 3 :
         if( ( a == 2 ) && t_7893( x, y, z ) ) result = 2.0 ;
         break ;
      case 4 :
         if( ( a == 0 ) && ( t_0467( x, y, z ) || t_9647( x, y, z ) ) )
         {
            result = 2.0 ;
         }
         else if( t_4158(x,y,z) || t_8549(x,y,z) )
         {
            result = -2.0 ;
         }
         else if( ( a == 2 ) && t_8794( x, y, z ) )
         {
            result = -2.0 ;
         }
         else if( ( a == 1 ) && t_6459( x, y, z ) )
         {
            result = -2.0  ;
         }
         break ;
      case 5 :
         if( ( a == 1 ) && t_4158( x, y, z ) )
         {
            result = 2.0 ;
         }
         else if( ( a == 0 ) && t_6529( x, y, z ) ) 
         {
            result = 2.0 ;
         }
         else if( ( a != 2 ) && ( t_8549( x, y, z ) || t_6459( x, y, z ) ) ) 
         {
            result = 2.0 ;
         }
         break ;
      case 6 :
         if( ( a == 1 ) && t_0467( x, y, z ) )
         {
            result = 2.0;
         }
         else if( t_6529( x, y, z ) )
         {
            result = -2.0 ;
         }
         else if( ( a != 1 ) && ( t_6459( x, y, z ) || t_9647( x, y, z ) ) ) 
         {
            result = -2.0 ;
         }
         break ;
      case 7 :
         if( ( a == 2 ) && t_0467( x, y, z ) )
         {
            result = 2.0 ;
         }
         else if( t_7893( x, y, z ) )
         {
            result = -2.0 ;
         }
         else if( ( a != 2 ) && ( t_8794( x, y, z ) || t_9647( x, y, z ) ) )
         {
            result = -2.0 ;
         }
         break ;
      case 8 :
         if( ( a == 2 ) && t_4158( x, y, z ) )
         {
            result = 2.0 ;
         }
         else if( ( a == 0 ) && t_7893( x, y, z ) )
         {
            result = 2.0 ;
         }
         else if( ( a != 1 ) && ( t_8794( x, y, z ) || t_8549( x, y, z ) ) )
         {
            result = 2.0 ;
         }
         break ;
      case 9 :
         if( ( a == 2 ) && ( t_6529( x, y, z ) || t_6459( x, y, z ) ) )
         {
            result = 2.0 ;
         }
         else if( ( a == 1 ) && ( t_7893( x, y, z ) || t_8794( x, y, z ) ) )
         {
            result = 2.0 ;
         }
         else if( ( a == 0 ) && t_8549( x, y, z ) ) 
         {
            result = -2.0 ;
         }
         else if( t_9647( x, y, z ) )
         {
            result = 2.0 ;
         }
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_3D_P1isoP2_10nodes:: d2N_local( size_t node,
                                    size_t a, size_t b,
                                    GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_P1isoP2_10nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   return( 0.0 ) ;
}

//----------------------------------------------------------------------
bool
PDE_3D_P1isoP2_10nodes:: t_0467( double x, double y, double z ) const
//----------------------------------------------------------------------
{
   return( x + y + z <= 0.5 ) ;
}

//----------------------------------------------------------------------
bool
PDE_3D_P1isoP2_10nodes:: t_4158( double x, double y, double z ) const
//----------------------------------------------------------------------
{
   return( x >= 0.5 ) ;
}

//----------------------------------------------------------------------
bool
PDE_3D_P1isoP2_10nodes:: t_6529( double x, double y, double z ) const
//----------------------------------------------------------------------
{
   return( y >= 0.5 ) ;
}

//----------------------------------------------------------------------
bool
PDE_3D_P1isoP2_10nodes:: t_7893( double x, double y, double z ) const
//----------------------------------------------------------------------
{
   return( z >= 0.5 ) ;
}

//----------------------------------------------------------------------
bool
PDE_3D_P1isoP2_10nodes:: t_8794( double x, double y, double z ) const
//----------------------------------------------------------------------
{
   return( ( z <= 0.5 ) && ( x + z >= 0.5 ) && ( x + y <= 0.5 ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_3D_P1isoP2_10nodes:: t_8549( double x, double y, double z ) const
//----------------------------------------------------------------------
{
   return( ( x <= 0.5 ) && ( x + y >= 0.5 ) && ( x + z >= 0.5 ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_3D_P1isoP2_10nodes:: t_6459( double x, double y, double z ) const
//----------------------------------------------------------------------
{
   return( ( y <= 0.5 ) && ( x + z <= 0.5 ) && ( x + y >= 0.5 ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_3D_P1isoP2_10nodes:: t_9647( double x, double y, double z ) const
//----------------------------------------------------------------------
{
   return( ( x + y <= 0.5 ) && ( x + z <= 0.5 ) && ( x + y + z >= 0.5 ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_3D_P1isoP2_10nodes:: reference_polyhedron_POST(
                               GE_ReferencePolyhedron const* result ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( PDE_ReferenceElement::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceTetrahedron::object() ) ;
   return( true ) ;
}

