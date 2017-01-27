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

#include <PDE_3D_Q1isoQ2_27nodes.hh>

#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceCube.hh>

PDE_3D_Q1isoQ2_27nodes const* 
PDE_3D_Q1isoQ2_27nodes:: REGISTRATOR = new PDE_3D_Q1isoQ2_27nodes() ;

//----------------------------------------------------------------------
PDE_3D_Q1isoQ2_27nodes:: ~PDE_3D_Q1isoQ2_27nodes( void )
//----------------------------------------------------------------------
{
   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
PDE_3D_Q1isoQ2_27nodes:: PDE_3D_Q1isoQ2_27nodes ( void )
//----------------------------------------------------------------------
   : PDE_ReferenceElement( "PDE_3D_Q1isoQ2_27nodes",
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
PDE_3D_Q1isoQ2_27nodes:: N_local( size_t node,
                                  GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_Q1isoQ2_27nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;
   
   double result = 0.0 ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   switch( node )
   {
      case 0 :
         if( x<0.5 && y<0.5 && z<0.5 ) 
            result = (1.0-2.0*x) * (1.0-2.0*y) * (1.0-2.0*z) ;
         break ;
      case 1 :
         if( y<0.5 && z<0.5 )
         {
            if( x<0.5 ) result = 2.0*x       * (1.0-2.0*y) * (1.0-2.0*z) ;
            else        result = 2.0*(1.0-x) * (1.0-2.0*y) * (1.0-2.0*z) ;
         }
         break ;
      case 2 :
         if( x>0.5 && y<0.5 && z<0.5 ) 
            result = (2.0*x-1.0) * (1.0-2.0*y) * (1.0-2.0*z) ;
         break ;
      case 3 :
         if( x>0.5 && z<0.5 )
         {
            if( y<0.5 ) result = 2.0*y       * (2.0*x-1.0) * (1.0-2.0*z) ;
            else        result = 2.0*(1.0-y) * (2.0*x-1.0) * (1.0-2.0*z) ;
         }
         break ;
      case 4 :
         if( z<0.5 )
         {
            if( x<0.5 )
            {
               if( y<0.5 ) result = 4.0*x*y       * (1.0-2.0*z) ;
               else        result = 4.0*x*(1.0-y) * (1.0-2.0*z) ;
            }
            else
            {
               if( y<0.5 ) result = 4.0*y*(1.0-x)       * (1.0-2.0*z) ;
               else        result = 4.0*(1.0-x)*(1.0-y) * (1.0-2.0*z) ;
            }
         }
         break ;
      case 5 :
         if( x<0.5 && z<0.5 )
         {
            if( y<0.5 ) result = 2.0*y       * (1.0-2.0*x) * (1.0-2.0*z) ;
            else        result = 2.0*(1.0-y) * (1.0-2.0*x) * (1.0-2.0*z) ;
         }
         break ;
      case 6 :
         if( x<0.5 && y>0.5 && z<0.5 ) 
            result = (1.0-2.0*x) * (2.0*y-1.0) * (1.0-2.0*z) ;
         break ;
      case 7 :
         if( y>0.5 && z<0.5 )
         {
            if( x<0.5 ) result = 2.0*x       * (2.0*y-1.0) * (1.0-2.0*z) ;
            else        result = 2.0*(1.0-x) * (2.0*y-1.0) * (1.0-2.0*z) ;
            
         }
         break ;
      case 8 :
         if( x>0.5 && y>0.5 && z<0.5 ) 
            result = (2.0*x-1.0) * (2.0*y-1.0) * (1.0-2.0*z) ;
         break ;
      case 9 :
         if( y<0.5 && x<0.5 )
         {
            if( z<0.5 ) result = 2.0*z       * (1.0-2.0*y) * (1.0-2.0*x) ;
            else        result = 2.0*(1.0-z) * (1.0-2.0*y) * (1.0-2.0*x) ;
         }
         break ;
      case 10 :
         if( y<0.5 )
         {
            if( z<0.5 )
            {
               if( x<0.5 ) result = 4.0*z*x       * (1.0-2.0*y) ;
               else        result = 4.0*z*(1.0-x) * (1.0-2.0*y) ;
            }
            else
            {
               if( x<0.5 ) result = 4.0*x*(1.0-z)       * (1.0-2.0*y) ;
               else        result = 4.0*(1.0-z)*(1.0-x) * (1.0-2.0*y) ;
            }
         }
         break ;
      case 11 :
         if( y<0.5 && x>0.5 )
         {
            if( z<0.5 ) result = 2.0*z       * (1.0-2.0*y) * (2.0*x-1.0) ;
            else        result = 2.0*(1.0-z) * (1.0-2.0*y) * (2.0*x-1.0) ;
         }
         break ;
      case 12 :
         if( x>0.5 )
         {
            if( z<0.5 )
            {
               if( y<0.5 ) result = 4.0*z*y       * (2.0*x-1.0) ;
               else        result = 4.0*z*(1.0-y) * (2.0*x-1.0) ;
            }
            else
            {
               if( y<0.5 ) result = 4.0*y*(1.0-z)       * (2.0*x-1.0) ;
               else        result = 4.0*(1.0-z)*(1.0-y) * (2.0*x-1.0) ;
            }
         }
         break ;
      case 13 :
         if( z<0.5 )
         {
            if( x<0.5 )
            {
               if( y<0.5 ) result = 8.0*x*y*z ;
               else        result = 8.0*x*(1.0-y)*z ;
            }
            else
            {
               if( y<0.5 ) result = 8.0*y*(1.0-x)*z ;
               else        result = 8.0*(1.0-x)*(1.0-y)*z ;
            }
         }
         else
         {
            if( x<0.5 )
            {
               if( y<0.5 ) result = 8.0*x*y*(1.0-z) ;
               else        result = 8.0*x*(1.0-y)*(1.0-z) ;
            }
            else
            {
               if( y<0.5 ) result = 8.0*y*(1.0-x)*(1.0-z) ;
               else        result = 8.0*(1.0-x)*(1.0-y)*(1.0-z) ;
            }
            
         }
         break ;
      case 14 :
         if( x<0.5 )
         {
            if( z<0.5 )
            {
               if( y<0.5 ) result = 4.0*z*y       * (1.0-2.0*x) ;
               else        result = 4.0*z*(1.0-y) * (1.0-2.0*x) ;
            }
            else
            {
               if( y<0.5 ) result = 4.0*y*(1.0-z)       * (1.0-2.0*x) ;
               else        result = 4.0*(1.0-z)*(1.0-y) * (1.0-2.0*x) ;
            }
         }
         break ;
      case 15 :
         if( y>0.5 && x<0.5 )
         {
            if( z<0.5 ) result = 2.0*z       * (2.0*y-1.0) * (1.0-2.0*x) ;
            else        result = 2.0*(1.0-z) * (2.0*y-1.0) * (1.0-2.0*x) ;
         }
         break ;
      case 16 :
         if( y>0.5 )
         {
            if( z<0.5 )
            {
               if( x<0.5 ) result = 4.0*z*x       * (2.0*y-1.0) ;
               else        result = 4.0*z*(1.0-x) * (2.0*y-1.0) ;
            }
            else
            {
               if( x<0.5 ) result = 4.0*x*(1.0-z)       * (2.0*y-1.0) ;
               else        result = 4.0*(1.0-z)*(1.0-x) * (2.0*y-1.0) ;
            }
         }
         break ;
      case 17 :
         if( y>0.5 && x>0.5 )
         {
            if( z<0.5 ) result = 2.0*z       * (2.0*y-1.0) * (2.0*x-1.0) ;
            else        result = 2.0*(1.0-z) * (2.0*y-1.0) * (2.0*x-1.0) ;
         }
         break ;
      case 18 :
         if( x<0.5 && y<0.5 && z>0.5 ) 
            result = (1.0-2.0*x) * (1.0-2.0*y) * (2.0*z-1.0) ;
         break ;
      case 19 :
         if( y<0.5 && z>0.5 )
         {
            if( x<0.5 ) result = 2.0*x       * (1.0-2.0*y) * (2.0*z-1.0) ;
            else        result = 2.0*(1.0-x) * (1.0-2.0*y) * (2.0*z-1.0) ;
         }
         break ;
      case 20 :
         if( x>0.5 && y<0.5 && z>0.5 ) 
            result = (2.0*x-1.0) * (1.0-2.0*y) * (2.0*z-1.0) ;
         break ;
      case 21 :
         if( x>0.5 && z>0.5 )
         {
            if( y<0.5 ) result = 2.0*y       * (2.0*x-1.0) * (2.0*z-1.0) ;
            else        result = 2.0*(1.0-y) * (2.0*x-1.0) * (2.0*z-1.0) ;
         }
         break ;
      case 22 :
         if( z>0.5 )
         {
            if( x<0.5 )
            {
               if( y<0.5 ) result = 4.0*x*y       * (2.0*z-1.0) ;
               else        result = 4.0*x*(1.0-y) * (2.0*z-1.0) ;
            }
            else
            {
               if( y<0.5 ) result = 4.0*y*(1.0-x)       * (2.0*z-1.0) ;
               else        result = 4.0*(1.0-x)*(1.0-y) * (2.0*z-1.0) ;
            }
         }
         break ;
      case 23 :
         if( x<0.5 && z>0.5 )
         {
            if( y<0.5 ) result = 2.0*y       * (1.0-2.0*x) * (2.0*z-1.0) ;
            else        result = 2.0*(1.0-y) * (1.0-2.0*x) * (2.0*z-1.0) ;
         }
         break ;
      case 24 :
         if( x<0.5 && y>0.5 && z>0.5 ) 
            result = (1.0-2.0*x) * (2.0*y-1.0) * (2.0*z-1.0) ;
         break ;
      case 25 :
         if( y>0.5 && z>0.5 )
         {
            if( x<0.5 ) result = 2.0*x       * (2.0*y-1.0) * (2.0*z-1.0) ;
            else        result = 2.0*(1.0-x) * (2.0*y-1.0) * (2.0*z-1.0) ;
            
         }
         break ;
      case 26 :
         if( x>0.5 && y>0.5 && z>0.5 ) 
            result = (2.0*x-1.0) * (2.0*y-1.0) * (2.0*z-1.0) ;
         break ; 
   }
   return( result ) ;
}



//----------------------------------------------------------------------
double
PDE_3D_Q1isoQ2_27nodes:: dN_local( size_t node, size_t a,
                                   GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_Q1isoQ2_27nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double result = 0.0 ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   switch( node )
   {
      case 0 :
         if( x<0.5 && y<0.5 && z<0.5 )
            // bf=(1.0-2.0*x)*(1.0-2.0*y)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 :
                  result = -2.0*(1.0-2.0*y)*(1.0-2.0*z) ;
                  break ;
               case 1 :
                  result = -2.0*(1.0-2.0*x)*(1.0-2.0*z) ;
                  break ;
               case 2 :
                  result = -2.0*(1.0-2.0*x)*(1.0-2.0*y) ;
                  break ;
            }
         break ;
      case 1 :
         if( y<0.5 && z<0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x*(1.0-2.0*y)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 :
                     result = 2.0*(1.0-2.0*y)*(1.0-2.0*z) ;
                     break ;
                  case 1 :
                     result = -4.0*x*(1.0-2.0*z) ;
                     break ;
                  case 2 :
                     result = -4.0*x*(1.0-2.0*y) ;
                     break ;
               }
            else // bf = 2.0*(1.0-x)*(1.0-2.0*y)*(1.0-2.0*z)
               switch( a )
               {
                  case 0 :
                     result = -2.0*(1.0-2.0*y)*(1.0-2.0*z) ;
                     break ;
                  case 1 :
                     result = -4.0*(1.0-x)*(1.0-2.0*z) ;
                     break ;
                  case 2 :
                     result = -4.0*(1.0-x)*(1.0-2.0*y) ;
                     break ;
               }
         }
         break ;
      case 2 :
         if( x>0.5 && y<0.5 && z<0.5 ) 
            // bf = (2.0*x-1.0)*(1.0-2.0*y)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 :
                  result =  2.0*(1.0-2.0*y)*(1.0-2.0*z) ;
                  break ;
               case 1 :
                  result = -2.0*(2.0*x-1.0)*(1.0-2.0*z) ;
                  break ;
               case 2 :
                  result = -2.0*(2.0*x-1.0)*(1.0-2.0*y) ;
                  break ;
            }
         break ;
      case 3 :
         if( x>0.5 && z<0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y*(2.0*x-1.0)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 :
                  result =  4.0*y*(1.0-2.0*z) ;
                  break ;
               case 1 :
                  result =  2.0*(2.0*x-1.0)*(1.0-2.0*z) ;
                  break ;
               case 2 :
                  result = -4.0*y*(2.0*x-1.0) ;
                  break ;
            }
            else // bf = 2.0*(1.0-y)*(2.0*x-1.0)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 :
                  result =  4.0*(1.0-y)*(1.0-2.0*z) ;
                  break ;
               case 1 :
                  result = -2.0*(2.0*x-1.0)*(1.0-2.0*z) ;
                  break ;
               case 2 :
                  result = -4.0*(1.0-y)*(2.0*x-1.0) ;
                  break ;
            }
         }
         break ;
      case 4 :
         if( z<0.5 )
         {
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 4.0*x*y*(1.0-2.0*z)
                  switch( a )
                  {
                     case 0 :
                        result = 4.0*y*(1.0-2.0*z) ;
                        break ;
                     case 1 :
                        result = 4.0*x*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = -8.0*x*y ;
                        break ;
                  }
               else // bf = 4.0*x*(1.0-y)*(1.0-2.0*z)
                  switch( a )
                  {
                     case 0 :
                        result = 4.0*(1.0-y)*(1.0-2.0*z) ;
                        break ;
                     case 1 :
                        result = -4.0*x*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = -8.0*x*(1.0-y) ;
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 4.0*y*(1.0-x)*(1.0-2.0*z) ;
                  switch( a )
                  {
                     case 0 :
                        result = -4.0*y*(1.0-2.0*z) ;
                        break ;
                     case 1 :
                        result =  4.0*(1.0-x)*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = -8.0*y*(1.0-x) ;
                        break ;
                  }
               else // bf = 4.0*(1.0-x)*(1.0-y)*(1.0-2.0*z) ;
                  switch( a )
                  {
                     case 0 :
                        result = -4.0*(1.0-y)*(1.0-2.0*z) ;
                        break ;
                     case 1 :
                        result = -4.0*(1.0-x)*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = -8.0*(1.0-x)*(1.0-y) ;
                        break ;
                  }
            }
         }
         break ;
      case 5 :
         if( x<0.5 && z<0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y*(1.0-2.0*x)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 :
                     result = -4.0*y*(1.0-2.0*z) ;
                     break ;
                  case 1 :
                     result =  2.0*(1.0-2.0*x)*(1.0-2.0*z) ;
                     break ;
                  case 2 :
                     result = -4.0*y*(1.0-2.0*x) ;
                     break ;
               }
            else //       result = 2.0*(1.0-y)*(1.0-2.0*x)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 :
                     result = -4.0*(1.0-y)*(1.0-2.0*z) ;
                     break ;
                  case 1 :
                     result = -2.0*(1.0-2.0*x)*(1.0-2.0*z) ;
                     break ;
                  case 2 :
                     result = -4.0*(1.0-y)*(1.0-2.0*x) ;
                     break ;
               }
         }
         break ;
      case 6 :
         if( x<0.5 && y>0.5 && z<0.5 ) 
            // bf = (1.0-2.0*x)*(2.0*y-1.0)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 :
                  result = -2.0*(2.0*y-1.0)*(1.0-2.0*z) ;
                  break ;
               case 1 :
                  result =  2.0*(1.0-2.0*x)*(1.0-2.0*z) ;
                  break ;
               case 2 :
                  result = -2.0*(1.0-2.0*x)*(2.0*y-1.0) ;
                  break ;
            }
         break ;
      case 7 :
         if( y>0.5 && z<0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x*(2.0*y-1.0)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 :
                     result =  2.0*(2.0*y-1.0)*(1.0-2.0*z) ;
                     break ;
                  case 1 :
                     result =  4.0*x*(1.0-2.0*z) ;
                     break ;
                  case 2 :
                     result = -4.0*x*(2.0*y-1.0) ;
                     break ;
               }
            else  // bf = 2.0*(1.0-x)*(2.0*y-1.0)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 :
                     result = -2.0*(2.0*y-1.0)*(1.0-2.0*z) ;
                     break ;
                  case 1 :
                     result =  4.0*(1.0-x)*(1.0-2.0*z) ;
                     break ;
                  case 2 :
                     result = -4.0*(1.0-x)*(2.0*y-1.0) ;
                     break ;
               }            
         }
         break ;
      case 8 :
         if( x>0.5 && y>0.5 && z<0.5 ) 
            // bf = (2.0*x-1.0)*(2.0*y-1.0)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 :
                  result =  2.0*(2.0*y-1.0)*(1.0-2.0*z) ;
                  break ;
               case 1 :
                  result =  2.0*(2.0*x-1.0)*(1.0-2.0*z) ;
                  break ;
               case 2 :
                  result = -2.0*(2.0*x-1.0)*(2.0*y-1.0) ;
                  break ;
            }
         break ;
      case 9 :
         if( y<0.5 && x<0.5 )
         {
            if( z<0.5 ) // bf = 2.0*z*(1.0-2.0*y)*(1.0-2.0*x) ;
               switch( a )
               {
                  case 0 :
                     result = -4.0*z*(1.0-2.0*y) ;
                     break ;
                  case 1 :
                     result = -4.0*z*(1.0-2.0*x) ;
                     break ;
                  case 2 :
                     result = 2.0*(1.0-2.0*y)*(1.0-2.0*x) ;
                     break ;
               }
            else // bf = 2.0*(1.0-z)*(1.0-2.0*y)*(1.0-2.0*x) ;
               switch( a )
               {
                  case 0 :
                     result = -4.0*(1.0-z)*(1.0-2.0*y) ;
                     break ;
                  case 1 :
                     result = -4.0*(1.0-z)*(1.0-2.0*x) ;
                     break ;
                  case 2 :
                     result = -2.0*(1.0-2.0*y)*(1.0-2.0*x) ;
                     break ;
               }
         }
         break ;
      case 10 :
         if( y<0.5 )
         {
            if( z<0.5 )
            {
               if( x<0.5 ) // bf = 4.0*z*x*(1.0-2.0*y) ;
                  switch( a )
                  {
                     case 0 :
                        result =  4.0*z*(1.0-2.0*y) ;
                        break ;
                     case 1 :
                        result = -8.0*z*x ;
                        break ;
                     case 2 :
                        result = 4.0*x*(1.0-2.0*y) ;
                        break ;
                  }
               else // result = 4.0*z*(1.0-x)*(1.0-2.0*y) ;
                  switch( a )
                  {
                     case 0 :
                        result = -4.0*z*(1.0-2.0*y) ;
                        break ;
                     case 1 :
                        result = -8.0*z*(1.0-x) ;
                        break ;
                     case 2 :
                        result = 4.0*(1.0-x)*(1.0-2.0*y) ;
                        break ;
                  }
            }
            else
            {
               if( x<0.5 ) // bf = 4.0*x*(1.0-z)*(1.0-2.0*y) ;
                  switch( a )
                  {
                     case 0 :
                        result =  4.0*(1.0-z)*(1.0-2.0*y) ;
                        break ;
                     case 1 :
                        result = -8.0*x*(1.0-z) ;
                        break ;
                     case 2 :
                        result = -4.0*x*(1.0-2.0*y) ;
                        break ;
                  }
               else // bf = 4.0*(1.0-z)*(1.0-x)*(1.0-2.0*y) ;
                  switch( a )
                  {
                     case 0 :
                        result = -4.0*(1.0-z)*(1.0-2.0*y) ;
                        break ;
                     case 1 :
                        result = -8.0*(1.0-z)*(1.0-x) ;
                        break ;
                     case 2 :
                        result = -4.0*(1.0-x)*(1.0-2.0*y) ;
                        break ;
                  }
            }
         }
         break ;
      case 11 :
         if( y<0.5 && x>0.5 )
         {
            if( z<0.5 ) // bf = 2.0*z*(1.0-2.0*y)*(2.0*x-1.0) ;
               switch( a )
               {
                  case 0 :
                     result =  4.0*z*(1.0-2.0*y) ;
                     break ;
                  case 1 :
                     result = -4.0*z*(2.0*x-1.0) ;
                     break ;
                  case 2 :
                     result = 2.0*(1.0-2.0*y)*(2.0*x-1.0) ;
                     break ;
               }
            else  // bf = 2.0*(1.0-z)*(1.0-2.0*y)*(2.0*x-1.0) ;
               switch( a )
               {
                  case 0 :
                     result =  4.0*(1.0-z)*(1.0-2.0*y) ;
                     break ;
                  case 1 :
                     result = -4.0*(1.0-z)*(2.0*x-1.0) ;
                     break ;
                  case 2 :
                     result = -2.0*(1.0-2.0*y)*(2.0*x-1.0) ;
                     break ;
               }
         }
         break ;
      case 12 :
         if( x>0.5 )
         {
            if( z<0.5 )
            {
               if( y<0.5 ) // bf = 4.0*z*y*(2.0*x-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result = 8.0*z*y ;
                        break ;
                     case 1 :
                        result = 4.0*z*(2.0*x-1.0) ;
                        break ;
                     case 2 :
                        result = 4.0*y*(2.0*x-1.0) ;
                        break ;
                  }
               else // bf = 4.0*z*(1.0-y)*(2.0*x-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result =  8.0*z*(1.0-y) ;
                        break ;
                     case 1 :
                        result = -4.0*z*(2.0*x-1.0) ;
                        break ;
                     case 2 :
                        result =  4.0*(1.0-y)*(2.0*x-1.0) ;
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 4.0*y*(1.0-z)*(2.0*x-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result =  8.0*y*(1.0-z) ;
                        break ;
                     case 1 :
                        result =  4.0*(1.0-z)*(2.0*x-1.0) ;
                        break ;
                     case 2 :
                        result = -4.0*y*(2.0*x-1.0) ;
                        break ;
                  }
               else // bf = 4.0*(1.0-z)*(1.0-y)*(2.0*x-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result =  8.0*(1.0-z)*(1.0-y) ;
                        break ;
                     case 1 :
                        result = -4.0*(1.0-z)*(2.0*x-1.0) ;
                        break ;
                     case 2 :
                        result = -4.0*(1.0-y)*(2.0*x-1.0) ;
                        break ;
                  }
            }
         }
         break ;
      case 13 :
         if( z<0.5 )
         {
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 8.0*x*y*z ;
                  switch( a )
                  {
                     case 0 :
                        result = 8.0*y*z ;
                        break ;
                     case 1 :
                        result = 8.0*x*z ;
                        break ;
                     case 2 :
                        result = 8.0*x*y ;
                        break ;
                  }
               else  // result = 8.0*x*(1.0-y)*z ;
                  switch( a )
                  {
                     case 0 :
                        result =  8.0*(1.0-y)*z ;
                        break ;
                     case 1 :
                        result = -8.0*x*z ;
                        break ;
                     case 2 :
                        result =  8.0*x*(1.0-y) ;
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 8.0*y*(1.0-x)*z ;
                  switch( a )
                  {
                     case 0 :
                        result = -8.0*y*z ;
                        break ;
                     case 1 :
                        result =  8.0*(1.0-x)*z ;
                        break ;
                     case 2 :
                        result =  8.0*y*(1.0-x) ;
                        break ;
                  }
               else // bf = 8.0*(1.0-x)*(1.0-y)*z ;
                  switch( a )
                  {
                     case 0 :
                        result = -8.0*(1.0-y)*z ;
                        break ;
                     case 1 :
                        result = -8.0*(1.0-x)*z ;
                        break ;
                     case 2 :
                        result =  8.0*(1.0-x)*(1.0-y) ;
                        break ;
                  }
            }
         }
         else
         {
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 8.0*x*y*(1.0-z) ;
                  switch( a )
                  {
                     case 0 :
                        result =  8.0*y*(1.0-z) ;
                        break ;
                     case 1 :
                        result =  8.0*x*(1.0-z) ;
                        break ;
                     case 2 :
                        result = -8.0*x*y ;
                        break ;
                  }
               else // bf = 8.0*x*(1.0-y)*(1.0-z) ;
                  switch( a )
                  {
                     case 0 :
                        result =  8.0*(1.0-y)*(1.0-z) ;
                        break ;
                     case 1 :
                        result = -8.0*x*(1.0-z) ;
                        break ;
                     case 2 :
                        result = -8.0*x*(1.0-y) ;
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 8.0*y*(1.0-x)*(1.0-z) ;
                  switch( a )
                  {
                     case 0 :
                        result = -8.0*y*(1.0-z) ;
                        break ;
                     case 1 :
                        result =  8.0*(1.0-x)*(1.0-z) ;
                        break ;
                     case 2 :
                        result = -8.0*y*(1.0-x) ;
                        break ;
                  }
               else  // bf = 8.0*(1.0-x)*(1.0-y)*(1.0-z) ;
                  switch( a )
                  {
                     case 0 :
                        result = -8.0*(1.0-y)*(1.0-z) ;
                        break ;
                     case 1 :
                        result = -8.0*(1.0-x)*(1.0-z) ;
                        break ;
                     case 2 :
                        result = -8.0*(1.0-x)*(1.0-y) ;
                        break ;
                  }
            }
            
         }
         break ;
      case 14 :
         if( x<0.5 )
         {
            if( z<0.5 )
            {
               if( y<0.5 ) // bf = 4.0*z*y*(1.0-2.0*x) ;
                  switch( a )
                  {
                     case 0 :
                        result = -8.0*z*y ;
                        break ;
                     case 1 :
                        result = 4.0*z*(1.0-2.0*x) ;
                        break ;
                     case 2 :
                        result = 4.0*y*(1.0-2.0*x) ;
                        break ;
                  }
               else  // result = 4.0*z*(1.0-y)*(1.0-2.0*x) ;
                  switch( a )
                  {
                     case 0 :
                        result = -8.0*z*(1.0-y) ;
                        break ;
                     case 1 :
                        result = -4.0*z*(1.0-2.0*x) ;
                        break ;
                     case 2 :
                        result =  4.0*(1.0-y)*(1.0-2.0*x) ;
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 4.0*y*(1.0-z)*(1.0-2.0*x) ;
                  switch( a )
                  {
                     case 0 :
                        result = -8.0*y*(1.0-z) ;
                        break ;
                     case 1 :
                        result =  4.0*(1.0-z)*(1.0-2.0*x) ;
                        break ;
                     case 2 :
                        result = -4.0*y*(1.0-2.0*x) ;
                        break ;
                  }
               else // bf = 4.0*(1.0-z)*(1.0-y)*(1.0-2.0*x) ;
                  switch( a )
                  {
                     case 0 :
                        result = -8.0*(1.0-z)*(1.0-y) ;
                        break ;
                     case 1 :
                        result = -4.0*(1.0-z)*(1.0-2.0*x) ;
                        break ;
                     case 2 :
                        result = -4.0*(1.0-y)*(1.0-2.0*x) ;
                        break ;
                  }
            }
         }
         break ;
      case 15 :
         if( y>0.5 && x<0.5 )
         {
            if( z<0.5 ) // bf = 2.0*z*(2.0*y-1.0)*(1.0-2.0*x) ;
               switch( a )
               {
                  case 0 :
                     result = -4.0*z*(2.0*y-1.0) ;
                     break ;
                  case 1 :
                     result =  4.0*z*(1.0-2.0*x) ;
                     break ;
                  case 2 :
                     result =  2.0*(2.0*y-1.0)*(1.0-2.0*x) ;
                     break ;
               }
            else  // bf = 2.0*(1.0-z)*(2.0*y-1.0)*(1.0-2.0*x) ;
               switch( a )
               {
                  case 0 :
                     result = -4.0*(1.0-z)*(2.0*y-1.0) ;
                     break ;
                  case 1 :
                     result =  4.0*(1.0-z)*(1.0-2.0*x) ;
                     break ;
                  case 2 :
                     result = -2.0*(2.0*y-1.0)*(1.0-2.0*x) ;
                     break ;
               }
         }
         break ;
      case 16 :
         if( y>0.5 )
         {
            if( z<0.5 )
            {
               if( x<0.5 ) // result = 4.0*z*x*(2.0*y-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result = 4.0*z*(2.0*y-1.0) ;
                        break ;
                     case 1 :
                        result = 8.0*z*x ;
                        break ;
                     case 2 :
                        result = 4.0*x*(2.0*y-1.0) ;
                        break ;
                  }
               else // bf = 4.0*z*(1.0-x)*(2.0*y-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result = -4.0*z*(2.0*y-1.0) ;
                        break ;
                     case 1 :
                        result = 8.0*z*(1.0-x) ;
                        break ;
                     case 2 :
                        result = 4.0*(1.0-x)*(2.0*y-1.0) ;
                        break ;
                  }
            }
            else
            {
               if( x<0.5 ) // bf = 4.0*x*(1.0-z)*(2.0*y-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result =  4.0*(1.0-z)*(2.0*y-1.0) ;
                        break ;
                     case 1 :
                        result =  8.0*x*(1.0-z) ;
                        break ;
                     case 2 :
                        result = -4.0*x*(2.0*y-1.0) ;
                        break ;
                  }
               else // bf = 4.0*(1.0-z)*(1.0-x)*(2.0*y-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result = -4.0*(1.0-z)*(2.0*y-1.0) ;
                        break ;
                     case 1 :
                        result =  8.0*(1.0-z)*(1.0-x) ;
                        break ;
                     case 2 :
                        result = -4.0*(1.0-x)*(2.0*y-1.0) ;
                        break ;
                  }
            }
         }
         break ;
      case 17 :
         if( y>0.5 && x>0.5 )
         {
            if( z<0.5 ) // bf = 2.0*z*(2.0*y-1.0)*(2.0*x-1.0) ;
               switch( a )
               {
                  case 0 :
                     result = 4.0*z*(2.0*y-1.0) ;
                     break ;
                  case 1 :
                     result = 4.0*z*(2.0*x-1.0) ;
                     break ;
                  case 2 :
                     result = 2.0*(2.0*y-1.0)*(2.0*x-1.0) ;
                     break ;
               }
            else // bf = 2.0*(1.0-z)*(2.0*y-1.0)*(2.0*x-1.0) ;
               switch( a )
               {
                  case 0 :
                     result =  4.0*(1.0-z)*(2.0*y-1.0) ;
                     break ;
                  case 1 :
                     result =  4.0*(1.0-z)*(2.0*x-1.0) ;
                     break ;
                  case 2 :
                     result = -2.0*(2.0*y-1.0)*(2.0*x-1.0) ;
                     break ;
               }
         }
         break ;
      case 18 :
         if( x<0.5 && y<0.5 && z>0.5 ) 
            // bf = (1.0-2.0*x)*(1.0-2.0*y)*(2.0*z-1.0) ;
            switch( a )
            {
               case 0 :
                  result = -2.0*(1.0-2.0*y)*(2.0*z-1.0) ;
                  break ;
               case 1 :
                  result = -2.0*(1.0-2.0*x)*(2.0*z-1.0) ;
                  break ;
               case 2 :
                  result = 2.0*(1.0-2.0*x)*(1.0-2.0*y) ;
                  break ;
            }
         break ;
      case 19 :
         if( y<0.5 && z>0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x*(1.0-2.0*y)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 :
                     result =  2.0*(1.0-2.0*y)*(2.0*z-1.0) ;
                     break ;
                  case 1 :
                     result = -4.0*x*(2.0*z-1.0) ;
                     break ;
                  case 2 :
                     result =  4.0*x*(1.0-2.0*y) ;
                     break ;
               }
            else // result = 2.0*(1.0-x)*(1.0-2.0*y)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 :
                     result = -2.0*(1.0-2.0*y)*(2.0*z-1.0) ;
                     break ;
                  case 1 :
                     result = -4.0*(1.0-x)*(2.0*z-1.0) ;
                     break ;
                  case 2 :
                     result =  4.0*(1.0-x)*(1.0-2.0*y) ;
                     break ;
               }
         }
         break ;
      case 20 :
         if( x>0.5 && y<0.5 && z>0.5 ) 
            // bf = (2.0*x-1.0)*(1.0-2.0*y)*(2.0*z-1.0) ;
            switch( a )
            {
               case 0 :
                  result =  2.0*(1.0-2.0*y)*(2.0*z-1.0) ;
                  break ;
               case 1 :
                  result = -2.0*(2.0*x-1.0)*(2.0*z-1.0) ;
                  break ;
               case 2 :
                  result =  2.0*(2.0*x-1.0)*(1.0-2.0*y) ;
                  break ;
            }
         break ;
      case 21 :
         if( x>0.5 && z>0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y*(2.0*x-1.0)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 :
                     result = 4.0*y*(2.0*z-1.0) ;
                     break ;
                  case 1 :
                     result = 2.0*(2.0*x-1.0)*(2.0*z-1.0) ;
                     break ;
                  case 2 :
                     result = 4.0*y*(2.0*x-1.0) ;
                     break ;
               }
            else  // bf = 2.0*(1.0-y)*(2.0*x-1.0)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 :
                     result =  4.0*(1.0-y)*(2.0*z-1.0) ;
                     break ;
                  case 1 :
                     result = -2.0*(2.0*x-1.0)*(2.0*z-1.0) ;
                     break ;
                  case 2 :
                     result =  4.0*(1.0-y)*(2.0*x-1.0) ;
                     break ;
               }
         }
         break ;
      case 22 :
         if( z>0.5 )
         {
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 4.0*x*y*(2.0*z-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result = 4.0*y*(2.0*z-1.0) ;
                        break ;
                     case 1 :
                        result = 4.0*x*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result = 8.0*x*y ;
                        break ;
                  }
               else  // bf = 4.0*x*(1.0-y)*(2.0*z-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result =  4.0*(1.0-y)*(2.0*z-1.0) ;
                        break ;
                     case 1 :
                        result = -4.0*x*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result =  8.0*x*(1.0-y) ;
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 4.0*y*(1.0-x)*(2.0*z-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result = -4.0*y*(2.0*z-1.0) ;
                        break ;
                     case 1 :
                        result =  4.0*(1.0-x)*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result =  8.0*y*(1.0-x) ;
                        break ;
                  }
               else // result = 4.0*(1.0-x)*(1.0-y)*(2.0*z-1.0) ;
                  switch( a )
                  {
                     case 0 :
                        result = -4.0*(1.0-y)*(2.0*z-1.0) ;
                        break ;
                     case 1 :
                        result = -4.0*(1.0-x)*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result =  8.0*(1.0-x)*(1.0-y) ;
                        break ;
                  }
            }
         }
         break ;
      case 23 :
         if( x<0.5 && z>0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y*(1.0-2.0*x)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 :
                     result = -4.0*y*(2.0*z-1.0) ;
                     break ;
                  case 1 :
                     result =  2.0*(1.0-2.0*x)*(2.0*z-1.0) ;
                     break ;
                  case 2 :
                     result =  4.0*y*(1.0-2.0*x) ;
                     break ;
               }
            else // bf = 2.0*(1.0-y)*(1.0-2.0*x)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 :
                     result = -4.0*(1.0-y)*(2.0*z-1.0) ;
                     break ;
                  case 1 :
                     result = -2.0*(1.0-2.0*x)*(2.0*z-1.0) ;
                     break ;
                  case 2 :
                     result =  4.0*(1.0-y)*(1.0-2.0*x) ;
                     break ;
               }
         }
         break ;
      case 24 :
         if( x<0.5 && y>0.5 && z>0.5 ) 
            // bf = (1.0-2.0*x)*(2.0*y-1.0)*(2.0*z-1.0) ;
            switch( a )
            {
               case 0 :
                  result = -2.0*(2.0*y-1.0)*(2.0*z-1.0) ;
                  break ;
               case 1 :
                  result =  2.0*(1.0-2.0*x)*(2.0*z-1.0) ;
                  break ;
               case 2 :
                  result =  2.0*(1.0-2.0*x)*(2.0*y-1.0) ;
                  break ;
            }
         break ;
      case 25 :
         if( y>0.5 && z>0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x*(2.0*y-1.0)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 :
                     result = 2.0*(2.0*y-1.0)*(2.0*z-1.0) ;
                     break ;
                  case 1 :
                     result = 4.0*x*(2.0*z-1.0) ;
                     break ;
                  case 2 :
                     result = 4.0*x*(2.0*y-1.0) ;
                     break ;
               }
            else // bf = 2.0*(1.0-x)*(2.0*y-1.0)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 :
                     result = -2.0*(2.0*y-1.0)*(2.0*z-1.0) ;
                     break ;
                  case 1 :
                     result = 4.0*(1.0-x)*(2.0*z-1.0) ;
                     break ;
                  case 2 :
                     result = 4.0*(1.0-x)*(2.0*y-1.0) ;
                     break ;
               }            
         }
         break ;
      case 26 :
         if( x>0.5 && y>0.5 && z>0.5 ) 
            // bf = (2.0*x-1.0)*(2.0*y-1.0)*(2.0*z-1.0) ;
            switch( a )
            {
               case 0 :
                  result = 2.0*(2.0*y-1.0)*(2.0*z-1.0) ;
                  break ;
               case 1 :
                  result = 2.0*(2.0*x-1.0)*(2.0*z-1.0) ;
                  break ;
               case 2 :
                  result = 2.0*(2.0*x-1.0)*(2.0*y-1.0) ;
                  break ;
            }
         break ; 
   }
   return( result ) ;
}


//----------------------------------------------------------------------
double
PDE_3D_Q1isoQ2_27nodes:: d2N_local( size_t node,
                                    size_t a, size_t b,
                                    GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_Q1isoQ2_27nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double result = 0.0 ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   switch( node )
   {
      case 0 :
         if( x<0.5 && y<0.5 && z<0.5 )
            // bf=(1.0-2.0*x)*(1.0-2.0*y)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 : // dbf = -2.0*(1.0-2.0*y)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 1 :
                        result = 4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = 4.0*(1.0-2.0*y) ;
                        break ;
                  }
                  break ;
               case 1 : // dbf = -2.0*(1.0-2.0*x)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 0 :
                        result = 4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = 4.0*(1.0-2.0*x)  ;
                        break ;
                  }
                  break ;
               case 2 : // dbf = -2.0*(1.0-2.0*x)*(1.0-2.0*y) ;
                  switch( b )
                  {
                     case 0 :
                        result = 4.0*(1.0-2.0*y) ;
                        break ;
                     case 1 :
                        result = 4.0*(1.0-2.0*x) ;
                        break ;
                  }
                  break ;
            }
         break ;
      case 1 :
         if( y<0.5 && z<0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x*(1.0-2.0*y)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 : // dbf = 2.0*(1.0-2.0*y)*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 1 :
                           result = -4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result = -4.0*(1.0-2.0*y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -4.0*x*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result =  8.0*x ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -4.0*x*(1.0-2.0*y) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(1.0-2.0*y) ;
                           break ;
                        case 1 :
                           result =  8.0*x ;
                           break ;
                     }
                     break ;
               }
            else // bf = 2.0*(1.0-x)*(1.0-2.0*y)*(1.0-2.0*z)
               switch( a )
               {
                  case 0 : // dbf = -2.0*(1.0-2.0*y)*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 1 :
                           result = 4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result = 4.0*(1.0-2.0*y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -4.0*(1.0-x)*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 0 :
                           result = 4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result = 8.0*(1.0-x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -4.0*(1.0-x)*(1.0-2.0*y) ;
                     switch( b )
                     {
                        case 0 :
                           result = 4.0*(1.0-2.0*y) ;
                           break ;
                        case 1 :
                           result = 8.0*(1.0-x) ;
                           break ;
                     }
                     break ;
               }
         }
         break ;
      case 2 :
         if( x>0.5 && y<0.5 && z<0.5 ) 
            // bf = (2.0*x-1.0)*(1.0-2.0*y)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 : // dbf =  2.0*(1.0-2.0*y)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 1 :
                        result =  -4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result =  -4.0*(1.0-2.0*y) ;
                        break ;
                  }
                  break ;
               case 1 : // dbf = -2.0*(2.0*x-1.0)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 0 :
                        result = -4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result =  4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
               case 2 : // dbf = -2.0*(2.0*x-1.0)*(1.0-2.0*y) ;
                  switch( b )
                  {
                     case 0 :
                        result = -4.0*(1.0-2.0*y) ;
                        break ;
                     case 1 :
                        result =  4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
            }
         break ;
      case 3 :
         if( x>0.5 && z<0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y*(2.0*x-1.0)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 : // dbf =  4.0*y*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 1 :
                        result =  4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = -8.0*y ;
                        break ;
                  }
                  break ;
               case 1 : // dbf =  2.0*(2.0*x-1.0)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 0 :
                        result =  4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = -4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
               case 2 : // dbf = -4.0*y*(2.0*x-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result = -8.0*y ;
                        break ;
                     case 1 :
                        result = -4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
            }
            else // bf = 2.0*(1.0-y)*(2.0*x-1.0)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 : // dbf =  4.0*(1.0-y)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 1 :
                        result =  -4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result =  -8.0*(1.0-y) ;
                        break ;
                  }
                  break ;
               case 1 : // dbf = -2.0*(2.0*x-1.0)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 0 :
                        result = -4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result =  4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
               case 2 : // dbf = -4.0*(1.0-y)*(2.0*x-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result = -8.0*(1.0-y) ;
                        break ;
                     case 1 :
                        result =  4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
            }
         }
         break ;
      case 4 :
         if( z<0.5 )
         {
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 4.0*x*y*(1.0-2.0*z)
                  switch( a )
                  {
                     case 0 : // dbf = 4.0*y*(1.0-2.0*z) ;
                        switch( b )
                        {
                           case 1 :
                              result =  4.0*(1.0-2.0*z) ;
                              break ;
                           case 2 :
                              result = -8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = 4.0*x*(1.0-2.0*z) ;
                        switch( b )
                        {
                           case 0 :
                              result =  4.0*(1.0-2.0*z) ;
                              break ;
                           case 2 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -8.0*x*y ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*y ;
                              break ;
                           case 1 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 4.0*x*(1.0-y)*(1.0-2.0*z)
                  switch( a )
                  {
                     case 0 : // dbf = 4.0*(1.0-y)*(1.0-2.0*z) ;
                        switch( b )
                        {
                           case 1 :
                              result = -4.0*(1.0-2.0*z) ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -4.0*x*(1.0-2.0*z) ;
                        switch( b )
                        {
                           case 0 :
                              result = -4.0*(1.0-2.0*z) ;
                              break ;
                           case 2 :
                              result =  8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -8.0*x*(1.0-y) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result =  8.0*x ;
                              break ;
                        }
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 4.0*y*(1.0-x)*(1.0-2.0*z) ;
                  switch( a )
                  {
                     case 0 : // dbf = -4.0*y*(1.0-2.0*z) ;
                        switch( b )
                        {
                           case 1 :
                              result = -4.0*(1.0-2.0*z) ;
                              break ;
                           case 2 :
                              result =  8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf =  4.0*(1.0-x)*(1.0-2.0*z) ;
                        switch( b )
                        {
                           case 0 :
                              result = -4.0*(1.0-2.0*z) ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -8.0*y*(1.0-x) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*y ;
                              break ;
                           case 1 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 4.0*(1.0-x)*(1.0-y)*(1.0-2.0*z) ;
                  switch( a )
                  {
                     case 0 : // dbf = -4.0*(1.0-y)*(1.0-2.0*z) ;
                        switch( b )
                        {
                           case 1 :
                              result = 4.0*(1.0-2.0*z) ;
                              break ;
                           case 2 :
                              result = 8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -4.0*(1.0-x)*(1.0-2.0*z) ;
                        switch( b )
                        {
                           case 0 :
                              result = 4.0*(1.0-2.0*z) ;
                              break ;
                           case 2 :
                              result = 8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -8.0*(1.0-x)*(1.0-y) ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result = 8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
            }
         }
         break ;
      case 5 :
         if( x<0.5 && z<0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y*(1.0-2.0*x)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 : // dbf = -4.0*y*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 1 :
                           result = -4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result =  8.0*y ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf =  2.0*(1.0-2.0*x)*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result = -4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -4.0*y*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result =  8.0*y ;
                           break ;
                        case 1 :
                           result = -4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
               }
            else // bf = 2.0*(1.0-y)*(1.0-2.0*x)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 : // dbf = -4.0*(1.0-y)*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 1 :
                           result = 4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result = 8.0*(1.0-y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -2.0*(1.0-2.0*x)*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 0 :
                           result = 4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result = 4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -4.0*(1.0-y)*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result = 8.0*(1.0-y) ;
                           break ;
                        case 1 :
                           result = 4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
               }
         }
         break ;
      case 6 :
         if( x<0.5 && y>0.5 && z<0.5 ) 
            // bf = (1.0-2.0*x)*(2.0*y-1.0)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 : // dbf = -2.0*(2.0*y-1.0)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 1 :
                        result = -4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result =  4.0*(2.0*y-1.0) ;
                        break ;
                  }
                  break ;
               case 1 : // dbf =  2.0*(1.0-2.0*x)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 0 :
                        result = -4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = -4.0*(1.0-2.0*x) ;
                        break ;
                  }
                  break ;
               case 2 : // dbf = -2.0*(1.0-2.0*x)*(2.0*y-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result =  4.0*(2.0*y-1.0) ;
                        break ;
                     case 1 :
                        result = -4.0*(1.0-2.0*x) ;
                        break ;
                  }
                  break ;
            }
         break ;
      case 7 :
         if( y>0.5 && z<0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x*(2.0*y-1.0)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 : // dbf =  2.0*(2.0*y-1.0)*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 1 :
                           result =  4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result = -4.0*(2.0*y-1.0) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf =  4.0*x*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 0 :
                           result =  4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result = -8.0*x ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -4.0*x*(2.0*y-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(2.0*y-1.0) ;
                           break ;
                        case 1 :
                           result = -8.0*x ;
                           break ;
                     }
                     break ;
               }
            else  // bf = 2.0*(1.0-x)*(2.0*y-1.0)*(1.0-2.0*z) ;
               switch( a )
               {
                  case 0 : // dbf = -2.0*(2.0*y-1.0)*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 1 :
                           result = -4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result =  4.0*(2.0*y-1.0) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf =  4.0*(1.0-x)*(1.0-2.0*z) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(1.0-2.0*z) ;
                           break ;
                        case 2 :
                           result = -8.0*(1.0-x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -4.0*(1.0-x)*(2.0*y-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result =  4.0*(2.0*y-1.0) ;
                           break ;
                        case 1 :
                           result = -8.0*(1.0-x) ;
                           break ;
                     }
                     break ;
               }            
         }
         break ;
      case 8 :
         if( x>0.5 && y>0.5 && z<0.5 ) 
            // bf = (2.0*x-1.0)*(2.0*y-1.0)*(1.0-2.0*z) ;
            switch( a )
            {
               case 0 : // dbf =  2.0*(2.0*y-1.0)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 1 :
                        result =  4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = -4.0*(2.0*y-1.0) ;
                        break ;
                  }
                  break ;
               case 1 : // dbf =  2.0*(2.0*x-1.0)*(1.0-2.0*z) ;
                  switch( b )
                  {
                     case 0 :
                        result =  4.0*(1.0-2.0*z) ;
                        break ;
                     case 2 :
                        result = -4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
               case 2 : // dbf = -2.0*(2.0*x-1.0)*(2.0*y-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result = -4.0*(2.0*y-1.0) ;
                        break ;
                     case 1 :
                        result = -4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
            }
         break ;
      case 9 :
         if( y<0.5 && x<0.5 )
         {
            if( z<0.5 ) // bf = 2.0*z*(1.0-2.0*y)*(1.0-2.0*x) ;
               switch( a )
               {
                  case 0 : // dbf = -4.0*z*(1.0-2.0*y) ;
                     switch( b )
                     {
                        case 1 :
                           result =  8.0*z ;
                           break ;
                        case 2 :
                           result = -4.0*(1.0-2.0*y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -4.0*z*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result =  8.0*z ;
                           break ;
                        case 2 :
                           result = -4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = 2.0*(1.0-2.0*y)*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(1.0-2.0*y) ;
                           break ;
                        case 1 :
                           result = -4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
               }
            else // bf = 2.0*(1.0-z)*(1.0-2.0*y)*(1.0-2.0*x) ;
               switch( a )
               {
                  case 0 : // dbf = -4.0*(1.0-z)*(1.0-2.0*y) ;
                     switch( b )
                     {
                        case 1 :
                           result = 8.0*(1.0-z) ;
                           break ;
                        case 2 :
                           result = 4.0*(1.0-2.0*y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -4.0*(1.0-z)*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result = 8.0*(1.0-z) ;
                           break ;
                        case 2 :
                           result = 4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -2.0*(1.0-2.0*y)*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result = 4.0*(1.0-2.0*y) ;
                           break ;
                        case 1 :
                           result = 4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
               }
         }
         break ;
      case 10 :
         if( y<0.5 )
         {
            if( z<0.5 )
            {
               if( x<0.5 ) // bf = 4.0*z*x*(1.0-2.0*y) ;
                  switch( a )
                  {
                     case 0 : // dbf =  4.0*z*(1.0-2.0*y) ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result =  4.0*(1.0-2.0*y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -8.0*z*x ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = 4.0*x*(1.0-2.0*y) ;
                        switch( b )
                        {
                           case 0 :
                              result =  4.0*(1.0-2.0*y) ;
                              break ;
                           case 1 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                  }
               else // result = 4.0*z*(1.0-x)*(1.0-2.0*y) ;
                  switch( a )
                  {
                     case 0 : // dbf = -4.0*z*(1.0-2.0*y) ;
                        switch( b )
                        {
                           case 1 :
                              result =  8.0*z ;
                              break ;
                           case 2 :
                              result = -4.0*(1.0-2.0*y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -8.0*z*(1.0-x) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*z ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = 4.0*(1.0-x)*(1.0-2.0*y) ;
                        switch( b )
                        {
                           case 0 :
                              result = -4.0*(1.0-2.0*y) ;
                              break ;
                           case 1 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
            }
            else
            {
               if( x<0.5 ) // bf = 4.0*x*(1.0-z)*(1.0-2.0*y) ;
                  switch( a )
                  {
                     case 0 : // dbf =  4.0*(1.0-z)*(1.0-2.0*y) ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -4.0*(1.0-2.0*y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -8.0*x*(1.0-z) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result =  8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -4.0*x*(1.0-2.0*y) ;
                        switch( b )
                        {
                           case 0 :
                              result = -4.0*(1.0-2.0*y) ;
                              break ;
                           case 1 :
                              result =  8.0*x ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 4.0*(1.0-z)*(1.0-x)*(1.0-2.0*y) ;
                  switch( a )
                  {
                     case 0 : // dbf = -4.0*(1.0-z)*(1.0-2.0*y) ;
                        switch( b )
                        {
                           case 1 :
                              result = 8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = 4.0*(1.0-2.0*y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -8.0*(1.0-z)*(1.0-x) ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = 8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -4.0*(1.0-x)*(1.0-2.0*y) ;
                        switch( b )
                        {
                           case 0 :
                              result = 4.0*(1.0-2.0*y) ;
                              break ;
                           case 1 :
                              result = 8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
            }
         }
         break ;
      case 11 :
         if( y<0.5 && x>0.5 )
         {
            if( z<0.5 ) // bf = 2.0*z*(1.0-2.0*y)*(2.0*x-1.0) ;
               switch( a )
               {
                  case 0 : // dbf =  4.0*z*(1.0-2.0*y) ;
                     switch( b )
                     {
                        case 1 :
                           result = -8.0*z ;
                           break ;
                        case 2 :
                           result =  4.0*(1.0-2.0*y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -4.0*z*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -8.0*z ;
                           break ;
                        case 2 :
                           result = -4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = 2.0*(1.0-2.0*y)*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result =  4.0*(1.0-2.0*y) ;
                           break ;
                        case 1 :
                           result = -4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
               }
            else  // bf = 2.0*(1.0-z)*(1.0-2.0*y)*(2.0*x-1.0) ;
               switch( a )
               {
                  case 0 : // dbf =  4.0*(1.0-z)*(1.0-2.0*y) ;
                     switch( b )
                     {
                        case 1 :
                           result = -8.0*(1.0-z) ;
                           break ;
                        case 2 :
                           result = -4.0*(1.0-2.0*y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -4.0*(1.0-z)*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -8.0*(1.0-z) ;
                           break ;
                        case 2 :
                           result =  4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -2.0*(1.0-2.0*y)*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(1.0-2.0*y) ;
                           break ;
                        case 1 :
                           result =  4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
               }
         }
         break ;
      case 12 :
         if( x>0.5 )
         {
            if( z<0.5 )
            {
               if( y<0.5 ) // bf = 4.0*z*y*(2.0*x-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf = 8.0*z*y ;
                        switch( b )
                        {
                           case 1 :
                              result = 8.0*z ;
                              break ;
                           case 2 :
                              result = 8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = 4.0*z*(2.0*x-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*z ;
                              break ;
                           case 2 :
                              result = 4.0*(2.0*x-1.0) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = 4.0*y*(2.0*x-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*y ;
                              break ;
                           case 1 :
                              result = 4.0*(2.0*x-1.0) ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 4.0*z*(1.0-y)*(2.0*x-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf =  8.0*z*(1.0-y) ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result =  8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -4.0*z*(2.0*x-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result = -4.0*(2.0*x-1.0) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf =  4.0*(1.0-y)*(2.0*x-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result = -4.0*(2.0*x-1.0) ;
                              break ;
                        }
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 4.0*y*(1.0-z)*(2.0*x-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf =  8.0*y*(1.0-z) ;
                        switch( b )
                        {
                           case 1 :
                              result =  8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf =  4.0*(1.0-z)*(2.0*x-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -4.0*(2.0*x-1.0) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -4.0*y*(2.0*x-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*y ;
                              break ;
                           case 1 :
                              result = -4.0*(2.0*x-1.0) ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 4.0*(1.0-z)*(1.0-y)*(2.0*x-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf =  8.0*(1.0-z)*(1.0-y) ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -4.0*(1.0-z)*(2.0*x-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result =  4.0*(2.0*x-1.0) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -4.0*(1.0-y)*(2.0*x-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result =  4.0*(2.0*x-1.0) ;
                              break ;
                        }
                        break ;
                  }
            }
         }
         break ;
      case 13 :
         if( z<0.5 )
         {
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 8.0*x*y*z ;
                  switch( a )
                  {
                     case 0 : // dbf = 8.0*y*z ;
                        switch( b )
                        {
                           case 1 :
                              result = 8.0*z ;
                              break ;
                           case 2 :
                              result = 8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = 8.0*x*z ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*z ;
                              break ;
                           case 2 :
                              result = 8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = 8.0*x*y ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*y ;
                              break ;
                           case 1 :
                              result = 8.0*x ;
                              break ;
                        }
                        break ;
                  }
               else  // result = 8.0*x*(1.0-y)*z ;
                  switch( a )
                  {
                     case 0 : // dbf =  8.0*(1.0-y)*z ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result =  8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -8.0*x*z ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf =  8.0*x*(1.0-y) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 8.0*y*(1.0-x)*z ;
                  switch( a )
                  {
                     case 0 : // dbf = -8.0*y*z ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result = -8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf =  8.0*(1.0-x)*z ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result =  8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf =  8.0*y*(1.0-x) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*y ;
                              break ;
                           case 1 :
                              result =  8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 8.0*(1.0-x)*(1.0-y)*z ;
                  switch( a )
                  {
                     case 0 : // dbf = -8.0*(1.0-y)*z ;
                        switch( b )
                        {
                           case 1 :
                              result =  8.0*z ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -8.0*(1.0-x)*z ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*z ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf =  8.0*(1.0-x)*(1.0-y) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
            }
         }
         else
         {
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 8.0*x*y*(1.0-z) ;
                  switch( a )
                  {
                     case 0 : // dbf =  8.0*y*(1.0-z) ;
                        switch( b )
                        {
                           case 1 :
                              result =  8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf =  8.0*x*(1.0-z) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -8.0*x*y ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*y ;
                              break ;
                           case 1 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 8.0*x*(1.0-y)*(1.0-z) ;
                  switch( a )
                  {
                     case 0 : // dbf =  8.0*(1.0-y)*(1.0-z) ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -8.0*x*(1.0-z) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result =  8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -8.0*x*(1.0-y) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result =  8.0*x ;
                              break ;
                        }
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 8.0*y*(1.0-x)*(1.0-z) ;
                  switch( a )
                  {
                     case 0 : // dbf = -8.0*y*(1.0-z) ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result =  8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf =  8.0*(1.0-x)*(1.0-z) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -8.0*y*(1.0-x) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*y ;
                              break ;
                           case 1 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
               else  // bf = 8.0*(1.0-x)*(1.0-y)*(1.0-z) ;
                  switch( a )
                  {
                     case 0 : // dbf = -8.0*(1.0-y)*(1.0-z) ;
                        switch( b )
                        {
                           case 1 :
                              result = 8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = 8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -8.0*(1.0-x)*(1.0-z) ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = 8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -8.0*(1.0-x)*(1.0-y) ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result = 8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
            }
            
         }
         break ;
      case 14 :
         if( x<0.5 )
         {
            if( z<0.5 )
            {
               if( y<0.5 ) // bf = 4.0*z*y*(1.0-2.0*x) ;
                  switch( a )
                  {
                     case 0 : // dbf = -8.0*z*y ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result = -8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = 4.0*z*(1.0-2.0*x) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*z  ;
                              break ;
                           case 2 :
                              result =  4.0*(1.0-2.0*x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = 4.0*y*(1.0-2.0*x) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*y ;
                              break ;
                           case 1 :
                              result =  4.0*(1.0-2.0*x) ;
                              break ;
                        }
                        break ;
                  }
               else  // result = 4.0*z*(1.0-y)*(1.0-2.0*x) ;
                  switch( a )
                  {
                     case 0 : // dbf = -8.0*z*(1.0-y) ;
                        switch( b )
                        {
                           case 1 :
                              result =  8.0*z ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -4.0*z*(1.0-2.0*x) ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*z ;
                              break ;
                           case 2 :
                              result = -4.0*(1.0-2.0*x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf =  4.0*(1.0-y)*(1.0-2.0*x) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result = -4.0*(1.0-2.0*x) ;
                              break ;
                        }
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 4.0*y*(1.0-z)*(1.0-2.0*x) ;
                  switch( a )
                  {
                     case 0 : // dbf = -8.0*y*(1.0-z) ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result =  8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf =  4.0*(1.0-z)*(1.0-2.0*x) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -4.0*(1.0-2.0*x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -4.0*y*(1.0-2.0*x) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*y ;
                              break ;
                           case 1 :
                              result = -4.0*(1.0-2.0*x) ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 4.0*(1.0-z)*(1.0-y)*(1.0-2.0*x) ;
                  switch( a )
                  {
                     case 0 : // dbf = -8.0*(1.0-z)*(1.0-y) ;
                        switch( b )
                        {
                           case 1 :
                              result = 8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = 8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // sbf = -4.0*(1.0-z)*(1.0-2.0*x) ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = 4.0*(1.0-2.0*x) ;
                              break ;
                        }
                        break ;
                     case 2 : //  dbf = -4.0*(1.0-y)*(1.0-2.0*x) ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result = 4.0*(1.0-2.0*x) ;
                              break ;
                        }
                        break ;
                  }
            }
         }
         break ;
      case 15 :
         if( y>0.5 && x<0.5 )
         {
            if( z<0.5 ) // bf = 2.0*z*(2.0*y-1.0)*(1.0-2.0*x) ;
               switch( a )
               {
                  case 0 : // dbf = -4.0*z*(2.0*y-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result = -8.0*z ;
                           break ;
                        case 2 :
                           result = -4.0*(2.0*y-1.0) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf =  4.0*z*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result = -8.0*z ;
                           break ;
                        case 2 :
                           result =  4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf =  2.0*(2.0*y-1.0)*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(2.0*y-1.0) ;
                           break ;
                        case 1 :
                           result =  4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
               }
            else  // bf = 2.0*(1.0-z)*(2.0*y-1.0)*(1.0-2.0*x) ;
               switch( a )
               {
                  case 0 : // dbf = -4.0*(1.0-z)*(2.0*y-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result = -8.0*(1.0-z) ;
                           break ;
                        case 2 :
                           result =  4.0*(2.0*y-1.0) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf =  4.0*(1.0-z)*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result = -8.0*(1.0-z) ;
                           break ;
                        case 2 :
                           result = -4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -2.0*(2.0*y-1.0)*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result =  4.0*(2.0*y-1.0) ;
                           break ;
                        case 1 :
                           result = -4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
               }
         }
         break ;
      case 16 :
         if( y>0.5 )
         {
            if( z<0.5 )
            {
               if( x<0.5 ) // result = 4.0*z*x*(2.0*y-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf = 4.0*z*(2.0*y-1.0) ;
                        switch( b )
                        {
                           case 1 :
                              result = 8.0*z ;
                              break ;
                           case 2 :
                              result = 4.0*(2.0*y-1.0) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = 8.0*z*x ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*z ;
                              break ;
                           case 2 :
                              result = 8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = 4.0*x*(2.0*y-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = 4.0*(2.0*y-1.0) ;
                              break ;
                           case 1 :
                              result = 8.0*x ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 4.0*z*(1.0-x)*(2.0*y-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf = -4.0*z*(2.0*y-1.0) ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result = -4.0*(2.0*y-1.0) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = 8.0*z*(1.0-x) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*z ;
                              break ;
                           case 2 :
                              result =  8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = 4.0*(1.0-x)*(2.0*y-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = -4.0*(2.0*y-1.0) ;
                              break ;
                           case 1 :
                              result =  8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
            }
            else
            {
               if( x<0.5 ) // bf = 4.0*x*(1.0-z)*(2.0*y-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf =  4.0*(1.0-z)*(2.0*y-1.0) ;
                        switch( b )
                        {
                           case 1 :
                              result =  8.0*(1.0-z);
                              break ;
                           case 2 :
                              result = -4.0*(2.0*y-1.0) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf =  8.0*x*(1.0-z) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -4.0*x*(2.0*y-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = -4.0*(2.0*y-1.0) ;
                              break ;
                           case 1 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                  }
               else // bf = 4.0*(1.0-z)*(1.0-x)*(2.0*y-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf = -4.0*(1.0-z)*(2.0*y-1.0) ;
                        switch( b )
                        {
                           case 1 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result =  4.0*(2.0*y-1.0) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf =  8.0*(1.0-z)*(1.0-x) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*(1.0-z) ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = -4.0*(1.0-x)*(2.0*y-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result =  4.0*(2.0*y-1.0) ;
                              break ;
                           case 1 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
            }
         }
         break ;
      case 17 :
         if( y>0.5 && x>0.5 )
         {
            if( z<0.5 ) // bf = 2.0*z*(2.0*y-1.0)*(2.0*x-1.0) ;
               switch( a )
               {
                  case 0 : // dbf = 4.0*z*(2.0*y-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result = 8.0*z ;
                           break ;
                        case 2 :
                           result = 4.0*(2.0*y-1.0) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = 4.0*z*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = 8.0*z ;
                           break ;
                        case 2 :
                           result = 4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = 2.0*(2.0*y-1.0)*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = 4.0*(2.0*y-1.0) ;
                           break ;
                        case 1 :
                           result = 4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
               }
            else // bf = 2.0*(1.0-z)*(2.0*y-1.0)*(2.0*x-1.0) ;
               switch( a )
               {
                  case 0 : // dbf =  4.0*(1.0-z)*(2.0*y-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result =  8.0*(1.0-z) ;
                           break ;
                        case 2 :
                           result = -4.0*(2.0*y-1.0) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf =  4.0*(1.0-z)*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result =  8.0*(1.0-z) ;
                           break ;
                        case 2 :
                           result = -4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = -2.0*(2.0*y-1.0)*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(2.0*y-1.0) ;
                           break ;
                        case 1 :
                           result = -4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
               }
         }
         break ;
      case 18 :
         if( x<0.5 && y<0.5 && z>0.5 ) 
            // bf = (1.0-2.0*x)*(1.0-2.0*y)*(2.0*z-1.0) ;
            switch( a )
            {
               case 0 : // dbf = -2.0*(1.0-2.0*y)*(2.0*z-1.0) ;
                  switch( b )
                  {
                     case 1 :
                        result =  4.0*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result = -4.0*(1.0-2.0*y) ;
                        break ;
                  }
                  break ;
               case 1 : // dbf = -2.0*(1.0-2.0*x)*(2.0*z-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result =  4.0*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result = -4.0*(1.0-2.0*x) ;
                        break ;
                  }
                  break ;
               case 2 : // dbf = 2.0*(1.0-2.0*x)*(1.0-2.0*y) ;
                  switch( b )
                  {
                     case 0 :
                        result = -4.0*(1.0-2.0*y) ;
                        break ;
                     case 1 :
                        result = -4.0*(1.0-2.0*x) ;
                        break ;
                  }
                  break ;
            }
         break ;
      case 19 :
         if( y<0.5 && z>0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x*(1.0-2.0*y)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 : // dbf =  2.0*(1.0-2.0*y)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result = -4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result =  4.0*(1.0-2.0*y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -4.0*x*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = -8.0*x ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf =  4.0*x*(1.0-2.0*y) ;
                     switch( b )
                     {
                        case 0 :
                           result =  4.0*(1.0-2.0*y) ;
                           break ;
                        case 1 :
                           result = -8.0*x ;
                           break ;
                     }
                     break ;
               }
            else // result = 2.0*(1.0-x)*(1.0-2.0*y)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 : // dbf = -2.0*(1.0-2.0*y)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result =  4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = -4.0*(1.0-2.0*y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -4.0*(1.0-x)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result =  4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = -8.0*(1.0-x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf =  4.0*(1.0-x)*(1.0-2.0*y) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(1.0-2.0*y) ;
                           break ;
                        case 1 :
                           result = -8.0*(1.0-x) ;
                           break ;
                     }
                     break ;
               }
         }
         break ;
      case 20 :
         if( x>0.5 && y<0.5 && z>0.5 ) 
            // bf = (2.0*x-1.0)*(1.0-2.0*y)*(2.0*z-1.0) ;
            switch( a )
            {
               case 0 : // dbf =  2.0*(1.0-2.0*y)*(2.0*z-1.0) ;
                  switch( b )
                  {
                     case 1 :
                        result = -4.0*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result =  4.0*(1.0-2.0*y) ;
                        break ;
                  }
                  break ;
               case 1 : // dbf = -2.0*(2.0*x-1.0)*(2.0*z-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result = -4.0*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result = -4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
               case 2 : // dbf =  2.0*(2.0*x-1.0)*(1.0-2.0*y) ;
                  switch( b )
                  {
                     case 0 :
                        result =  4.0*(1.0-2.0*y) ;
                        break ;
                     case 1 :
                        result = -4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
            }
         break ;
      case 21 :
         if( x>0.5 && z>0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y*(2.0*x-1.0)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 : // dbf = 4.0*y*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result = 4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = 8.0*y ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = 2.0*(2.0*x-1.0)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = 4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = 4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = 4.0*y*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = 8.0*y ;
                           break ;
                        case 1 :
                           result = 4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
               }
            else  // bf = 2.0*(1.0-y)*(2.0*x-1.0)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 : // dbf =  4.0*(1.0-y)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result = -4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result =  8.0*(1.0-y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -2.0*(2.0*x-1.0)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = -4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf =  4.0*(1.0-y)*(2.0*x-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result =  8.0*(1.0-y) ;
                           break ;
                        case 1 :
                           result = -4.0*(2.0*x-1.0) ;
                           break ;
                     }
                     break ;
               }
         }
         break ;
      case 22 :
         if( z>0.5 )
         {
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 4.0*x*y*(2.0*z-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf = 4.0*y*(2.0*z-1.0) ;
                        switch( b )
                        {
                           case 1 :
                              result = 4.0*(2.0*z-1.0) ;
                              break ;
                           case 2 :
                              result = 8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = 4.0*x*(2.0*z-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = 4.0*(2.0*z-1.0) ;
                              break ;
                           case 2 :
                              result = 8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf = 8.0*x*y ;
                        switch( b )
                        {
                           case 0 :
                              result = 8.0*y ;
                              break ;
                           case 1 :
                              result = 8.0*x ;
                              break ;
                        }
                        break ;
                  }
               else  // bf = 4.0*x*(1.0-y)*(2.0*z-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf =  4.0*(1.0-y)*(2.0*z-1.0) ;
                        switch( b )
                        {
                           case 1 :
                              result = -4.0*(2.0*z-1.0) ;
                              break ;
                           case 2 :
                              result =  8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -4.0*x*(2.0*z-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = -4.0*(2.0*z-1.0) ;
                              break ;
                           case 2 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf =  8.0*x*(1.0-y) ;
                        switch( b )
                        {
                           case 0 :
                              result =  8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result = -8.0*x ;
                              break ;
                        }
                        break ;
                  }
            }
            else
            {
               if( y<0.5 ) // bf = 4.0*y*(1.0-x)*(2.0*z-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf = -4.0*y*(2.0*z-1.0) ;
                        switch( b )
                        {
                           case 1 :
                              result = -4.0*(2.0*z-1.0) ;
                              break ;
                           case 2 :
                              result = -8.0*y ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf =  4.0*(1.0-x)*(2.0*z-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result = -4.0*(2.0*z-1.0) ;
                              break ;
                           case 2 :
                              result =  8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf =  8.0*y*(1.0-x) ;
                        switch( b )
                        {
                           case 0 :
                              result = -8.0*y ;
                              break ;
                           case 1 :
                              result =  8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
               else // result = 4.0*(1.0-x)*(1.0-y)*(2.0*z-1.0) ;
                  switch( a )
                  {
                     case 0 : // dbf = -4.0*(1.0-y)*(2.0*z-1.0) ;
                        switch( b )
                        {
                           case 1 :
                              result =  4.0*(2.0*z-1.0) ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-y) ;
                              break ;
                        }
                        break ;
                     case 1 : // dbf = -4.0*(1.0-x)*(2.0*z-1.0) ;
                        switch( b )
                        {
                           case 0 :
                              result =  4.0*(2.0*z-1.0) ;
                              break ;
                           case 2 :
                              result = -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                     case 2 : // dbf =  8.0*(1.0-x)*(1.0-y) ;
                        switch( b )
                        {
                           case 0 :
                              result =  -8.0*(1.0-y) ;
                              break ;
                           case 1 :
                              result =  -8.0*(1.0-x) ;
                              break ;
                        }
                        break ;
                  }
            }
         }
         break ;
      case 23 :
         if( x<0.5 && z>0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y*(1.0-2.0*x)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 : // dbf = -4.0*y*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result = -4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = -8.0*y ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf =  2.0*(1.0-2.0*x)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result =  4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf =  4.0*y*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result = -8.0*y ;
                           break ;
                        case 1 :
                           result =  4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
               }
            else // bf = 2.0*(1.0-y)*(1.0-2.0*x)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 : // dbf = -4.0*(1.0-y)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result =  4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = -8.0*(1.0-y) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = -2.0*(1.0-2.0*x)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result =  4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = -4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf =  4.0*(1.0-y)*(1.0-2.0*x) ;
                     switch( b )
                     {
                        case 0 :
                           result = -8.0*(1.0-y) ;
                           break ;
                        case 1 :
                           result = -4.0*(1.0-2.0*x) ;
                           break ;
                     }
                     break ;
               }
         }
         break ;
      case 24 :
         if( x<0.5 && y>0.5 && z>0.5 ) 
            // bf = (1.0-2.0*x)*(2.0*y-1.0)*(2.0*z-1.0) ;
            switch( a )
            {
               case 0 : // dbf = -2.0*(2.0*y-1.0)*(2.0*z-1.0) ;
                  switch( b )
                  {
                     case 1 :
                        result = -4.0*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result = -4.0*(2.0*y-1.0) ;
                        break ;
                  }
                  break ;
               case 1 : // dbf =  2.0*(1.0-2.0*x)*(2.0*z-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result = -4.0*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result =  4.0*(1.0-2.0*x) ;
                        break ;
                  }
                  break ;
               case 2 : // dbf =  2.0*(1.0-2.0*x)*(2.0*y-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result = -4.0*(2.0*y-1.0) ;
                        break ;
                     case 1 :
                        result =  4.0*(1.0-2.0*x) ;
                        break ;
                  }
                  break ;
            }
         break ;
      case 25 :
         if( y>0.5 && z>0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x*(2.0*y-1.0)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 : // dbf = 2.0*(2.0*y-1.0)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result = 4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = 4.0*(2.0*y-1.0) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = 4.0*x*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = 4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = 8.0*x ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = 4.0*x*(2.0*y-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = 4.0*(2.0*y-1.0) ;
                           break ;
                        case 1 :
                           result = 8.0*x ;
                           break ;
                     }
                     break ;
               }
            else // bf = 2.0*(1.0-x)*(2.0*y-1.0)*(2.0*z-1.0) ;
               switch( a )
               {
                  case 0 : // dbf = -2.0*(2.0*y-1.0)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 1 :
                           result = -4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result = -4.0*(2.0*y-1.0) ;
                           break ;
                     }
                     break ;
                  case 1 : // dbf = 4.0*(1.0-x)*(2.0*z-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(2.0*z-1.0) ;
                           break ;
                        case 2 :
                           result =  8.0*(1.0-x) ;
                           break ;
                     }
                     break ;
                  case 2 : // dbf = 4.0*(1.0-x)*(2.0*y-1.0) ;
                     switch( b )
                     {
                        case 0 :
                           result = -4.0*(2.0*y-1.0) ;
                           break ;
                        case 1 :
                           result =  8.0*(1.0-x) ;
                           break ;
                     }
                     break ;
               }            
         }
         break ;
      case 26 :
         if( x>0.5 && y>0.5 && z>0.5 ) 
            // bf = (2.0*x-1.0)*(2.0*y-1.0)*(2.0*z-1.0) ;
            switch( a )
            {
               case 0 : // dbf = 2.0*(2.0*y-1.0)*(2.0*z-1.0) ;
                  switch( b )
                  {
                     case 1 :
                        result = 4.0*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result = 4.0*(2.0*y-1.0) ;
                        break ;
                  }
                  break ;
               case 1 : // dbf = 2.0*(2.0*x-1.0)*(2.0*z-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result = 4.0*(2.0*z-1.0) ;
                        break ;
                     case 2 :
                        result = 4.0*(2.0*x-1.0) ;
                        break ;
                  }
                  break ;
               case 2 : // dbf = 2.0*(2.0*x-1.0)*(2.0*y-1.0) ;
                  switch( b )
                  {
                     case 0 :
                        result = 4.0*(2.0*y-1.0) ;
                        break ;
                     case 1 :
                        result = 4.0*(2.0*x-1.0) ;
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
PDE_3D_Q1isoQ2_27nodes:: reference_polyhedron_POST(
                               GE_ReferencePolyhedron const* result ) const 
//----------------------------------------------------------------------------
{
   PEL_ASSERT( PDE_ReferenceElement::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceCube::object() ) ;
   return( true ) ;
}

