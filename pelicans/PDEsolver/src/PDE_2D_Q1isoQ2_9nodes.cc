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

#include <PDE_2D_Q1isoQ2_9nodes.hh>

#include <PEL.hh>
#include <PEL_Error.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSquare.hh>

using std::string ;

PDE_2D_Q1isoQ2_9nodes const*
PDE_2D_Q1isoQ2_9nodes:: REGISTRATOR = new PDE_2D_Q1isoQ2_9nodes()  ;

//----------------------------------------------------------------------
PDE_2D_Q1isoQ2_9nodes:: ~PDE_2D_Q1isoQ2_9nodes( void )
//----------------------------------------------------------------------
{
   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
PDE_2D_Q1isoQ2_9nodes:: PDE_2D_Q1isoQ2_9nodes( void )
//----------------------------------------------------------------------
   : PDE_ReferenceElement( "PDE_2D_Q1isoQ2_9nodes", GE_ReferenceSquare::object() )
{
   append_node( GE_Point::create( this, 0.0, 0.0 ) ) ; // 0
   append_node( GE_Point::create( this, 1.0, 0.0 ) ) ; // 1
   append_node( GE_Point::create( this, 1.0, 1.0 ) ) ; // 2
   append_node( GE_Point::create( this, 0.0, 1.0 ) ) ; // 3
   append_node( GE_Point::create( this, 0.5, 0.0 ) ) ; // 4
   append_node( GE_Point::create( this, 1.0, 0.5 ) ) ; // 5
   append_node( GE_Point::create( this, 0.5, 1.0 ) ) ; // 6
   append_node( GE_Point::create( this, 0.0, 0.5 ) ) ; // 7
   append_node( GE_Point::create( this, 0.5, 0.5 ) ) ; // 8
}

//----------------------------------------------------------------------
double
PDE_2D_Q1isoQ2_9nodes:: N_local( size_t node,
                                 GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1isoQ2_9nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;
   
   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   double result = 0.0 ;
   switch( node )
   {
      case 0 :
         if( x<0.5 && y<0.5 ) result = (1.0-2.0*x) * (1.0-2.0*y) ;
         break ;
      case 1 :
         if( x>0.5 && y<0.5 ) result = (2.0*x-1.0) * (1.0-2.0*y) ;
         break ;
      case 2 :
         if( x>0.5 && y>0.5 ) result = (2.0*x-1.0) * (2.0*y-1.0) ;
         break ;
      case 3 :
         if( x<0.5 && y>0.5 ) result = (1.0-2.0*x) * (2.0*y-1.0) ;
         break ;
      case 4 :
         if( y<0.5 )
         {
            if( x<0.5 ) result = 2.0*x       * (1.0-2.0*y) ;
            else        result = 2.0*(1.0-x) * (1.0-2.0*y) ;
         }
         break ;
      case 5 :
         if( x>0.5 )
         {
            if( y<0.5 ) result = 2.0*y       * (2.0*x-1.0) ;
            else        result = 2.0*(1.0-y) * (2.0*x-1.0) ;
         }
         break ;
      case 6 :
         if( y>0.5 )
         {
            if( x<0.5 ) result = 2.0*x       * (2.0*y-1.0) ;
            else        result = 2.0*(1.0-x) * (2.0*y-1.0) ;
         }
         break ;
      case 7 :
         if( x<0.5 )
         {
            if( y<0.5 ) result = 2.0*y       * (1.0-2.0*x) ;
            else        result = 2.0*(1.0-y) * (1.0-2.0*x) ;
         }
         break ;
      case 8 :
         if( x<0.5 )
         {
            if( y<0.5 ) result = 4.0*x*y ;
            else        result = 4.0*x*(1.0-y) ;
         }
         else
         {
            if( y<0.5 ) result = 4.0*y*(1.0-x) ;
            else        result = 4.0*(1.0-x)*(1.0-y) ;
         }
         break ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_2D_Q1isoQ2_9nodes:: dN_local( size_t node, size_t a,
                                  GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1isoQ2_9nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   double result = 0.0 ;
   switch( node )
   {
      case 0 :
         if( x<0.5 && y<0.5 ) // bf = (1.0-2.0*x)*(1.0-2.0*y)
         {
            if( a==0 ) result = 4.0*y - 2.0 ;
            else       result = 4.0*x - 2.0 ;
         }
         break ;
      case 1 :
         if( x>0.5 && y<0.5 ) // bf = (2.0*x-1.0)*(1.0-2.0*y)
         {
            if( a==0 ) result = 2.0 - 4.0*y ;
            else       result = 2.0 - 4.0*x ;
         }
         break ;
      case 2 :
         if( x>0.5 && y>0.5 ) // bf = (2.0*x-1.0)*(2.0*y-1.0)
         {
            if( a==0 ) result = 4.0*y - 2.0 ;
            else       result = 4.0*x - 2.0 ;
         }
         break ;
      case 3 :
         if( x<0.5 && y>0.5 ) // bf = (1.0-2.0*x)*(2.0*y-1.0)
         {
            if( a==0 ) result = 2.0 - 4.0*y ;
            else       result = 2.0 - 4.0*x ;
         }
         break ;
      case 4 :
         if( y<0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x * (1.0-2.0*y)
            {
               if( a==0 ) result = 2.0 - 4.0*y ;
               else       result = -4.0*x ;
            }
            else  // bf = 2.0*(1.0-x) * (1.0-2.0*y)
            {
               if( a==0 ) result = 4.0*y - 2.0 ;
               else       result = 4.0*x - 4.0 ;
            }
         }
         break ;
      case 5 :
         if( x>0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y * (2.0*x-1.0)
            {
               if( a==0 ) result = 4.0*y ;
               else       result = 4.0*x - 2.0 ; 
            }
            else // bf = 2.0*(1.0-y) * (2.0*x-1.0) ;
            {
               if( a==0 ) result = 4.0 - 4.0*y ;
               else       result = 2.0 - 4.0*x ;
            }
         }
         break ;
      case 6 :
         if( y>0.5 )
         {
            if( x<0.5 ) // bf = 2.0*x * (2.0*y-1.0)
            {
               if( a==0 ) result = 4.0*y - 2.0 ;
               else       result = 4.0*x ;
            }
            else // bf = 2.0*(1.0-x) * (2.0*y-1.0)
            {
               if( a==0 ) result = 2.0 - 4.0*y ;
               else       result = 4.0 - 4.0*x ;
            }
            
         }
         break ;
      case 7 :
         if( x<0.5 )
         {
            if( y<0.5 ) // bf = 2.0*y * (1.0-2.0*x)
            {
               if( a==0 ) result = -4.0*y ;
               else       result = 2.0 - 4.0*x ;
            }
            else // bf = 2.0*(1.0-y) * (1.0-2.0*x)
            {
               if( a==0 ) result = 4.0*y - 4.0 ;
               else       result = 4.0*x - 2.0 ;
            }
            
         }
         break ;
      case 8 :
         if( x<0.5 )
         {
            if( y<0.5 ) // bf = 4.0*x*y
            {
               if( a==0 ) result = 4.0*y ;
               else       result = 4.0*x ;
            }
            else // bf = 4.0*x*(1.0-y)
            {
               if( a==0 ) result = 4.0 - 4.0*y ;
               else       result = -4.0*x ;
            }
         }
         else
         {
            if( y<0.5 ) // bf = 4.0*y*(1.0-x)
            {
               if( a==0 ) result = -4.0*y ;
               else       result = 4.0 - 4.0*x ;
            }
            else // bf = 4.0*(1.0-x)*(1.0-y)
            {
               if( a==0 ) result = 4.0*y - 4.0 ;
               else       result = 4.0*x - 4.0 ;
            }
         }
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_2D_Q1isoQ2_9nodes:: d2N_local( size_t node,
                                   size_t a, size_t b,
                                   GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1isoQ2_9nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double result = 0.0 ;
   
   if( a != b )
   {
      double x = pt_ref->coordinate( 0 ) ;
      double y = pt_ref->coordinate( 1 ) ;

      switch( node )
      {
         case 0 :
            if( x<0.5 && y<0.5 ) // bf = (1.0-2.0*x)*(1.0-2.0*y)
            {
               result = 4.0 ;
            }
            break ;
         case 1 :
            if( x>0.5 && y<0.5 ) // bf = (2.0*x-1.0)*(1.0-2.0*y)
            {
               result = -4.0 ;
            }
            break ;
         case 2 :
            if( x>0.5 && y>0.5 ) // bf = (2.0*x-1.0)*(2.0*y-1.0)
            {
               result = 4.0 ;
            }
            break ;
         case 3 :
            if( x<0.5 && y>0.5 ) // bf = (1.0-2.0*x)*(2.0*y-1.0)
            {
               result = -4.0 ;
            }
            break ;
         case 4 :
            if( y<0.5 )
            {
               if( x<0.5 ) // bf = 2.0*x * (1.0-2.0*y)
               {
                  result = -4.0 ;
               }
               else  // bf = 2.0*(1.0-x) * (1.0-2.0*y)
               {
                  result = 4.0 ;
               }
            }
            break ;
         case 5 :
            if( x>0.5 )
            {
               if( y<0.5 ) // bf = 2.0*y * (2.0*x-1.0)
               {
                  result = 4.0 ; 
               }
               else // bf = 2.0*(1.0-y) * (2.0*x-1.0) ;
               {
                  result = -4.0 ;
               }
            }
            break ;
         case 6 :
            if( y>0.5 )
            {
               if( x<0.5 ) // bf = 2.0*x * (2.0*y-1.0)
               {
                  result = 4.0 ;
               }
               else // bf = 2.0*(1.0-x) * (2.0*y-1.0)
               {
                  result = -4.0 ;
               }

            }
            break ;
         case 7 :
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 2.0*y * (1.0-2.0*x)
               {
                  result = -4.0 ;
               }
               else // bf = 2.0*(1.0-y) * (1.0-2.0*x)
               {
                  result = 4.0 ;
               }

            }
            break ;
         case 8 :
            if( x<0.5 )
            {
               if( y<0.5 ) // bf = 4.0*x*y
               {
                  result = 4.0 ;
               }
               else // bf = 4.0*x*(1.0-y)
               {
                  result = -4.0 ;
               }
            }
            else
            {
               if( y<0.5 ) // bf = 4.0*y*(1.0-x)
               {
                  result = -4.0 ;
               }
               else // bf = 4.0*(1.0-x)*(1.0-y)
               {
                  result = 4.0 ;
               }
            }
            break ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_2D_Q1isoQ2_9nodes:: reference_polyhedron_POST(
                             GE_ReferencePolyhedron const* result ) const 
//----------------------------------------------------------------------
{
   PEL_ASSERT( PDE_ReferenceElement::reference_polyhedron_POST( result ) ) ;
   PEL_ASSERT( result == GE_ReferenceSquare::object() ) ;
   return( true ) ;
}


