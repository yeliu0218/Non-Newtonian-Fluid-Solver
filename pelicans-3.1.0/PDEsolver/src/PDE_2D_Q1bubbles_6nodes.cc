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

#include <PDE_2D_Q1bubbles_6nodes.hh>

#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSquare.hh>

//----------------------------------------------------------------------
PDE_2D_Q1bubbles_6nodes*
PDE_2D_Q1bubbles_6nodes::REGISTRATOR = new PDE_2D_Q1bubbles_6nodes() ;
//----------------------------------------------------------------------

//----------------------------------------------------------------------
PDE_2D_Q1bubbles_6nodes:: PDE_2D_Q1bubbles_6nodes( void )
//----------------------------------------------------------------------
   : PDE_ReferenceElement( "PDE_2D_Q1bubbles_6nodes",
                           GE_ReferenceSquare::object() )
{
   append_node( GE_Point::create( this, 0.0, 0.0 ) ) ;
   append_node( GE_Point::create( this, 1.0, 0.0 ) ) ;
   append_node( GE_Point::create( this, 1.0, 1.0 ) ) ;
   append_node( GE_Point::create( this, 0.0, 1.0 ) ) ;
   append_node( GE_Point::create( this, 2.0/3.0, 1.0/3.0 ) ) ;
   append_node( GE_Point::create( this, 1.0/3.0, 2.0/3.0 ) ) ;
}

//----------------------------------------------------------------------
PDE_2D_Q1bubbles_6nodes:: ~PDE_2D_Q1bubbles_6nodes( void )
//----------------------------------------------------------------------
{
   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
double
PDE_2D_Q1bubbles_6nodes:: N_local( size_t node,
                                   GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1bubbles_6nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;

   double result = PEL::bad_double() ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   switch( node )
   {
      case 0 :
         result = (1.0-x)*(1.0-y) ;
         break ;
      case 1 :
         result = x*(1.0-y) ;
         break ;
      case 2 :
         result = x*y ;
         break ;
      case 3 :
         result = (1.0-x)*y ;
         break ;
      case 4 :
         if( x>=y )
         {
            result = 27.0*(1.0-x)*y*(x-y) ;
         }
         else
         {
            result = 0.0 ;
         }
         break ;
      case 5 :
         if( x<=y )
         {
            result = 27.0*x*(1.0-y)*(y-x) ;
         }
         else
         {
            result = 0.0 ;
         }
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_2D_Q1bubbles_6nodes:: dN_local( size_t node,
                                    size_t a,
                                    GE_Point const* pt_ref )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1bubbles_6nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double result = PEL::bad_double() ;
   
   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   
   switch( node )
   {
      case 0 : // bf = (1.0-x)*(1.0-y) ;
         switch( a )
         {
            case 0 :
               result = -(1.0-y) ;
               break ;
            case 1 :
               result = -(1.0-x) ;
               break ;
         }
         break ;
      case 1 : // bf = x*(1.0-y) ;
         switch( a )
         {
            case 0 :
               result = (1.0-y) ;
               break ;
            case 1 :
               result = -x ;
               break ;
         }
         break ;
      case 2 : // bf = x*y ;
         switch( a )
         {
            case 0 :
               result = y ;
               break ;
            case 1 :
               result = x ;
               break ;
         }
         break ;
      case 3 : // bf = (1.0-x)*y  ;
         switch( a )
         {
            case 0 :
               result = -y ;
               break ;
            case 1 :
               result = (1.0-x) ;
               break ;
         }
         break ;
      case 4 : // bf = 27.0*(1.0-x)*y*(x-y) ;
         if( x>=y )
         {
            switch( a )
            {
               case 0 :
                  result = 27.0*y*(1.0-2.0*x+y) ;
                  break ;
               case 1 :
                  result = 27.0*(1.0-x)*(x-2.0*y) ;
                  break ;
            }
         }
         else
         {
            result = 0.0 ;
         }
         break ;
      case 5 : // bf = 27.0*x*(1.0-y)*(y-x) ;
         if( x<=y )
         {
            switch( a )
            {
               case 0 :
                  result = 27.0*(1.0-y)*(-2.0*x+y) ;
                  break ;
               case 1 :
                  result = 27.0*x*(1.0+x-2.0*y) ;
                  break ;
            }
         }
         else
         {
            result = 0.0 ;
         }
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_2D_Q1bubbles_6nodes:: d2N_local( size_t node,
                                     size_t a, size_t b,
                                     GE_Point const* pt_ref )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1bubbles_6nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double result=PEL::bad_double() ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   switch( node )
   {
      case 0 : // bf = (1.0-x)*(1.0-y) ;
         switch( a )
         {
            case 0 : // dbf = -(1.0-y) ;
               switch( b )
               {
                  case 0 :
                     result = 0.0 ;
                     break ;
                  case 1 :
                     result = 1.0 ;
                     break ;
               }
               break ;
            case 1 : // dbf = -(1.0-x) ;
               switch( b )
               {
                  case 0 :
                     result = 1.0 ;
                     break ;
                  case 1 :
                     result = 0.0 ;
                     break ;
               }
               break ;
         }
         break ;
      case 1 : // bf = x*(1.0-y) ;
         switch( a )
         {
            case 0 : // dbf = (1.0-y) ;
               switch( b )
               {
                  case 0 :
                     result = 0.0 ;
                     break ;
                  case 1 :
                     result = -1.0 ;
                     break ;
               }
               break ;
            case 1 : // dbf = -x ;
               switch( b )
               {
                  case 0 :
                     result = -1.0 ;
                     break ;
                  case 1 :
                     result = 0.0 ;
                     break ;
               }
               break ;
         }
         break ;
      case 2 : // bf = x*y ;
         switch( a )
         {
            case 0 : // dbf = y ;
               switch( b )
               {
                  case 0 :
                     result = 0.0 ;
                     break ;
                  case 1 :
                     result = 1.0 ;
                     break ;
               }
               break ;
            case 1 : // dbf = x ;
               switch( b )
               {
                  case 0 :
                     result = 1.0 ;
                     break ;
                  case 1 :
                     result = 0.0 ;
                     break ;
               }
               break ;
         }
         break ;
      case 3 : // bf = (1.0-x)*y  ;
         switch( a )
         {
            case 0 : // dbf = -y ;
               switch( b )
               {
                  case 0 :
                     result = 0.0 ;
                     break ;
                  case 1 :
                     result = -1.0 ;
                     break ;
               }
               break ;
            case 1 : // dbf = (1.0-x) ;
               switch( b )
               {
                  case 0 :
                     result = -1.0 ;
                     break ;
                  case 1 :
                     result = 0.0 ;
                     break ;
               }
               break ;
         }
         break ;
      case 4 : // bf = 27.0*(1.0-x)*y*(x-y) ;
         if( x>=y )
         {
            switch( a )
            {
               case 0 : // dbf = 27.0*y*(1.0-2.0*x+y) ;
                  switch( b )
                  {
                     case 0 :
                        result = -54.0*y ;
                        break ;
                     case 1 :
                        result = 27.0*(1.0-2.0*x+2.0*y ) ;
                        break ;                  
                  }
                  break ;
               case 1 : // dbf = 27.0*(1.0-x)*(x-2.0*y) ;
                  switch( b )
                  {
                     case 0 :
                        result = 27.0*(1.0-2.0*x+2.0*y ) ;
                        break ;
                     case 1 :
                        result = -54.0*(1.0-x) ;
                        break ;                  
                  }
                  break ;
            }
         }
         else
         {
            result = 0.0 ;
         }
         break ;
      case 5 : // bf = 27.0*x*(1.0-y)*(y-x) ;
         if( x<=y )
         {
            switch( a )
            {
               case 0 : // dbf = 27.0*(1.0-y)*(-2.0*x+y) ;
                  switch( b )
                  {
                     case 0 :
                        result = -54.0*(1.0-y) ;
                        break ;
                     case 1 :
                        result = 27.0*(1.0+2.0*x-2.0*y ) ;
                        break ;                  
                  }
                  break ;
               case 1 : // dbf = 27.0*x*(1.0+x-2.0*y) ;
                  switch( b )
                  {
                     case 0 :
                        result = 27.0*(1.0+2.0*x-2.0*y ) ;
                        break ;
                     case 1 :
                        result = -54.0*x ;
                        break ;                  
                  }
                  break ;
            }
            break ;
         }
         else
         {
            result = 0.0 ;
         }
   } 
   return( result ) ;
}
