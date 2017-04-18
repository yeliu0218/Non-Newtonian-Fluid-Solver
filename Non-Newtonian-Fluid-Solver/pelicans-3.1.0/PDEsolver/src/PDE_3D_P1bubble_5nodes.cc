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

#include <PDE_3D_P1bubble_5nodes.hh>

#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceTetrahedron.hh>

//----------------------------------------------------------------------
PDE_3D_P1bubble_5nodes const*
 PDE_3D_P1bubble_5nodes::REGISTRATOR = new PDE_3D_P1bubble_5nodes() ;
//----------------------------------------------------------------------

//----------------------------------------------------------------------
PDE_3D_P1bubble_5nodes:: ~PDE_3D_P1bubble_5nodes( void )
//----------------------------------------------------------------------
{
   REGISTRATOR = 0 ;
}

//----------------------------------------------------------------------
PDE_3D_P1bubble_5nodes:: PDE_3D_P1bubble_5nodes( void )
//----------------------------------------------------------------------
   : PDE_ReferenceElement( "PDE_3D_P1bubble_5nodes",
                           GE_ReferenceTetrahedron::object() )
{
   append_node( GE_Point::create( this, 0.0,  0.0,  0.0  ) ) ;
   append_node( GE_Point::create( this, 1.0,  0.0,  0.0  ) ) ;
   append_node( GE_Point::create( this, 0.0,  1.0,  0.0  ) ) ;
   append_node( GE_Point::create( this, 0.0,  0.0,  1.0  ) ) ;
   append_node( GE_Point::create( this, 0.25, 0.25, 0.25 ) ) ;
}

//----------------------------------------------------------------------
double
PDE_3D_P1bubble_5nodes:: N_local( size_t node,
				  GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_P1bubble_5nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;
   
   double result = PEL::max_double() ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   switch( node )
   {
      case 0 :
         result = 1.0-x-y-z ;
         break ;
      case 1 :
         result = x ;
         break ;
      case 2 :
         result = y ;
         break ;
      case 3 :
         result = z ;
         break ;
      case 4 :
         result = 256.*x*y*z*(1.0-x-y-z) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_3D_P1bubble_5nodes:: dN_local( size_t node, size_t a,
                                   GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_P1bubble_5nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double result = PEL::max_double() ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   switch( node )
   {
      case 0 : // bf = 1.0-x-y-z ;
         switch( a )
         {
            case 0 :
               result = -1.0 ;
               break ;
            case 1 :
               result = -1.0 ;
               break ;
            case 2 :
               result = -1.0 ;
               break ;
         }
         break ;
      case 1 : // bf = x;
         switch( a )
         {
            case 0 :
               result = 1.0 ;
               break ;
            case 1 :
               result = 0.0 ;
               break ;
            case 2 :
               result = 0.0 ;
               break ;
         }
         break ;
      case 2 : // bf = y ;
         switch( a )
         {
            case 0 :
               result = 0.0 ;
               break ;
            case 1 :
               result = 1.0 ;
               break ;
            case 2 :
               result = 0.0 ;
               break ;
         }
         break ;
      case 3 : // bf = z  ;
         switch( a )
         {
            case 0 :
               result = 0.0 ;
               break ;
            case 1 :
               result = 0.0 ;
               break ;
            case 2 :
               result = 1.0 ;
               break ;
         }
         break ;
      case 4 : // bf = 256xyz(1.0-x-y-z)  ;
         switch( a )
         {
            case 0 :
               result =  256.*y*z*(1.0-2.*x-y-z) ;
               break ;
            case 1 :
               result =  256.*x*z*(1.0-x-2.*y-z) ;
               break ;
            case 2 :
               result =  256.*x*y*(1.0-x-y-2.*z) ;
               break ;
         }
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PDE_3D_P1bubble_5nodes:: d2N_local( size_t node,
				    size_t a, size_t b,
				    GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_3D_P1bubble_5nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double result = 0.0 ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double z = pt_ref->coordinate( 2 ) ;

   if( node==4 )
   {
      switch( a ) 
      {
	 case 0 : // dbf =  256yz(1.0-2x-y-z)
	    switch( b )
	    {
	       case 0 :
		  result = -512.*y*z ;
		  break ;
	       case 1 :
		  result = 256.*z*(1.-2.*(x+y)-z) ;
		  break ;
	       case 2 :
		  result = 256.*y*(1.-2.*(x+z)-y) ;
		  break ;
	    }
	    break ;
	 case 1 : // dbf =  256xz(1.0-x-2y-z)
	    switch( b )
	    {
	       case 0 :
		  result = 256.*z*(1.0-2.*(x+y)-z) ;
		  break ;
	       case 1 :
		  result = -512.*x*z ;
		  break ;
	       case 2 :
		  result = 256.*x*(1.0-2.*(z+y)-x) ;
		  break ;
	    }
	    break ;
	 case 2 : // dbf =  256xy(1.0-x-y-2z)
	    switch( b )
	    {
	       case 0 :
		  result = 256.*y*(1.0-2.*(x+z)-y) ;
		  break ;
	       case 1 :
		  result = 256.*x*(1.0-2.*(y+z)-x) ;
		  break ;
	       case 2 :
		  result = -512.*x*y ;
		  break ;
	    }
	    break ;
      }
   }

   return( result ) ;
}
