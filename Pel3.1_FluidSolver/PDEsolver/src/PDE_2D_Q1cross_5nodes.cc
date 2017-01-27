#include <PDE_2D_Q1cross_5nodes.hh>

#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSquare.hh>
#include <iostream>
using std::string ;
using std::endl ;

//----------------------------------------------------------------------
PDE_2D_Q1cross_5nodes* PDE_2D_Q1cross_5nodes::uniqueInstance = new PDE_2D_Q1cross_5nodes() ;
//----------------------------------------------------------------------



//----------------------------------------------------------------------
PDE_2D_Q1cross_5nodes:: ~PDE_2D_Q1cross_5nodes( void )
//----------------------------------------------------------------------
{
   uniqueInstance = 0 ;
}


//----------------------------------------------------------------------
PDE_2D_Q1cross_5nodes:: PDE_2D_Q1cross_5nodes( void )
//----------------------------------------------------------------------
   : PDE_ReferenceElement( "PDE_2D_Q1cross_5nodes", GE_ReferenceSquare::object() )
{
   append_node( GE_Point::create( this, 0.0, 0.0 ) ) ; // case 0
   append_node( GE_Point::create( this, 1.0, 0.0 ) ) ; // case 1
   append_node( GE_Point::create( this, 1.0, 1.0 ) ) ; // case 2
   append_node( GE_Point::create( this, 0.0, 1.0 ) ) ; // case 3
   append_node( GE_Point::create( this, 0.5, 0.5 ) ) ; // case 4
}



//----------------------------------------------------------------------
double
PDE_2D_Q1cross_5nodes:: N_local( size_t node,
                      GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1cross_5nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   double dd1 = x + y - 1.0 ;
   double dd2 = y - x ;

   double result = 0.0 ;

   switch( node )
   {
      case 0 :
    	if ( dd1 <= 0.0 )
	      result = 1.0-x-y ;
          break ;
      case 1 :
        if ( dd2 <= 0.0 )
	      result = x-y ;
          break ;
      case 2 :
        if ( dd1 >= 0.0 )
	      result = x+y-1.0 ;
          break ;
      case 3 :
     	if ( dd2 >= 0.0 )
	      result = y-x ;
          break ;
      case 4 :
     	if ( dd2 <= 0.0 ) {
	       if (dd1 <= 0.0)
	           result =  2.*y;
	       else
	           result = -2.*x;
     	 } else {
      	     if (dd1 <= 0.0)
	          result =  2.*x;
	         else
	           result = -2.*y;
     	 }
         break ;
   }

   return( result ) ;
}



//----------------------------------------------------------------------
double
PDE_2D_Q1cross_5nodes:: dN_local( size_t node,
                       size_t a,
                       GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1cross_5nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   double dd1 = x + y - 1.0 ;
   double dd2 = y - x ;

   double dbf = 0.0 ;

   switch( node )
   {
      case 0 : // result = 1-x-y
     	 if ( dd1 <= 0.0 ) dbf = -1.0 ;
      case 1 : // result = x-y
     	 if ( dd2 <= 0.0 )
     		 switch( a )
     		 {
				 case 0 :
					 dbf =  1.0 ;
					 break ;
				 case 1 :
					 dbf = -1.0;
					 break ;
     		 }
         break ;
      case 2 : // result = x+y-1 ;
    	 if ( dd1 >= 0.0 ) dbf = 1.0 ;
    	 break ;
      case 3 : // result = y-x ;
      	 if ( dd2 >= 0.0 )
      		 switch( a )
      		 {
      		 case 0 :
      			 dbf = -1.0 ;
      			 break ;
      		 case 1 :
      			 dbf =  1.0 ;
      			 break ;
      		 }
		 break ;
	 case 4 : // result =  {2y, -2x, 2x, -2y}
     	 if ( dd2 <= 0.0 ) {
     		 if (dd1 <= 0.0) // result =  2.*y;
          		 switch( a )
          		 {
          		 case 0 :
          			 dbf = 0.0 ;
          			 break ;
          		 case 1 :
          			 dbf = 2.0 ;
          			 break ;
          		 }
     		 else // result = -2.*x;
          		 switch( a )
          		 {
          		 case 0 :
          			 dbf = -2.0 ;
          			 break ;
          		 case 1 :
          			 dbf =  0.0 ;
          			 break ;
          		 }
     	 } else {
     		 if (dd1 <= 0.0) // result =  2.*x;
          		 switch( a )
          		 {
          		 case 0 :
          			 dbf = 2.0 ;
          			 break ;
          		 case 1 :
          			 dbf = 0.0 ;
          			 break ;
          		 }
     		 else // result = -2.*y;
          		 switch( a )
          		 {
          		 case 0 :
          			 dbf = 0.0 ;
          			 break ;
          		 case 1 :
          			 dbf = -2.0 ;
          			 break ;
          		 }
     	 }
   }

   return( dbf ) ;
}



//----------------------------------------------------------------------
double
PDE_2D_Q1cross_5nodes:: d2N_local( size_t node,
                        size_t a, size_t b,
                        GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1cross_5nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double d2bf = 0.0 ;
   return( d2bf ) ;
}
