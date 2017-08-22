#include <PDE_2D_Q1CR_3nodes.hh>

#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSquare.hh>
#include <iostream>

using std::string ;
using std::endl ;

//----------------------------------------------------------------------
PDE_2D_Q1CR_3nodes* PDE_2D_Q1CR_3nodes::uniqueInstance = new PDE_2D_Q1CR_3nodes() ;
//----------------------------------------------------------------------



//----------------------------------------------------------------------
PDE_2D_Q1CR_3nodes:: ~PDE_2D_Q1CR_3nodes( void )
//----------------------------------------------------------------------
{
   uniqueInstance = 0 ;
}


//----------------------------------------------------------------------
PDE_2D_Q1CR_3nodes:: PDE_2D_Q1CR_3nodes( void )
//----------------------------------------------------------------------
   : PDE_ReferenceElement( "PDE_2D_Q1CR_3nodes", GE_ReferenceSquare::object() )
{
   append_node( GE_Point::create( this, 0.25, 0.25 ) ) ;
   append_node( GE_Point::create( this, 0.75, 0.25 ) ) ;
   append_node( GE_Point::create( this, 0.5,  0.75 ) ) ;
}



//----------------------------------------------------------------------
double
PDE_2D_Q1CR_3nodes:: N_local( size_t node,
                      GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1CR_3nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;

   double bf = PEL::max_double() ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   switch( node )
   {
      case 0 :
            bf = -2.*x-y+1.75 ;
         break ;
      case 1 :
         bf = 2.*x-y-0.25 ;
         break ;
      case 2 :
         bf = 2.*y -0.5 ;
         break ;
   }
  // PEL::out() << "Q1CR: CASE " << node << " (x,y) = (" << x << ", " << y << "): bf =  " << bf << std::endl;

   return( bf ) ;
}



//----------------------------------------------------------------------
double
PDE_2D_Q1CR_3nodes:: dN_local( size_t node,
                       size_t a,
                       GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1CR_3nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double dbf = PEL::max_double() ;

   switch( node )
   {
      case 0 : // bf = -2.*x-y+2.25
         switch( a )
         {
            case 0 :
               dbf = -2. ;
               break ;
            case 1 :
               dbf = -1. ;
               break ;
         }
         break ;
      case 1 : // bf = 2.*x-y-0.25
         switch( a )
         {
            case 0 :
	       dbf = 2. ;
               break ;
            case 1 :
               dbf = -1. ;
               break ;
         }
         break ;
      case 2 : // bf =  2.*y -0.5 ;
         switch( a )
         {
            case 0 :
               dbf = 0. ;
               break ;
            case 1 :
               dbf = 2. ;
               break ;
         }
         break ;
   }

   return( dbf ) ;
}



//----------------------------------------------------------------------
double
PDE_2D_Q1CR_3nodes:: d2N_local( size_t node,
                        size_t a, size_t b,
                        GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q1CR_3nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double d2bf = 0.0 ;

   return( d2bf ) ;
}
