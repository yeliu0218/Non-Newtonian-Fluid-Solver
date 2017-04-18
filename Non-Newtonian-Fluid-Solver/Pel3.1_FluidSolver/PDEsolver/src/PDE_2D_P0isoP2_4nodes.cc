#include <PDE_2D_P0isoP2_4nodes.hh>
#include <PEL.hh>
#include <GE_Point.hh>
#include <GE_ReferenceTriangle.hh>
using std::string ;


PDE_2D_P0isoP2_4nodes const* 
PDE_2D_P0isoP2_4nodes::REGISTRATOR = new PDE_2D_P0isoP2_4nodes() ;


/* Destructor 
-------------*/
PDE_2D_P0isoP2_4nodes:: ~PDE_2D_P0isoP2_4nodes( void )
{
   REGISTRATOR = 0 ;
   
}




/* Constructor without argument 
-------------------------------*/
PDE_2D_P0isoP2_4nodes:: PDE_2D_P0isoP2_4nodes( void )
   : PDE_ReferenceElement( "PDE_2D_P0isoP2_4nodes", 
                           GE_ReferenceTriangle::object() )
{
   append_node( GE_Point::create( this, 0.1667, 0.1667 ) ) ;
   append_node( GE_Point::create( this, 0.6667, 0.1667 ) ) ;
   append_node( GE_Point::create( this, 0.3333, 0.3333 ) ) ;
   append_node( GE_Point::create( this, 0.1667, 0.6667 ) ) ;
   
}




/* Return the basis function at a given point
---------------------------------------------*/
double PDE_2D_P0isoP2_4nodes:: N_local( size_t node,
                                 GE_Point const* pt_ref ) const
{
   PEL_LABEL( "PDE_2D_P0isoP2_4nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;
   double dd = x + y - 0.5 ;

   double result = 0.0 ;

   switch( node )
   {
      case 0 :
         if ( dd <= 0.0 ) result = 1.0 ;
         break ;
      case 1 :
         if ( x >= 0.5 ) result = 1.0 ;
         break ;
      case 2 :
         if (( dd >= 0.0 )&&( x <= 0.5 )&&( y <= 0.5 )) result = 1.0 ;
         break ;
      case 3 :
         if( y >= 0.5 ) result = 1.0 ;
         break ;
   }
   return( result ) ;
   
}




/* Return the derivative of a basis function at a given point
-------------------------------------------------------------*/
double PDE_2D_P0isoP2_4nodes:: dN_local( size_t node,
                                  size_t a,
                                  GE_Point const* pt_ref )  const
{
   PEL_LABEL( "PDE_2D_P0isoP2_4nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double result = 0.0 ;

   return( result ) ;
   
}




/* Return the second derivative of a basis function at a given point
--------------------------------------------------------------------*/
double PDE_2D_P0isoP2_4nodes:: d2N_local( size_t node,
                                   size_t a, 
                                   size_t b,
                                   GE_Point const* pt_ref )  const
{
   PEL_LABEL( "PDE_2D_P0isoP2_4nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double result = 0.0 ;

   return( result ) ;
   
}
