#include <PDE_2D_P1DC_3nodes.hh>
#include <PEL.hh>
#include <GE_Point.hh>
#include <GE_ReferenceTriangle.hh>
using std::string ;

PDE_2D_P1DC_3nodes* 
PDE_2D_P1DC_3nodes::uniqueInstance = new PDE_2D_P1DC_3nodes() ;

// !!! Do not set tolerance less than 1e-6 !!! 
const double tol = 1.e-6;

/* Destructor 
-------------*/
PDE_2D_P1DC_3nodes:: ~PDE_2D_P1DC_3nodes( void )
{
   uniqueInstance = 0 ;
}




/* Constructor without argument 
-------------------------------*/
PDE_2D_P1DC_3nodes:: PDE_2D_P1DC_3nodes( void )
 : PDE_ReferenceElement( "PDE_2D_P1DC_3nodes", GE_ReferenceTriangle::object() )
{
   append_node( GE_Point::create( this, tol, tol ) ) ;
   append_node( GE_Point::create( this, 1.0-tol, tol ) ) ;
   append_node( GE_Point::create( this, tol, 1.0-tol ) ) ;
}




/* Return the basis function at a given point
---------------------------------------------*/
double PDE_2D_P1DC_3nodes:: N_local( size_t node,
                      GE_Point const* pt_ref ) const
{
   PEL_LABEL( "PDE_2D_P1DC_3nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;
   double bf = PEL::max_double() ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   switch( node )
   {
      case 0 :
         bf = 1.0-x-y ;
         break ;
      case 1 :
         bf = x ;
         break ;
      case 2 :
         bf = y ;
         break ;
   }
   return( bf ) ;
}




/* Return the derivative of a basis function at a given point
-------------------------------------------------------------*/
double PDE_2D_P1DC_3nodes:: dN_local( size_t node,
                       size_t a,
                       GE_Point const* pt_ref )  const
{
   PEL_LABEL( "PDE_2D_P1DC_3nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double dbf = PEL::max_double() ;

   switch( node )
   {
      case 0 : // bf = 1.0-x-y ;
         switch( a )
         {
            case 0 :
               dbf = -1.0 ;
               break ;
            case 1 :
               dbf = -1.0 ;
               break ;
         }
         break ;
      case 1 : // bf = x ;
         switch( a )
         {
            case 0 :
               dbf = 1.0 ;
               break ;
            case 1 :
               dbf = 0.0 ;
               break ;
         }
         break ;
      case 2 : // bf = y ;
         switch( a )
         {
            case 0 :
               dbf = 0.0 ;
               break ;
            case 1 :
               dbf = 1.0 ;
               break ;
         }
   }
   return( dbf ) ;
}




/* Return the second derivative of a basis function at a given point
--------------------------------------------------------------------*/
double PDE_2D_P1DC_3nodes:: d2N_local( size_t node,
                        size_t a, size_t b,
                        GE_Point const* pt_ref )  const
{
   PEL_LABEL( "PDE_2D_P1DC_3nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double d2bf=0.0 ;
   return( d2bf ) ;
}
