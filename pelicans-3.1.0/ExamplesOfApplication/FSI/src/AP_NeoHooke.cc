#include <AP_NeoHooke.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <AP_KinematicState.hh>

#include <iostream>

using std::cout ; using std::endl ;

//---------------------------------------------------------------------------
AP_NeoHooke*
AP_NeoHooke:: create( PEL_Object* a_owner,
                      PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NeoHooke:: create" ) ;

   AP_NeoHooke* result = new AP_NeoHooke( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result  ) ;
}

//---------------------------------------------------------------------------
AP_NeoHooke:: AP_NeoHooke( PEL_Object* a_owner, 
                           PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : AP_ConstitutiveLaw( a_owner )
   , C1( exp->double_data( "c_1" ) )
{
}

//---------------------------------------------------------------------------
AP_NeoHooke:: ~AP_NeoHooke( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
AP_NeoHooke:: update_S_dSdE( AP_KinematicState const* st,
                             doubleArray2D& Pio2,
                             doubleArray4D& dPio2dGL ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_NeoHooke:: update_S_dSdE" ) ;
   PEL_CHECK_PRE( update_S_dSdE_PRE( st, Pio2, dPio2dGL ) ) ;
 
   size_t dim = 3 ;

   double jac = st->detF() ;

   doubleArray2D const& rcg = st->rCG() ;
   doubleArray2D const& inv_rcg = st->inv_rCG() ;

   double jac23 = jac*jac ;
   jac23 = PEL::pow( jac23, -1./3. ) ;

   double trace_rcg = 0.0 ;
   for( size_t i=0 ; i<dim ; ++i )
   {
      trace_rcg += rcg( i, i ) ;
   }

   double a0 = 2. * C1 * jac23 * trace_rcg ;
   
   for( size_t i=0 ; i<dim ; ++i )
   {
      for( size_t j=0 ; j<dim ; ++j )
      {
         double xx = - a0 * inv_rcg( i, j ) / 3.0 ;
 
         if( i == j ) 
            xx += 2. * C1 * jac23 ;

         Pio2( i, j ) = xx ;

         for( size_t k=0 ; k<dim ; ++k )
         {
            for( size_t l=0 ; l<dim ; ++l )
            {
               dPio2dGL( i, j, k, l ) = 2. * a0 
                                  * inv_rcg( i, j ) * inv_rcg( k, l )/ 9.
                             + a0 * ( inv_rcg( i, k ) * inv_rcg( j, l )
                                    + inv_rcg( i, l ) * inv_rcg( j, k ) )/ 3. ;
               if( i == j ) 
                  dPio2dGL( i, j, k, l ) -= 4.*C1*jac23 * inv_rcg( k, l )/ 3. ;
               if( k == l )
                  dPio2dGL( i, j, k, l ) -= 4.*C1*jac23 * inv_rcg( i, j )/ 3. ;
            }
	 }
      }
  } 
}
