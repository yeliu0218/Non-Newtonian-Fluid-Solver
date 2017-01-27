#include <AP_StVenantKirchhoff.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <AP_KinematicState.hh>

#include <iostream>

//---------------------------------------------------------------------------
AP_StVenantKirchhoff*
AP_StVenantKirchhoff:: create( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_StVenantKirchhoff:: create" ) ;

   AP_StVenantKirchhoff* result = new AP_StVenantKirchhoff( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result  ) ;
}
//---------------------------------------------------------------------------
AP_StVenantKirchhoff:: AP_StVenantKirchhoff( PEL_Object* a_owner, 
                                             PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : AP_ConstitutiveLaw( a_owner )
   , YOUNG( exp->double_data( "Young_modulus" ) )
   , POISSON( exp->double_data( "Poisson_coefficient" ) ) 
{
}

//---------------------------------------------------------------------------
AP_StVenantKirchhoff:: ~AP_StVenantKirchhoff( void )
//---------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
AP_StVenantKirchhoff:: update_S_dSdE( AP_KinematicState const* st,
                                      doubleArray2D& Pio2,
                                      doubleArray4D& dPio2dGL ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_StVenantKirchhoff:: update_S_dSdE" ) ;
   PEL_CHECK_PRE( update_S_dSdE_PRE( st, Pio2, dPio2dGL ) ) ;

   size_t dim = 3 ;

   double mu = 0.5 * YOUNG / ( 1. + POISSON ) ;
   double lambda = 2. * mu * POISSON / ( 1. - 2. * POISSON ) ;

   doubleArray2D const& gl = st->GL() ;

   double trace_gl = 0.0 ;
   for( size_t i=0 ; i<dim ; ++i )
   {
      trace_gl += gl( i, i ) ;
   }

   for( size_t i=0 ; i<dim ; ++i )
   {
      for( size_t j=0 ; j<dim ; ++j )
      {
         for( size_t k=0 ; k<dim ; ++k )
         {
            for( size_t l=0 ; l<dim ; ++l )
            {
               double xx = 0.0 ;
               if( i == j && k == l ) xx += lambda ;
               if( i == k && j == l ) xx += mu ;
               if( i == l && j == k ) xx += mu ;
               dPio2dGL( i, j, k, l ) = xx ;
            }
         }
         double yy = 2. * mu * gl( i, j ) ;
         if( i == j ) yy += lambda * trace_gl ;
         Pio2( i, j ) = yy ;
      }
   }
}
