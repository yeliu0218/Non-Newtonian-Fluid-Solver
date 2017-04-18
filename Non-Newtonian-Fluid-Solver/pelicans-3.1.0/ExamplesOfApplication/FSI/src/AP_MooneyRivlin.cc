#include <AP_MooneyRivlin.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <AP_KinematicState.hh>

#include <iostream>

//---------------------------------------------------------------------------
AP_MooneyRivlin*
AP_MooneyRivlin:: create( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_MooneyRivlin:: create" ) ;

   AP_MooneyRivlin* result = new AP_MooneyRivlin( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result  ) ;
}
//---------------------------------------------------------------------------
AP_MooneyRivlin:: AP_MooneyRivlin( PEL_Object* a_owner, 
                                   PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : AP_ConstitutiveLaw( a_owner )
   , C1( exp->double_data( "c_1" ) )
   , C2( exp->double_data( "c_2" ) )
{
}

//---------------------------------------------------------------------------
AP_MooneyRivlin:: ~AP_MooneyRivlin( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
AP_MooneyRivlin:: update_S_dSdE( AP_KinematicState const* st,
                                 doubleArray2D& Pio2,
                                 doubleArray4D& dPio2dGL ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_MooneyRivlin:: update_S_dSdE" ) ;
   PEL_CHECK_PRE( update_S_dSdE_PRE( st, Pio2, dPio2dGL ) ) ;
 
   size_t dim = 3 ;

   double jac = st->detF() ;

   doubleArray2D const& rcg = st->rCG() ;
   doubleArray2D const& inv_rcg = st->inv_rCG() ;

   double jac23 = jac * jac ;
   jac23 = PEL::pow( jac23, -1./3. ) ;
   double jac43 = jac23 * jac23 ;

   double trace_rcg = 0.0 ;
   double trace_rcg2 = 0.0 ;
   for( size_t i=0 ; i<dim ; ++i )
   {
      trace_rcg += rcg( i, i ) ;

      for( size_t j=0 ; j<dim ; ++j )
      {
         trace_rcg2 += rcg( i, j ) * rcg( j, i ) ; 
      }
   }

   double I_2 = 0.5 * ( trace_rcg * trace_rcg - trace_rcg2 ) ;
   double a0 = 2. * C1 * jac23 * trace_rcg ;
   double a1 = 2. * C2 * jac43 ;
   double a2 = a1 * trace_rcg ;

   for( size_t i=0 ; i<dim ; ++i )
   {
      for( size_t j=0 ; j<dim ; ++j )
      {
         double xx = - a0 * inv_rcg( i, j )/ 3.0 ;
         xx -= a1 * ( rcg( i, j ) + 2. * I_2 * inv_rcg( i, j )/ 3.0 ) ;

         if( i == j )
         { 
            xx += 2. * C1 * jac23 ;
            xx += a2 ;
         }
         Pio2( i, j ) = xx ;         

         for( size_t k=0 ; k<dim ; ++k )
         {
            for( size_t l=0 ; l<dim ; ++l )
            {
               dPio2dGL( i, j, k, l ) = 2. * a0 *
                            inv_rcg( i, j ) * inv_rcg( k, l )/ 9. +
                            a0 * ( inv_rcg( i, k ) * inv_rcg( j, l )
                                 + inv_rcg( i, l ) * inv_rcg( j, k ) )/ 3. ;
               dPio2dGL( i, j, k, l ) += 8. * a1 * I_2 *
                            inv_rcg( i, j ) * inv_rcg( k, l )/ 9. +
                            2. * a1 * I_2 *
                            ( inv_rcg( i, k ) * inv_rcg( j, l ) +
                              inv_rcg( i, l ) * inv_rcg( j, k ) )/ 3. +
                            4. * a1 * 
                            ( rcg( i, j ) * inv_rcg( k, l ) +
                              inv_rcg( i, j ) * rcg( k, l ) )/ 3. ;
               if( i == j )
               { 
                  dPio2dGL( i, j, k, l ) -= 4. * C1 * jac23 * 
                                                 inv_rcg( k, l )/ 3. ;
                  dPio2dGL( i, j, k, l ) -= 4. * a2 * inv_rcg( k, l )/ 3. ;
               }
               if( k == l )
               {
                  dPio2dGL( i, j, k, l ) -= 4. * C1 * jac23 *
                                                 inv_rcg( i, j )/ 3. ;
                  dPio2dGL( i, j, k, l ) -= 4. * a2 * inv_rcg( i, j )/ 3. ;
               }
               if( i == j && k == l )
                  dPio2dGL( i, j, k, l ) += 2. * a1 ;
               if( i == k && j == l )
                  dPio2dGL( i, j, k, l ) -= a1 ;
               if( i == l && j == k )
                  dPio2dGL( i, j, k, l ) -= a1 ;               
            }
	 }
      }
   }   
}

