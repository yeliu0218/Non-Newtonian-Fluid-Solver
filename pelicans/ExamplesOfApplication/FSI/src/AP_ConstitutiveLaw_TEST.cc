#include <AP_ConstitutiveLaw_TEST.hh>

#include <PEL_ModuleExplorer.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>

#include <AP_ConstitutiveLaw.hh>
#include <AP_KinematicState.hh>

#include <ios>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::string ; using std::ostringstream ;

AP_ConstitutiveLaw_TEST const* 
AP_ConstitutiveLaw_TEST::PROTOTYPE = new AP_ConstitutiveLaw_TEST() ;

//---------------------------------------------------------------------------
AP_ConstitutiveLaw_TEST:: AP_ConstitutiveLaw_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "AP_ConstitutiveLaw", "AP_ConstitutiveLaw_TEST" )
{
}

//---------------------------------------------------------------------------
AP_ConstitutiveLaw_TEST:: ~AP_ConstitutiveLaw_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
AP_ConstitutiveLaw_TEST:: process_one_test( PEL_ModuleExplorer const* exp ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_ConstitutiveLaw_TEST:: process_one_test" ) ;

   size_t nb_dims = exp->int_data( "nb_space_dimensions" ) ;
   double d_eps = exp->double_data( "dbl_epsilon" ) ;
   double d_min = exp->double_data( "dbl_minimum" ) ;
   double hh = exp->double_data( "hh" ) ;
   bool eq ;

   AP_KinematicState* st = AP_KinematicState::create( 0, nb_dims ) ;

   PEL_ModuleExplorer* ee = 
                       exp->create_subexplorer( 0, "AP_ConstitutiveLaw" ) ;
   AP_ConstitutiveLaw* law = AP_ConstitutiveLaw::create( 0, ee ) ;
   string nn = ee->string_data( "concrete_name" ) ;
   ee->destroy() ; ee = 0 ;

   doubleArray2D const& g_disp =
                        exp->doubleArray2D_data( "displacement_gradient" ) ;

   doubleArray2D Pio2( 3, 3 ) ;
   doubleArray4D dPio2dGL( 3, 3, 3, 3 ) ;
   st->set_state( g_disp ) ;
   doubleArray2D gl = st->GL() ;
   law->update_S_dSdE( st, Pio2, dPio2dGL ) ;

   doubleArray2D new_Pio2( 3, 3 ) ;
   doubleArray4D dummy( 3, 3, 3, 3 ) ;

   bool ok = true ;
   for( size_t k=0 ; k<nb_dims ; ++k )
   {
      for( size_t l=0 ; l<nb_dims ; ++l )
      {
         doubleArray2D new_g_disp = g_disp ;
         new_g_disp( k, l ) += hh ;
         st->set_state( new_g_disp ) ;
         doubleArray2D new_gl = st->GL() ;
         law->update_S_dSdE( st, new_Pio2, dummy ) ;
         
         for( size_t i=0 ; i<nb_dims ; ++i )
         {
            for( size_t j=0 ; j<nb_dims ; ++j )
            {
               double xx1 = ( new_Pio2(i,j) - Pio2(i,j) ) / hh ;
               double xx2 = 0.0 ;
               for( size_t p=0 ; p<nb_dims ; ++p )
               {
                  for( size_t q=0 ; q<nb_dims ; ++q )
                  {
                     xx2 += dPio2dGL( i, j, p, q ) * 
                            ( new_gl( p, q ) - gl( p, q ) ) / hh ;
                  }
               }
               eq = PEL::double_equality( xx1, xx2, d_eps, d_min ) ;
               if( !eq ) display_error( xx1, xx2 ) ;
               ok = ok && eq ;
            }
         }
      }
   }
   {
      ostringstream mesg ;
      mesg << " update_S_dSdE (" << nb_dims << "D)" ;
      notify_one_test_result( nn + mesg.str(), ok ) ;
   }
   st->destroy() ;
   law->destroy() ;
}

//----------------------------------------------------------------------------
void
AP_ConstitutiveLaw_TEST:: display_error( double xx_1, double xx_2 ) const
//----------------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = cout.flags() ;
   cout.setf( ios_base::uppercase | ios_base::scientific ) ;

   cout << std::setprecision( 10 ) << std::setw( 20 ) << xx_1
        << std::setprecision( 10 ) << std::setw( 20 ) << xx_2
        << endl ;

   cout.flags( original_flags ) ;
}
