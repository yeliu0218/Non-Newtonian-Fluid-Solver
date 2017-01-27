#include <AP_KinematicState_TEST.hh>

#include <PEL_ModuleExplorer.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>

#include <AP_KinematicState.hh>

#include <ios>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::string ; using std::ostringstream ;

AP_KinematicState_TEST const* 
AP_KinematicState_TEST::PROTOTYPE = new AP_KinematicState_TEST() ;

//---------------------------------------------------------------------------
AP_KinematicState_TEST:: AP_KinematicState_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "AP_KinematicState", "AP_KinematicState_TEST" )
{
}

//---------------------------------------------------------------------------
AP_KinematicState_TEST:: ~AP_KinematicState_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
AP_KinematicState_TEST:: process_one_test( PEL_ModuleExplorer const* exp ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState_TEST:: process_one_test" ) ;

   size_t nb_dims = exp->int_data( "nb_space_dimensions" ) ;
   D_EPS = exp->double_data( "dbl_epsilon" ) ;
   D_MIN = exp->double_data( "dbl_minimum" ) ;
   bool eq ;

   AP_KinematicState* st = AP_KinematicState::create( 0, nb_dims ) ;

   doubleArray2D const& g_disp =
                        exp->doubleArray2D_data( "displacement_gradient" ) ;

   st->set_state( g_disp ) ;

   doubleArray2D const& rcg = st->rCG() ;
   doubleArray2D const& inv_rcg = st->inv_rCG() ;
   doubleArray2D const& gl = st->GL() ;
   doubleArray2D const& F = st->DG() ;
   doubleArray2D const& inv_F = st->inv_F() ;

   bool ok_0 = ok_rCG_inv_rCG( rcg, inv_rcg, nb_dims ) ;
   {
      ostringstream mesg ;
      mesg << "rCG*inv_RCG=Id (" << nb_dims << "D)" ; 
      notify_one_test_result( mesg.str(), ok_0 ) ;
   }

   if( nb_dims == 2 )
   {
      ok_0 = ok_rCG_inv_rCG( rcg, inv_rcg, 3 ) ;
      notify_one_test_result( "rCG*inv_RCG=Id (3D)", ok_0 ) ;
   }

   bool ok_1 = ok_rCG_GL( rcg, gl, nb_dims ) ;
   {
      ostringstream mesg ;
      mesg << "rCG=Id+2*GL (" << nb_dims << "D)" ; 
      notify_one_test_result( mesg.str(), ok_1 ) ;
   }

   if( nb_dims == 2 )
   {
      ok_1 = ok_rCG_GL( rcg, gl, 3 ) ;
      notify_one_test_result( "rCG=Id+2*GL (3D)", ok_1 ) ;
   }

   bool ok_2 = ok_F_inv_F( F, inv_F, nb_dims ) ;
   {
      ostringstream mesg ;
      mesg << "F*inv_F=Id (" << nb_dims << "D)" ; 
      notify_one_test_result( mesg.str(), ok_2 ) ;
   }

   if( nb_dims == 2 )
   {
      ok_2 = ok_F_inv_F( F, inv_F, 3 ) ;
      notify_one_test_result( "F*inv_F=Id (3D)", ok_2 ) ;
   }   

   if( exp->has_entry( "J" ) )
   {
      double expected_J = exp->double_data( "J" ) ;
      eq = PEL::double_equality( st->detF(), expected_J, D_EPS, D_MIN ) ;
      if( !eq ) display_error( expected_J, st->detF() ) ;
      notify_one_test_result( "detF", eq ) ;
   }

   st->destroy() ;
}

//----------------------------------------------------------------------------
bool
AP_KinematicState_TEST:: ok_rCG_inv_rCG( doubleArray2D const& rcg,
                                         doubleArray2D const& inv_rcg,
                                         size_t dim ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState_TEST:: ok_rCG_inv_rCG" ) ;
   bool eq = false ;

   bool result = true ;
   for( size_t i=0 ; i<dim ; ++i )
   {
      for( size_t j=0 ; j<dim ; ++j )
      {
         double xx = 0.0 ;
         for( size_t k=0 ; k<dim ; ++k )
         {
            xx += rcg( i, k ) * inv_rcg( k, j ) ;
         }
         double expected = ( i==j ? 1.0 : 0.0 ) ;
         eq = PEL::double_equality( xx, expected, D_EPS, D_MIN ) ;
         if( !eq ) display_error( xx, expected ) ;
         result = result && eq ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
AP_KinematicState_TEST:: ok_rCG_GL( doubleArray2D const& rcg,
                                    doubleArray2D const& gl,
                                    size_t dim ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState_TEST:: ok_rCG_GL" ) ;
   bool eq = false ;

   bool result = true ;
   for( size_t i=0 ; i<dim ; ++i )
   {
      for( size_t j=0 ; j<dim ; ++j )
      {
         double xx1 = rcg( i, j ) ;
         double xx2 = ( i==j ? 1.0 : 0.0 ) ;
         xx2 += 2.0 *  gl( i, j ) ;
         eq = PEL::double_equality( xx1, xx2, D_EPS, D_MIN ) ;
         if( !eq ) display_error( xx1, xx2 ) ;
         result = result && eq ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
AP_KinematicState_TEST:: ok_F_inv_F( doubleArray2D const& F,
                                     doubleArray2D const& inv_F,
                                     size_t dim ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState_TEST:: ok_F_inv_F" ) ;
   bool eq = false ;

   bool result = true ;
   for( size_t i=0 ; i<dim ; ++i )
   {
      for( size_t j=0 ; j<dim ; ++j )
      {
         double xx = 0.0 ;
         for( size_t k=0 ; k<dim ; ++k )
         {
            xx += inv_F( i, k ) * F( k, j ) ;
         }
         double expected = ( i==j ? 1.0 : 0.0 ) ;
         eq = PEL::double_equality( xx, expected, D_EPS, D_MIN ) ;
         if( !eq ) display_error( xx, expected ) ;
         result = result && eq ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------------
void
AP_KinematicState_TEST:: display_error( double xx_1, double xx_2 ) const
//----------------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = cout.flags() ;
   cout.setf( ios_base::uppercase | ios_base::scientific ) ;

   cout << std::setprecision( 10 ) << std::setw( 20 ) << xx_1
        << std::setprecision( 10 ) << std::setw( 20 ) << xx_2
        << endl ;

   cout.flags( original_flags ) ;
}
