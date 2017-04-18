#include <AP_1Dcheck.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <FE_TimeIterator.hh>

#include <fstream>
#include <iostream>

using std::cout ;
using std::endl ;

AP_1Dcheck const* AP_1Dcheck::PROTOTYPE = new AP_1Dcheck() ;

//----------------------------------------------------------------------------
AP_1Dcheck:: AP_1Dcheck( void )
//----------------------------------------------------------------------------
   : PEL_Application( "AP_1Dcheck" )
{
}

//----------------------------------------------------------------------------
AP_1Dcheck*
AP_1Dcheck:: create_replica( PEL_Object* a_owner, 
                             PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_1Dcheck:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   AP_1Dcheck* result = new AP_1Dcheck( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
AP_1Dcheck:: AP_1Dcheck( PEL_Object* a_owner, 
                         PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , GEOMETRY( exp->string_data( "geometry" ) )
   , a( exp->double_data( "a" ) )
   , c0( exp->double_data( "c_init" ) )
   , b( exp->double_data( "b" ) )
   , RHO( exp->double_data( "fluid_density" ) )
   , VISC( exp->has_entry( "fluid_viscosity" ) ? 
           exp->double_data( "fluid_viscosity" ) : PEL::bad_double() ) 
   , YOUNG( exp->double_data( "Young_modulus" ) )
   , POIS( exp->has_entry( "Poisson_coefficient" ) ? 
           exp->double_data( "Poisson_coefficient" ) : PEL::bad_double() ) 
   , tauN( exp->double_data( "imposed_pressure" ) )
   , OMEGA( exp->has_entry( "relaxation_coefficient" ) ?
            exp->double_data( "relaxation_coefficient" ) : PEL::bad_double() )
   , TOLERANCE( exp->double_data( "tolerance" ) )
   , MY_EPS( exp->double_data( "dbl_epsilon" ) )
   , MY_MIN( exp->double_data( "dbl_minimum" ) )
   , TIME_ORDER( exp->int_data( "time_order" ) )
   , OUTPUT_FILE( exp->string_data( "output_file_name" ) )
{
   PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "FE_TimeIterator" ) ;
   TIME_IT = FE_TimeIterator::create( this, ee ) ;
   ee->destroy() ; ee = 0 ;
}

//----------------------------------------------------------------------------
AP_1Dcheck:: ~AP_1Dcheck( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
AP_1Dcheck:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_1Dcheck:: run" ) ;
   
   std::ofstream os( OUTPUT_FILE.c_str() ) ;
   if( !os ) PEL_Error::object()->raise_file_handling( OUTPUT_FILE, "open" ) ;

   if( GEOMETRY == "cartesian" )
   {
      run_cartesian( os ) ;
   }
   else if( GEOMETRY == "axisymetrical" )
   {
      run_axisymetrical( os ) ;
   }
   os.close() ;
}

//----------------------------------------------------------------------------
void
AP_1Dcheck:: run_cartesian( std::ostream& os )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_1Dcheck:: run_cartesian" ) ;

   double A = 0.0 ;
   double A_old = 0.0 ;
   double D = 0.0 ;
   double D_exp = 0.0 ;
   double D_exp_exp = 0.0 ;
   double D_point ;
   double c = c0 ;
   double c_exp = c0 ;
   double c_exp_exp = c0 ;
   double c_point ;
   double relax = PEL::bad_double() ;

   bool first_it = true ;
   for( TIME_IT->start() ; !TIME_IT->is_finished() ; TIME_IT->go_next_time() )
   {
      PEL::out() << endl << "++++++ ITERATION = " 
                 << TIME_IT->iteration_number() 
                 << " *** TIME = "      << TIME_IT->time() 
                 << " *** TIME STEP = " << TIME_IT->time_step() 
                 << " ++++++" << endl << endl ;
      double dt = TIME_IT->time_step() ;

      double ddd = PEL::max_double() ;
      bool stop = false ;
      size_t iter = 0 ;
      do
      {
         ++iter ;
         c = c0 + A*(1.-c0/b ) ;

         if( first_it || TIME_ORDER==1 )
         {
            c_point = (c-c_exp)/dt ;
         }
         else
         {
            c_point = (3.0*c - 4.0*c_exp + c_exp_exp)/2.0/dt ;
         }

         D = c_point ;

         A_old = A ;

         if( first_it || TIME_ORDER==1 )
         {
            D_point = (D-D_exp)/dt ;
         }
         else
         {
            D_point = (3.0*D - 4.0*D_exp + D_exp_exp)/2.0/dt ;
         }

         A = RHO * D_point * (c-a) + tauN ;
         A *= -b/YOUNG ;

         stop = PEL::double_equality( c-c0, A*(1.-c0/b), MY_EPS, MY_MIN ) ;
         ddd = PEL::abs( c - c0 - A*(1.-c0/b)) ;
//         stop = ( ddd < TOLERANCE ) ;

         if( !stop )
         {
            relax = relaxation_coefficient( A_old, A ) ;
            A = relax * A + (1.-relax) * A_old ;
         }
         PEL::out() << iter 
                    << "     " << ddd 
                    << "   c=" << c 
                    << "   D=" << D ;
         if( !stop ) PEL::out() << "  relax=" << relax ;
         PEL::out() << std::endl ;
      } while( !stop ) ;

      first_it = false ;
      c_exp_exp = c_exp ;
      c_exp = c ;
      D_exp_exp = D_exp ;
      D_exp = D ;

      os << TIME_IT->time() << "  " << A*(1.-c0/b ) << endl ;
   }
}

//----------------------------------------------------------------------------
void
AP_1Dcheck:: run_axisymetrical( std::ostream& os )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_1Dcheck:: run_axisymetrical" ) ;

   double A = 0.0 ;
   double A_old = 0.0 ;
   double D = 0.0 ;
   double D_exp = 0.0 ;
   double D_exp_exp = 0.0 ;
   double D_point ;
   double c = c0 ;
   double c_exp = c0 ;
   double c_exp_exp = c0 ;
   double c_point ;
   double relax = PEL::bad_double() ;

   bool first_it = true ;
   for( TIME_IT->start() ; !TIME_IT->is_finished() ; TIME_IT->go_next_time() )
   {
      PEL::out() << endl << "++++++ ITERATION = " 
                 << TIME_IT->iteration_number() 
                 << " *** TIME = "      << TIME_IT->time() 
                 << " *** TIME STEP = " << TIME_IT->time_step() 
                 << " ++++++" << endl << endl ;
      double dt = TIME_IT->time_step() ;

      double ddd = PEL::max_double() ;
      bool stop = false ;
      size_t iter = 0 ;
      do
      {
         ++iter ;
         c = c0 + A*(c0-b*b/c0 ) ;

         if( first_it || TIME_ORDER==1 )
         {
            c_point = (c-c_exp)/dt ;
         }
         else
         {
            c_point = (3.0*c - 4.0*c_exp + c_exp_exp)/2.0/dt ;
         }

         D = c * c_point ;

         A_old = A ;

         if( first_it || TIME_ORDER==1 )
         {
            D_point = (D-D_exp)/dt ;
         }
         else
         {
            D_point = (3.0*D - 4.0*D_exp + D_exp_exp)/2.0/dt ;
         }

         double pp = - RHO * ( D_point*(PEL::log(c)-PEL::log(a)) + 
                               D*D/2.*(1./c/c-1./a/a) )
                     - tauN - 2.0*VISC*D/a/a ;
         double coef = YOUNG/(1.0+POIS) * (1.+b*b/c0/c0 )
                     + 2.0*YOUNG*POIS/(1.+POIS)/(1.-2.*POIS) ;
         A = ( - pp - 2.0*VISC*D/c/c ) / coef  ;

         stop = PEL::double_equality( c-c0, A*(c0-b*b/c0), MY_EPS, MY_MIN ) ;
         ddd = PEL::abs( c - c0 - A*(c0-b*b/c0)) ;
//         stop = ( ddd < TOLERANCE ) ;

         if( !stop )
         {
            relax = relaxation_coefficient( A_old, A ) ;
            A = relax * A + (1.-relax) * A_old ;
         }
         PEL::out() << iter 
                    << "     " << ddd 
                    << "   c=" << c 
                    << "   c_point=" << c_point ;
         if( !stop ) PEL::out() << "  relax=" << relax ;
         PEL::out() << std::endl ;
      } while( !stop ) ;

      first_it = false ;
      c_exp_exp = c_exp ;
      c_exp = c ;
      D_exp_exp = D_exp ;
      D_exp = D ;

      os << TIME_IT->time() << "  " << A*(c0-b*b/c0 ) << endl ;
   }
}

//--------------------------------------------------------------------------
double
AP_1Dcheck:: relaxation_coefficient( double A_old, double A ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "AP_1Dcheck:: relaxation_coefficient" ) ;

   static double result = PEL::bad_double() ;

   static double DeltaA     = PEL::bad_double() ;
   static double DeltaA_old = PEL::bad_double() ;

   static double mu = PEL::bad_double() ;

   if( OMEGA == PEL::bad_double() )
   {
      if( result == PEL::bad_double() )
      {
         DeltaA = A_old - A ;
         mu = 0.0 ;
         result = 1.0 ;
      }
      else
      {
         DeltaA_old = DeltaA ;
         DeltaA = A_old - A;
//                   cout << "DeltaA_old=" << DeltaA_old << endl ;
//                   cout << "DeltaA=" << DeltaA << endl ;
         mu = mu + (mu-1.)*DeltaA/(DeltaA_old-DeltaA) ;
         result = 1. - mu ;
      }
   }
   else
   {
      result = OMEGA ;
   }
   return( result ) ;
}
