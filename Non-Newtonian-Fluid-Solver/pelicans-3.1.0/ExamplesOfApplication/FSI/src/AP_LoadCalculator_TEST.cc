#include <AP_LoadCalculator_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <LA_SeqVector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_SetOfDomains.hh>

#include <FE.hh>
#include <FE_OneStepIteration.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>
#include <FE_SurfaceForce.hh>

#include <AP_LoadCalculator.hh>

#include <ios>
#include <iomanip>
#include <iostream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::string ;

AP_LoadCalculator_TEST const* 
AP_LoadCalculator_TEST::PROTOTYPE = new AP_LoadCalculator_TEST() ;

//---------------------------------------------------------------------------
AP_LoadCalculator_TEST:: AP_LoadCalculator_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "AP_LoadCalculator", "AP_LoadCalculator_TEST" )
{
}

//---------------------------------------------------------------------------
AP_LoadCalculator_TEST:: ~AP_LoadCalculator_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
AP_LoadCalculator_TEST:: process_one_test( PEL_ModuleExplorer const* exp ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator_TEST:: process_one_test" ) ;

   // -- similar to instantiations in FE_StepByStepProgression
   FE::set_geometry( FE::cartesian ) ;

   PEL_ModuleExplorer* ee = 
                       exp->create_subexplorer( 0, "PDE_SetOfDomains" ) ;
   PDE_SetOfDomains* sdoms = PDE_SetOfDomains::create( 0, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "FE_SetOfParameters" ) ;
   FE_SetOfParameters* prms = FE_SetOfParameters::create( sdoms, sdoms, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "FE_OneStepIteration" ) ;
   FE_OneStepIteration* one_it = 
                        FE_OneStepIteration::make( sdoms, sdoms, prms, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "FE_TimeIterator" ) ;
   FE_TimeIterator* t_it = FE_TimeIterator::create( 0, ee ) ;
   ee->destroy() ; ee=0 ;

   // -- perform one iteration... similar to FE_StepByStepProgression
   one_it->do_before_time_stepping( t_it ) ;
   t_it->start() ;
   one_it->do_before_inner_iterations_stage( t_it ) ;
   one_it->do_inner_iterations_stage( t_it ) ;
   one_it->do_after_inner_iterations_stage( t_it ) ;

   ee = exp->create_subexplorer( 0, "load_test" ) ;
   AP_LoadCalculator* lc = 
          AP_LoadCalculator::object( ee->string_data( "load_calculator" ) ) ;
   FE_SurfaceForce* sf = 
          FE_SurfaceForce::object( ee->string_data( "surface_force" ) ) ;
   double d_eps = ee->double_data( "dbl_epsilon" ) ;
   double d_min = ee->double_data( "dbl_minimum" ) ;
   ee->destroy() ; ee=0 ;

   bool eq = false ;

   sf->compute_force( t_it ) ;

   PDE_DiscreteField const* fv = lc->fluid_velocity() ;
   bool ok_0 = true ;
   doubleVector fl( fv->nb_components() ) ;
   for( size_t ic=0 ; ic<fv->nb_components() ; ++ic )
   {
      fl( ic ) = 0.0 ;
      for( size_t n=0 ; n<fv->nb_nodes() ; ++n )
      {
         fl( ic ) += lc->fluid_load( n, ic ) ;
      }
      eq = PEL::double_equality( fl( ic ), sf->force( ic ), d_eps, d_min ) ;
      if( !eq ) display_error( fl( ic ), sf->force( ic ) ) ;
      ok_0 = ok_0 && eq ;
   }

   notify_one_test_result( "fluid_load / surface force", ok_0 ) ;

   PDE_DiscreteField const* sd = lc->structure_displacement() ;
   PDE_LinkDOF2Unknown* sd_link = 
        PDE_LinkDOF2Unknown::create( 0, sd, 
                                     "sequence_of_the_components", true ) ;
   LA_SeqVector* load = LA_SeqVector::create( 0, sd_link->unknown_vector_size() ) ;
   lc->update_load( sd_link, load ) ;

   bool ok_1 = true ;
   bool ok_2 = true ;
   for( size_t ic=0 ; ic<sd->nb_components() ; ++ic )
   {
      double fs = 0.0 ;
      for( size_t n=0 ; n<sd->nb_nodes() ; ++n )
      {
         if( sd_link->DOF_is_unknown( n, ic ) )
         {
            size_t idx = sd_link->unknown_linked_to_DOF( n, ic ) ;
            ok_1 = ok_1 && ( lc->solid_load( n, ic ) == load->item( idx ) ) ;
         }
         fs += lc->solid_load( n, ic ) ;
      }
      eq = PEL::double_equality( fs, fl(ic), d_eps, d_min ) ;
      if( !eq ) display_error( fs, fl(ic) ) ;
      ok_2 = ok_2 && eq ;
   }

   notify_one_test_result( "update_load / solid_load", ok_1 ) ;
   notify_one_test_result( "fluid_load / solid_load", ok_2 ) ;

   load->destroy() ;
   sd_link->destroy() ;
   sdoms->destroy() ;
   t_it->destroy() ;
}

//----------------------------------------------------------------------------
void
AP_LoadCalculator_TEST:: display_error( double xx_1, double xx_2 ) const
//----------------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = cout.flags() ;
   cout.setf( ios_base::uppercase | ios_base::scientific ) ;

   cout << std::setprecision( 10 ) << std::setw( 20 ) << xx_1
        << std::setprecision( 10 ) << std::setw( 20 ) << xx_2
        << endl ;

   cout.flags( original_flags ) ;
}
