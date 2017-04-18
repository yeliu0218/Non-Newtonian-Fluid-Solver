#include <AP_FluidStructureMaster.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_System.hh>
#include <PEL_Root.hh>
#include <PEL_Timer.hh>
#include <PEL_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_SetOfPoints.hh>

#include <LA_SeqVector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_InterfaceAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfDomains.hh>

#include <FE.hh>
#include <FE_OneStepIteration.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <intVector.hh>
#include <doubleVector.hh>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::ostringstream ;
using std::string ;

size_t AP_FluidStructureMaster::graphics_level = 0 ;

AP_FluidStructureMaster const* 
AP_FluidStructureMaster:: PROTOTYPE = new AP_FluidStructureMaster() ;

//----------------------------------------------------------------------
AP_FluidStructureMaster:: AP_FluidStructureMaster( void )
//----------------------------------------------------------------------
   : PEL_Application( "AP_FluidStructureMaster" )
   , TIME_IT( 0 )
   , graphics_times( 0 )
   , graphics_next_time( 0.0 )
   , saver_times( 0 )
   , saver_next_time( 0 )
   , overall( 0 )
   , SAVEFG( true )
   , VERBOSE( true )
{
}

//----------------------------------------------------------------------
AP_FluidStructureMaster*
AP_FluidStructureMaster:: create( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   AP_FluidStructureMaster* result = 
                              new AP_FluidStructureMaster( a_owner, exp ) ;
   result->register_storable_objects() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
AP_FluidStructureMaster*
AP_FluidStructureMaster:: create_replica( 
                                      PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   AP_FluidStructureMaster* result = 
                              new AP_FluidStructureMaster( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
AP_FluidStructureMaster:: AP_FluidStructureMaster( 
                                            PEL_Object* a_owner,
                                            PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , TIME_IT( 0 )
   , SDOMS( 0 )
   , F_DOM( 0 )
   , S_DOM( 0 )
   , PRMS( 0 )
   , FGRID( 0 )
   , FLUID( 0 )
   , LOADC( 0 )
   , SOLID( 0 )
   , FGD( 0 )
   , FV( 0 )
   , FP( 0 )
   , SD( 0 )
   , SV( 0 )
   , L_NEW( PEL::bad_index() )
   , L_EXPLICIT( PEL::bad_index() )
   , L_EXPLICIT_EXPLICIT( PEL::bad_index() )
   , L_OLD( PEL::bad_index() )
   , VERTS_INI( 0 )
   , VERTS_NEW( 0 )
   , bFE( 0 )
   , cFE( 0 )
   , INTERF_COL( 0 ) 
   , TT( Weak )
   , N_d_INTERF( PEL::bad_index() )
   , DELTA_d_INTERF( 0 )
   , DELTA_d_INTERF_OLD( 0 )
   , OMEGA( PEL::bad_double() )
   , TOL_EPS( PEL::bad_double() )
   , TOL_MIN( PEL::bad_double() )
   , SAVE_ALL( false )
   , graphics_times( 0 )
   , graphics_next_time( PEL::max_double() )
   , saver_times( 0 )
   , saver_next_time( PEL::max_double() )
   , overall( PEL_Timer::create( this ) )
   , SAVEFG( true )
   , DBL_EPS_GRID_POS( 1.E-4 )
   , DBL_MIN_GRID_POS( 1.E-8 )
   , VERBOSE( true )
{
   PEL_LABEL( "AP_FluidStructureMaster:: AP_FluidStructureMaster" ) ;

   if( exp->has_entry( "graphics_output_times" ) )
   {
      graphics_times = exp->doubleVector_data( "graphics_output_times" ) ;
   }

   if( exp->has_entry( "state_saving_times" ) )
   {
      saver_times = exp->doubleVector_data( "state_saving_times" ) ;
   }
   if( exp->has_entry( "save_grid_and_fields_for_postprocessing" ) )
   {
      SAVEFG = exp->bool_data( "save_grid_and_fields_for_postprocessing" ) ;
   }

   if( exp->has_entry( "verbose" ) ) VERBOSE = exp->bool_data( "verbose" ) ;

   FE::set_geometry( FE::cartesian ) ;

   PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "FE_TimeIterator" ) ;
   TIME_IT = FE_TimeIterator::create( this, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "PDE_SetOfDomains" ) ;
   SDOMS = PDE_SetOfDomains::create( this, ee ) ;
   ee->destroy() ;

   F_DOM = SDOMS->domain( "fluid" ) ;
   FGD = F_DOM->set_of_discrete_fields()->item( "fluid_grid_displacement" ) ;
   FV = F_DOM->set_of_discrete_fields()->item( "fluid_velocity" ) ;
   FP = F_DOM->set_of_discrete_fields()->item( "fluid_pressure" ) ;

   S_DOM = SDOMS->domain( "structure" ) ;
   SD = S_DOM->set_of_discrete_fields()->item( "structure_displacement" ) ;
   SV = S_DOM->set_of_discrete_fields()->item( "structure_velocity" ) ;

   bFE = F_DOM->create_LocalFEbound( this ) ;
   bFE->require_field_calculation( FV, PDE_LocalFE::node ) ;
   bFE->require_field_calculation( FGD, PDE_LocalFE::node ) ;
   bFE->require_field_calculation( FGD, PDE_LocalFE::N ) ; //??? pour les tests
   bFE->require_field_calculation( SD, PDE_LocalFE::N ) ;

   cFE = F_DOM->create_LocalFEcell( this ) ;
   cFE->require_field_calculation( FGD, PDE_LocalFE::N ) ;

   //???? tout construire avec les bons domaines
   ee = exp->create_subexplorer( 0, "FE_SetOfParameters" ) ;
   PRMS = FE_SetOfParameters::create( this, SDOMS, ee ) ;
   ee->destroy() ;
   
   ee = exp->create_subexplorer( 0, "fluid_grid_solver" ) ;
   FGRID = FE_OneStepIteration::make( this, SDOMS, PRMS, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "fluid_solver" ) ;
   FLUID = FE_OneStepIteration::make( this, SDOMS, PRMS, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "load_calculator" ) ;
   LOADC = FE_OneStepIteration::make( this, SDOMS, PRMS, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "structure_solver" ) ;
   SOLID = FE_OneStepIteration::make( this, SDOMS, PRMS, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "fluid_structure_coupling" ) ;
   INTERF_COL = GE_Color::object( 
      ee->string_data( "color_of_fluid_interface_with_structure" ) ) ;
   L_NEW = ee->int_data( "level_of_new" ) ;
   L_EXPLICIT = ee->int_data( "level_of_explicit" ) ;
   L_EXPLICIT_EXPLICIT = ee->int_data( "level_of_explicit_explicit" ) ;
   L_OLD = ee->int_data( "level_of_old_in_internal_iteration" ) ;
   BC_FV_INTERF = 
         ee->int_data( "boundary_condition_of_fluid_velocity_at_interface" ) ;
   if( BC_FV_INTERF!=1 && BC_FV_INTERF!=2 )
      PEL_Error::object()->raise_bad_data_value( ee, 
         "boundary_condition_of_fluid_velocity_at_interface", "   1\n   2" ) ;

   PEL_ModuleExplorer const* eee = ee->create_subexplorer( 0, "strategy" ) ;
   std::string const& type = eee->string_data( "type" ) ;
   if( type == "fixed_point" )
   {
      TT = FixedPoint ;
      TOL_EPS = eee->double_data( "tol_epsilon" ) ;
      TOL_MIN = eee->double_data( "tol_minimum" ) ;
      OMEGA = eee->double_data( "relaxation_coefficient" ) ;
      if( eee->has_entry( "save_all_iterations" ) )
      {
         SAVE_ALL = eee->bool_data( "save_all_iterations" ) ;
      }
   }
   else if( type == "fixed_point_Aitken" )
   {
      TT = FixedPointAitken ;
      TOL_EPS = eee->double_data( "tol_epsilon" ) ;
      TOL_MIN = eee->double_data( "tol_minimum" ) ;
      if( eee->has_entry( "save_all_iterations" ) )
      {
         SAVE_ALL = eee->bool_data( "save_all_iterations" ) ;
      }
   }
   else if( type == "weak" )
   {
      TT = Weak ;
   }
   else
      PEL_Error::object()->raise_bad_data_value( eee, "type", 
                 "  \"fixed_point\"\n  \"fixed_point_Aitken\"\n  \"weak\"" ) ;
   eee->destroy() ;
   if( ee->has_module( "grid_position_check" ) )
   {
      eee = ee->create_subexplorer( 0, "grid_position_check" ) ;
      DBL_EPS_GRID_POS = eee->double_data( "dbl_epsilon" ) ;
      DBL_MIN_GRID_POS = eee->double_data( "dbl_minimum" ) ;
      eee->destroy() ;
   }
   ee->destroy() ;

   size_t nb_dims = F_DOM->nb_space_dimensions() ;
   check_field_nb_components( FGD, nb_dims ) ;
   check_field_nb_components( FV, nb_dims ) ;
   check_field_nb_components( FP, 1 ) ;
   check_field_nb_components( SD, nb_dims ) ;
   check_field_nb_components( SV, nb_dims ) ;

   check_field_storage_depth( FGD, L_NEW ) ;
   check_field_storage_depth( FGD, L_EXPLICIT ) ;
   check_field_storage_depth( FV, L_NEW ) ;
   check_field_storage_depth( FV, L_EXPLICIT ) ;
   check_field_storage_depth( FV, L_EXPLICIT_EXPLICIT ) ;
   check_field_storage_depth( FGD, L_EXPLICIT_EXPLICIT ) ;
   check_field_storage_depth( FP, L_NEW ) ;
   check_field_storage_depth( FP, L_EXPLICIT ) ;
   check_field_storage_depth( SD, L_NEW ) ;
   check_field_storage_depth( SD, L_EXPLICIT ) ;
   check_field_storage_depth( SV, L_NEW ) ;
   check_field_storage_depth( SV, L_EXPLICIT ) ;
   check_field_storage_depth( SD, L_OLD ) ;
   check_field_storage_depth( SV, L_OLD ) ;

   GE_SetOfPoints const* vvvs = F_DOM->set_of_vertices() ;
   VERTS_INI = PEL_Vector::create( this, 
                                   F_DOM->set_of_vertices()->nb_points() ) ;
   VERTS_NEW = PEL_Vector::create( this, 
                                   F_DOM->set_of_vertices()->nb_points() ) ;
   for( size_t iv=0 ; iv<vvvs->nb_points() ; ++iv )
   {
      VERTS_INI->set_at( iv, vvvs->point( iv )->create_clone( VERTS_INI ) ) ;
      VERTS_NEW->set_at( iv, GE_Point::create( VERTS_NEW, nb_dims ) ) ;
   }

   FE_OneStepIteration::increase_indent() ;
   if( VERBOSE && PEL_Exec::communicator()->rank() == 0 )
   {
      PEL::out() << "   fluid_grid_solver" << endl ;
      FGRID->print( PEL::out(), 4 ) ;
      PEL::out() << "   fluid_solver" << endl ;
      FLUID->print( PEL::out(), 4 ) ;
      PEL::out() << "   load_calculator" << endl ;
      LOADC->print( PEL::out(), 4 ) ;
      PEL::out() << "   structure_solver" << endl ;
      SOLID->print( PEL::out(), 4 ) ;
   }

   overall->start() ;
}

//----------------------------------------------------------------------
AP_FluidStructureMaster:: ~AP_FluidStructureMaster( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
FE_TimeIterator const* 
AP_FluidStructureMaster:: time_iterator( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: time_iterator" ) ;

   FE_TimeIterator const* result = TIME_IT ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_SetOfParameters const* 
AP_FluidStructureMaster:: set_of_parameters( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: set_of_parameters" ) ;

   FE_SetOfParameters const* result = PRMS ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
AP_FluidStructureMaster:: add_storable_objects( PEL_ListIdentity* list )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: add_storable_objects" ) ;
   list->extend( TIME_IT ) ;

   for( size_t i=0 ; i<SDOMS->nb_domains() ; ++i )
   {
      list->extend( SDOMS->domain( i ) ) ;
   }
   for( size_t i=0 ; i<SDOMS->nb_interfaces() ; ++i )
   {
      list->extend( SDOMS->interface( i ) ) ;
   }
   FLUID->register_storable_objects( list ) ;
   SOLID->register_storable_objects( list ) ;
   LOADC->register_storable_objects( list ) ;
}

//----------------------------------------------------------------------
void
AP_FluidStructureMaster:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: run" ) ; 

   double xx ;

   std::ofstream osfp( "fixed_point.txt" ) ;
   if( !osfp ) 
      PEL_Error::object()->raise_file_handling( "fixed_point.txt", "open" ) ;
   osfp << "#" << std::setw( 9 ) << "time" 
        << std::setw( 20 ) << "nb_fixed_point_its" << std::endl ;

   check_times( "graphics_output_times", graphics_times ) ;
   check_times( "state_saving_times", saver_times ) ;
   
   graphics_next_time = next_time( graphics_times ) ;
   saver_next_time = next_time( saver_times ) ;
   
   PEL::out() << endl << "++++++ TIME STEPPING " 
        << " *** INITIAL TIME = " << TIME_IT->initial_time() 
        << " *** FINAL TIME = "   << TIME_IT->final_time() 
        << " ++++++" << endl ;

   FGRID->do_before_time_stepping( TIME_IT ) ;
   FLUID->do_before_time_stepping( TIME_IT ) ;
   LOADC->do_before_time_stepping( TIME_IT ) ;
   SOLID->do_before_time_stepping( TIME_IT ) ;

   save_for_post_processing() ;
   if( VERBOSE ) PEL::out() << endl ;

   double time = 0.0 ;
   double dt = 0.0 ;
   for( TIME_IT->start() ; !TIME_IT->is_finished() ; TIME_IT->go_next_time() )
   {
      time = TIME_IT->time() ;
      dt = TIME_IT->time_step() ;

      if( VERBOSE ) PEL::out() << endl << "++++++ ITERATION = " 
                               << TIME_IT->iteration_number() 
                               << " *** TIME = "      << time
                               << " *** TIME STEP = " << dt 
                               << " ++++++" << endl << endl ;

      if( VERBOSE) PEL::out() << "   do_before_inner_iterations_stage..." 
                              << endl ;
      FGRID->do_before_inner_iterations_stage( TIME_IT ) ;
      FLUID->do_before_inner_iterations_stage( TIME_IT ) ;
      LOADC->do_before_inner_iterations_stage( TIME_IT ) ;
      SOLID->do_before_inner_iterations_stage( TIME_IT ) ;

      size_t nb_iter = 0 ;
      bool stop = true ;
      double relax = PEL::bad_double() ;

      // extrapolation du déplacement de la structure
      if( VERBOSE ) PEL::out() << "   extrapolation of structure displacement"
                               << endl ;
      for( size_t n=0 ; n<SD->nb_nodes() ; ++n )
         for( size_t ic=0 ; ic<SD->nb_components() ; ++ic )
         {
            if( !SD->DOF_has_imposed_value( n, ic ) )
            {
               xx =      SD->DOF_value( L_EXPLICIT, n, ic ) +
                    dt * SV->DOF_value( L_EXPLICIT, n, ic ) ;
               SD->set_DOF_value( L_NEW, n, xx, ic ) ;
            }
         }
      
      do
      {
         nb_iter++ ;

         if( VERBOSE ) 
            PEL::out() << "   " << nb_iter 
                       << ".1 update of grid displacement at "
                       << "interface \"" << INTERF_COL->name()
                       << "\"" << endl ;
         update_FGD_at_interf() ;

         // calcul du déplacement de la grille
         if( VERBOSE ) 
            PEL::out() << "   " << nb_iter 
                       << ".2 computation of grid displacement increment"
                       << endl ;
         FGRID->do_inner_iterations_stage( TIME_IT ) ;

         // déplacement de la grille
         if( VERBOSE )
            PEL::out() << "   " << nb_iter
                       << ".3 grid move" << endl ;
         change_fluid_grid_position() ;

         if( VERBOSE ) 
            PEL::out() << "   " << nb_iter 
                       << ".4 update of fluid velocity at interface \"" 
                       << INTERF_COL->name() << "\""  << endl ;
         update_FV_at_interf( dt ) ;

         if( VERBOSE )
            PEL::out() << "   " << nb_iter 
                       << ".5 fluid solver" << endl ;
         FLUID->do_inner_iterations_stage( TIME_IT ) ;

         if( VERBOSE )
            PEL::out() << "   " << nb_iter 
                       << ".6 load calculation" << endl ;
         LOADC->do_inner_iterations_stage( TIME_IT ) ;

         SD->copy_DOFs_value( L_NEW, L_OLD ) ;
         SV->copy_DOFs_value( L_NEW, L_OLD ) ;

         if( VERBOSE )
            PEL::out() << "   " << nb_iter  
                       << ".7 structure solver" << endl ;
         SOLID->do_inner_iterations_stage( TIME_IT ) ;

         if( TT != Weak ) 
         {
            if( VERBOSE ) 
               PEL::out() << "   " << nb_iter  
                          << ".8 convergence check" << endl ;
            stop = check_convergence() ;

            if( VERBOSE ) 
               PEL::out() << "   " << nb_iter  
                          << ".9 compute relaxation coefficient" << endl ;
            relax = relaxation_coefficient() ;
         }

         if( !stop && ( relax != PEL::bad_double() ) )
         {
            if( VERBOSE )
               PEL::out() << "   " << nb_iter 
                          << ".10 structure relaxation" << endl ;
            for( size_t n=0 ; n<SD->nb_nodes() ; ++n )
               for( size_t ic=0 ; ic<SD->nb_components() ; ++ic )
               {
                  xx =       relax * SD->DOF_value( L_NEW, n, ic ) +
                       (1.0-relax) * SD->DOF_value( L_OLD, n, ic ) ;
                  SD->set_DOF_value( L_NEW, n, xx, ic ) ;

                  xx =       relax * SV->DOF_value( L_NEW, n, ic ) +
                       (1.0-relax) * SV->DOF_value( L_OLD, n, ic ) ;
                  SV->set_DOF_value( L_NEW, n, xx, ic ) ;
               }
               PEL::out() << "   - - - - - - - - - - " << endl ;
         }

         if( SAVE_ALL ) save_for_post_processing() ;

      } while( !stop ) ;

      FGRID->do_after_inner_iterations_stage( TIME_IT ) ;
      FLUID->do_after_inner_iterations_stage( TIME_IT ) ;
      LOADC->do_after_inner_iterations_stage( TIME_IT ) ;
      SOLID->do_after_inner_iterations_stage( TIME_IT ) ;

      FGD->copy_DOFs_value( L_EXPLICIT, L_EXPLICIT_EXPLICIT ) ;
      FGD->copy_DOFs_value( L_NEW, L_EXPLICIT ) ;

      FV->copy_DOFs_value( L_EXPLICIT, L_EXPLICIT_EXPLICIT ) ;
      FV->copy_DOFs_value( L_NEW, L_EXPLICIT ) ;

      FP->copy_DOFs_value( L_NEW, L_EXPLICIT ) ;

      SD->copy_DOFs_value( L_NEW, L_EXPLICIT ) ;
      SV->copy_DOFs_value( L_NEW, L_EXPLICIT ) ;

      if( !SAVE_ALL ) save_for_post_processing() ;

      save_for_restart() ;

      if( VERBOSE ) print_memory_info() ;
      osfp << std::setw( 10 ) << TIME_IT->time() 
           << std::setw( 20 ) << nb_iter << std::endl ;
   }

   if( VERBOSE )
      PEL::out() << endl << "++++++ TIME STEPPING COMPLETED " 
           << " ++++++" << endl << endl ;

   FGRID->do_after_time_stepping() ;
   FLUID->do_after_time_stepping() ;
   LOADC->do_after_time_stepping() ;
   SOLID->do_after_time_stepping() ;

   save_for_post_processing() ;
   
   overall->stop() ;
   
//    if( !VERBOSE ) 
//    {
//       PEL::out() << endl << "++++++ " << TIME_IT->iteration_number() 
//            << " ITERATIONS ++++++" << endl << endl ;
//    }
   PEL::out() << endl << "++++++ TIMERS  ++++++" << endl ;

   PEL::out() << endl << "total time : " ;
   overall->print( PEL::out(), 0 ) ;
   PEL::out() << endl << endl ;
   FE_OneStepIteration::print_standard_times( PEL::out(), 1 ) ;
   FGRID->print_additional_times( PEL::out(), 1 ) ;
   FLUID->print_additional_times( PEL::out(), 1 ) ;
   LOADC->print_additional_times( PEL::out(), 1 ) ;
   SOLID->print_additional_times( PEL::out(), 1 ) ;
   PEL::out() << endl ;
   
   osfp.close() ;
}

//----------------------------------------------------------------------------
void
AP_FluidStructureMaster:: update_FGD_at_interf( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: update_FGD_at_interf" ) ;

   size_t n = 0 ;

   boolVector done( FGD->nb_nodes() ) ;

   size_t nb_DOFs_updated = 0 ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   { 
      GE_Color const* color = bFE->color() ;
      if( color->is_matching( INTERF_COL ) )
      {
         for( size_t loc_n=0 ; loc_n<bFE->nb_local_nodes( FGD ) ; ++loc_n )
         {
            if( bFE->local_node_is_in_mesh( FGD, loc_n ) )
	    {
	       n = bFE->global_node( FGD, loc_n ) ;
	       if( ! done( n ) )
	       {
	          done( n ) = true ;
 	          for( size_t ic=0 ; ic<FGD->nb_components() ; ++ic )
		  {
                     if( !FGD->DOF_has_imposed_value( n, ic ) )
                     {
                        PEL_Error::object()->raise_plain( "error 0" ) ;
                     }
                     else
                     {
                        // on suppose les éléments conformes
                        // on suppose les memes noeuds pour DeltaD et FGD 
	                GE_Point const* pt = 
                                 bFE->local_node_location( FGD, loc_n ) ;
                        bFE->set_calculation_point( pt ) ;
                        PEL_ASSERT( bFE->nb_local_nodes( SD ) != 0 ) ;
                        double d_new = bFE->value_at_pt( SD, L_NEW, ic ) ;
                        FGD->set_DOF_imposed_value( n, d_new, ic ) ;
                        FGD->set_DOF_value( L_NEW, n, d_new, ic ) ;
                        nb_DOFs_updated++ ;
                     }
	          }
	       }
	    }
	 }
      }
   }
   PEL::out() << "      " << nb_DOFs_updated << " DOFs (imposed and level="
              << L_NEW << ") of \"" 
              << FGD->name() << "\" updated" << endl ;

   if( N_d_INTERF == PEL::bad_index() )
   {
      N_d_INTERF = nb_DOFs_updated ;
      PEL_ASSERT( DELTA_d_INTERF_OLD==0 && DELTA_d_INTERF==0 ) ;
      DELTA_d_INTERF = LA_SeqVector::create( this, N_d_INTERF ) ;
      DELTA_d_INTERF_OLD = LA_SeqVector::create( this, N_d_INTERF ) ;
   }
   else
   {
      PEL_ASSERT( nb_DOFs_updated == N_d_INTERF ) ;
   }
}

//----------------------------------------------------------------------------
void
AP_FluidStructureMaster:: update_FV_at_interf( double dt )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: update_FV_at_interf" ) ;

   size_t n = 0 ;

   boolVector done( FV->nb_nodes() ) ;

   size_t nb_DOFs_updated = 0 ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   { 
      GE_Color const* color = bFE->color() ;
      if( color->is_matching( INTERF_COL ) )
      {
         for( size_t loc_n=0 ; loc_n<bFE->nb_local_nodes( FV ) ; ++loc_n )
         {
	    if( bFE->local_node_is_in_mesh( FV, loc_n ) )
	    {
	       n = bFE->global_node( FV, loc_n ) ;
               GE_Point const* pt = bFE->local_node_location( FV, loc_n ) ;
	       if( ! done( n ) )
	       {
	          done( n ) = true ;
//                  cout << "n=" << n << endl ;
 	          for( size_t ic=0 ; ic<FV->nb_components() ; ++ic )
		  {
                     if( !FV->DOF_has_imposed_value( n, ic ) )
                     {
                        cout << "--------- bound----------" << endl ;
                        bFE->print( cout, 3 ) ;
                        cout << endl ;
                        cout << "--------- node-----------" << endl ;
                        cout << "   node : " << n << endl ;
                        pt->print( cout, 3 ) ;
                        cout << endl ;
                        PEL_Error::object()->raise_plain( "error update FV" ) ;
                     }
                     else
                     {
                        // on suppose les éléments conformes
                        bFE->set_calculation_point( pt ) ;
                        PEL_ASSERT( bFE->nb_local_nodes( SD ) != 0 ) ;
                        double d_new = 
                           bFE->value_at_pt( FGD, L_NEW, ic ) ;
                        double d_exp = 
                           bFE->value_at_pt( FGD, L_EXPLICIT, ic ) ;
                        double xx = PEL::bad_double() ;
                        if( BC_FV_INTERF == 2 )
                        {
                           double d_exp_exp = 
                             bFE->value_at_pt( FGD, L_EXPLICIT_EXPLICIT, ic ) ;
                           xx = ( 3.0*d_new - 4.0*d_exp + d_exp_exp )
                                / 2.0 / dt ;
                        }
                        else if( BC_FV_INTERF == 1 )
                        {
                           xx = ( d_new - d_exp ) / dt ;
                        }
//                        std::cout << "imposed velo=" << xx << std::endl ;
                        FV->set_DOF_imposed_value( n, xx, ic ) ;
                        FV->set_DOF_value( L_NEW, n, xx, ic ) ;
                        nb_DOFs_updated++ ;
                     }
	          }
	       }
	    }
	 }
      }
   }
   PEL::out() << "      " << nb_DOFs_updated << " DOFs (imposed and level="
              << L_NEW << ") of \"" 
              << FV->name() << "\" updated" << endl ;
}

//----------------------------------------------------------------------
void
AP_FluidStructureMaster:: change_fluid_grid_position( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: change_fluid_grid_position" ) ;
   std::string const banner =
            "*** AP_FluidStructureMaster : \""+FGD->name()+"\"" ;
   // Le choix des deux constantes suivantes est sujet à caution

   size_t nb_dims = F_DOM->nb_space_dimensions() ;
   GE_SetOfPoints const* vvvs = F_DOM->set_of_vertices() ;
   
   size_t const nbvs = vvvs->nb_points() ;
   for( size_t iv=0 ; iv<nbvs ; ++iv )
   {
      GE_Point* pt = static_cast<GE_Point*>( VERTS_NEW->at( iv ) ) ;
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         pt->set_coordinate( ic, PDE_ResultSaver::undefined_value() ) ;
      }
   }

   boolVector done( vvvs->nb_points() ) ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      GE_Mpolyhedron const* poly = cFE->polyhedron() ;

      for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
      {
         GE_Point const* pt = poly->vertex( i ) ;
         PEL_CHECK( vvvs->has( pt ) ) ;
         size_t iv = vvvs->index( pt ) ;
         GE_Point* pos_ini = static_cast<GE_Point*>( VERTS_INI->at( iv ) ) ;
         GE_Point* pos_new = static_cast<GE_Point*>( VERTS_NEW->at( iv ) ) ;

         cFE->set_calculation_point( pt ) ;
         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            double val = cFE->value_at_pt( FGD, L_NEW, ic ) +
                         pos_ini->coordinate( ic ) ;

            // the value at a vertex obtained from different meshes should
            // be the same
            PDE_ResultSaver::check_value_consistency_at_vertex(
               banner, pt, pos_new->coordinate( ic ), val, 
               DBL_EPS_GRID_POS, DBL_MIN_GRID_POS ) ;
//             std::cout << "iv=" << iv << " ic=" << ic << " val=" << val << std::endl ;
            pos_new->set_coordinate( ic, val ) ;
         }
         done( iv ) = true ;
      }
   }

   //???? couteux et probablement inutile ???
   for( size_t iv=0 ; iv<vvvs->nb_points() ; ++iv ) PEL_ASSERT( done( iv ) ) ;
   
   vvvs->take_position_and_update( VERTS_NEW ) ;
}

//----------------------------------------------------------------------
bool
AP_FluidStructureMaster:: check_convergence( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: check_convergence" ) ;

   DELTA_d_INTERF->synchronize() ;
   DELTA_d_INTERF_OLD->set( DELTA_d_INTERF ) ;

   size_t n = 0 ;

   double dd = 0.0 ;
   size_t nb_pts = 0 ;

   boolVector done( FGD->nb_nodes() ) ;

   bool result = true ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   { 
      GE_Color const* color = bFE->color() ;
      if( color->is_matching( INTERF_COL ) )
      {
         for( size_t loc_n=0 ; loc_n<bFE->nb_local_nodes( FGD ) ; ++loc_n )
         {
	    if( bFE->local_node_is_in_mesh( FGD, loc_n ) )
	    {
	       n = bFE->global_node( FGD, loc_n ) ;
	       if( ! done( n ) )
	       {
	          done( n ) = true ;
 	          for( size_t ic=0 ; ic<FGD->nb_components() ; ++ic )
		  {
  	             GE_Point const* pt = 
                                     bFE->local_node_location( FGD, loc_n ) ;
                     bFE->set_calculation_point( pt ) ;
                     PEL_ASSERT( bFE->nb_local_nodes( SD ) != 0 ) ;
                     double d_solid = bFE->value_at_pt( SD, L_NEW, ic ) ;
                     double d_grid = FGD->DOF_value( L_NEW, n, ic ) ;
                     
                     bool eq = check_closeness( d_solid, d_grid ) ;

//                      std::cout << "(" << eq << ") --- ic=" << ic 
//                                << " d_solid=" << d_solid
//                                << " d_grid=" << d_grid 
//                                << " diff=" << d_solid-d_grid << endl ;

                     result = result && eq ;

                     dd += (d_solid-d_grid)*(d_solid-d_grid) ;
//???? on parcours tjrs dans le meme ordre
                     DELTA_d_INTERF->set_item( nb_pts, d_grid-d_solid ) ;
                     nb_pts++ ;
	          }
	       }
	    }
	 }
      }
   }
   dd = PEL::sqrt( dd ) / nb_pts ;

   PEL_ASSERT( nb_pts == N_d_INTERF ) ;

   if( VERBOSE )
   {
      PEL::out() << "      2_norm( d_fluid - d_solid )_interf = " 
                 << dd << "   (" << nb_pts << ")" << endl ;
   }
//   bool result = ( dd < TOLER ? true : false ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
AP_FluidStructureMaster:: check_closeness( double d_solid, 
                                           double d_grid ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: check_closeness" ) ;

   bool result = false ;

   double abs_x = PEL::abs( d_solid ) ;
   double abs_y = PEL::abs( d_grid ) ;

   if( abs_y < TOL_MIN ) // test for closeness to 0.0
   {
      // |y| is undistinguishable from 0.0
      if( abs_x < TOL_MIN)
      {
         // |x| is undistinguishable from 0.0
         result = true ;
      }
      else
      {
	 result = false ;
      }
   }
   else if( abs_x > TOL_MIN ) // test for closeness to each other
   {
      if( PEL::abs( d_solid - d_grid ) < TOL_EPS )
      {
         result = true ;
      }
      else
      {
         result = false ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
AP_FluidStructureMaster:: relaxation_coefficient( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: relaxation_coefficient" ) ;

   PEL_ASSERT( TT==FixedPoint || TT==FixedPointAitken ) ;

   static double result = PEL::bad_double() ;

   static double mu = PEL::bad_double() ;

   if( TT == FixedPointAitken )
   {
      if( result == PEL::bad_double() )
      {
         mu = 0.0 ;
         result = 1.0 ;
      }
      else
      {
         DELTA_d_INTERF_OLD->sum( DELTA_d_INTERF, -1.0 ) ;
         
         mu = mu + (mu-1.)*DELTA_d_INTERF_OLD->dot( DELTA_d_INTERF )/
                           DELTA_d_INTERF_OLD->dot( DELTA_d_INTERF_OLD ) ;
         result = 1. - mu ;
      }
   }
   else if( TT == FixedPoint )
   {
      result = OMEGA ;
   }

   if( VERBOSE )
   {
      PEL::out() << "      relax = " << result << endl ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
AP_FluidStructureMaster:: save_for_post_processing( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: save_for_post_processing" ) ;

   static size_t last_iter_saved = PEL::bad_index() ;
   size_t const iter =
      TIME_IT->is_started() ? TIME_IT->iteration_number() : 0 ;
   bool save_iteration =
      !TIME_IT->is_started() ||
      TIME_IT->is_started() && TIME_IT->is_finished() ||
      greater_or_equal( TIME_IT->time(), graphics_next_time ) ;
   
   if( SAVE_ALL || 
       ( save_iteration && iter!=last_iter_saved ) )
   {
      for( size_t i=0 ; i<SDOMS->nb_domains() ; ++i )
      {
         do_save_for_post( SDOMS->domain( i )->result_saver() )  ;
      }
      last_iter_saved = iter ;
      graphics_next_time = next_time( graphics_times ) ;
   }
}

//----------------------------------------------------------------------
void
AP_FluidStructureMaster:: do_save_for_post( PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: do_save_for_post" ) ;

   rs->start_cycle() ;

   PEL::out() << endl << "   +++ SAVE FOR POSTPROCESSING" 
        << " *** CYCLE = " << rs->cycle_number()
        << " *** TIME = "  << TIME_IT->time() 
        << " ++++++" << endl << endl ;

   if( rs->attached_domain() == F_DOM )
   {
      if( TIME_IT->is_started() ) rs->save_grid() ;
      FGRID->do_additional_savings( TIME_IT, rs ) ;
      FLUID->do_additional_savings( TIME_IT, rs ) ;
      LOADC->do_additional_savings( TIME_IT, rs ) ;
   }
   if( rs->attached_domain() == S_DOM )
   {
      SOLID->do_additional_savings( TIME_IT, rs ) ;
   }

   if( SAVEFG )
   {
      if( !TIME_IT->is_started() ) rs->save_grid() ;
      rs->save_fields( graphics_level ) ;
   }
   rs->save_variable( TIME_IT->time(), "TIME" ) ;

   rs->terminate_cycle() ;
}

//----------------------------------------------------------------------
bool
AP_FluidStructureMaster:: greater_or_equal( double t1, double t2 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: greater_or_equal" ) ;
   bool result = ( t2 > 0.0  ? ( t1 >= t2*( 1.0 - 1.E-8 ) )
                   : ( t1 >= t2*( 1.0 + 1.E-8 ) ) ) ;
   return result ;
}

//----------------------------------------------------------------------
double
AP_FluidStructureMaster:: next_time( doubleVector const& dates ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: next_time" ) ;
   double result = PEL::max_double() ;
   for( size_t i=0 ; i<dates.size() ; i++ )
   {
      if( !greater_or_equal( TIME_IT->time(), dates(i) ) )
      {
         result = dates(i) ;
         break ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
void
AP_FluidStructureMaster:: check_times( std::string const& list_name,
                                        doubleVector const& dates ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: check_times" ) ;
   size_t n = dates.size() ;
   if( n>0 )
   {
      bool ok = greater_or_equal( dates(0), TIME_IT->initial_time() )
         && greater_or_equal( TIME_IT->final_time(), dates(n-1) ) ;
      for( size_t i=0 ; i<n-1 ; i++ )
      {
         ok = ok && dates(i+1)>dates(i) ;
      }
      if( !ok )
      {
         std::ostringstream mesg ;
         mesg << "AP_FluidStructureMaster :" << endl ;
         mesg << "   Invalid data of keyword : \""<< list_name <<"\"" << endl ;
         mesg << "   The time sequence must " << endl ;
         mesg << "      - start after the initial time ;" << endl ;
         mesg << "      - stop before the final time ;" << endl ;
         mesg << "      - be strictly increasing." ;
         PEL_Error::object()->raise_plain( mesg.str() ) ;
      }     
   }
}

//----------------------------------------------------------------------
void
AP_FluidStructureMaster:: save_for_restart( void )
//----------------------------------------------------------------------
{
   bool save_iteration = greater_or_equal( TIME_IT->time(), saver_next_time ) ;
   if( save_iteration )
   {
      write_storable_objects() ;
      saver_next_time = next_time( saver_times ) ;
   }
}

//----------------------------------------------------------------------
void
AP_FluidStructureMaster:: print_memory_info( void )
//----------------------------------------------------------------------
{
   static size_t mo = 1024*1024 ;
   static size_t go = 1024*1024*1024 ;
   // Memory in octets :
   size_t const mem = PEL_System::used_memory() ;
         
   PEL::out() << endl
              << "   Memory usage : "
              << mem
              << " octets" ;
   if( mem>go )
   {
      PEL::out() << " ( " << ( (double) mem )/go << " Go )" ;
   }
   else if( mem>mo )
   {
      PEL::out() << " ( " << ( (double) mem)/mo << " Mo )" ;
   }
   PEL::out() << endl
              << "   Number of objects : "
              << PEL_Object::GetNumberOf_PEL_objects()
              << std::endl ;
}

//----------------------------------------------------------------------------
void
AP_FluidStructureMaster:: check_field_nb_components( 
                                          PDE_DiscreteField const* ff,
                                          size_t required_nb_cmps ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: check_field_nb_components" ) ;
   PEL_CHECK( ff != 0 ) ;

   if( !( required_nb_cmps == ff->nb_components() ) )
   {
      ostringstream mesg ;
      mesg << "The field : \"" << ff->name() << "\"" << endl ;
      mesg << "should have " << required_nb_cmps << " components" << endl ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( ff->storage_depth() > required_nb_cmps ) ;
}

//----------------------------------------------------------------------
void
AP_FluidStructureMaster:: check_field_storage_depth( 
                                          PDE_DiscreteField const* ff,
                                          size_t requested_level ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: check_field_storage_depth" ) ;
   PEL_CHECK( ff != 0 ) ;

   if( !( requested_level < ff->storage_depth() ) )
   {
      ostringstream mesg ;
      mesg << "field : \"" << ff->name() << "\"" << endl ;
      mesg << "   the level of storage number " 
           << requested_level << " will be" << endl ;
      mesg << "   accessed in class " << type_name() << endl ;
      mesg << "   whereas the storage depth is only " 
           << ff->storage_depth() ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK( ff->storage_depth() > requested_level ) ;
}

//----------------------------------------------------------------------
bool
AP_FluidStructureMaster:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Application::invariant() ) ;

   return( true ) ;
}

//----------------------------------------------------------------------------
void
AP_FluidStructureMaster:: check_displacement_continuity( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_FluidStructureMaster:: check_displacement_continuity" ) ;

   size_t n = 0 ;

   boolVector done( FGD->nb_nodes() ) ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   { 
      GE_Color const* color = bFE->color() ;
      if( color->is_matching( INTERF_COL ) )
      {
         for( size_t loc_n=0 ; loc_n<bFE->nb_local_nodes( FGD ) ; ++loc_n )
         {
	    if( bFE->local_node_is_in_mesh( FGD, loc_n ) )
	    {
 	       GE_Point const* pt = bFE->local_node_location( FGD, loc_n ) ;
	       n = bFE->global_node( FGD, loc_n ) ;
               pt->print( cout, 0 ) ; cout << endl ;
	       if( ! done( n ) )
	       {
	          done( n ) = true ;
 	          for( size_t ic=0 ; ic<FGD->nb_components() ; ++ic )
		  {
                     if( !FGD->DOF_has_imposed_value( n, ic ) )
                     {
                        PEL_Error::object()->raise_plain( "error 0" ) ;
                     }
                     else
                     {
                        bFE->set_calculation_point( pt ) ;
                        PEL_ASSERT( bFE->nb_local_nodes( SD ) != 0 ) ;
                        double d_solid = bFE->value_at_pt( SD, L_NEW, ic ) ;
                        double d_grid = FGD->DOF_value( L_NEW, n, ic ) ;
                        cout << "   d_solid(" << ic << ")=" << d_solid
                             << "   d_grid(" << ic << ")=" << d_grid 
                             << "   delta=" << d_solid-d_grid << endl ;
                     }
	          }
	       }
	    }
	 }
      }
   }
}

