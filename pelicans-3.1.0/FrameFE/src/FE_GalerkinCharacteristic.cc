/*
 *  Copyright 1995-2010 by IRSN
 *
 *  This software is an application framework, with a set of integrated  
 *  reusable components, whose purpose is to simplify the task of developing 
 *  softwares of numerical mathematics and scientific computing.
 * 
 *  This software is governed by the CeCILL-C license under French law and 
 *  abiding by the rules of distribution of free software. You can use, modify 
 *  and/or redistribute the software under the terms of the CeCILL-C license  
 *  as circulated by CEA, CNRS and INRIA at the following URL 
 *  "http://www.cecill.info". 
 *
 *  As a counterpart to the access to the source code and rights to copy,  
 *  modify and redistribute granted by the license, users are provided only 
 *  with a limited warranty and the software's author, the holder of the  
 *  economic rights, and the successive licensors have only limited liability. 
 *
 *  In this respect, the user's attention is drawn to the risks associated  
 *  with loading, using, modifying and/or developing or reproducing the  
 *  software by the user in light of its specific status of free software,
 *  that may mean that it is complicated to manipulate, and that also  
 *  therefore means that it is reserved for developers and experienced 
 *  professionals having in-depth computer knowledge. Users are therefore 
 *  encouraged to load and test the software's suitability as regards their 
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and, more generally, to use and operate it in the same 
 *  conditions as regards security. 
 *
 *  The fact that you are presently reading this means that you have had 
 *  knowledge of the CeCILL-C license and that you accept its terms.
 */

#include <FE_GalerkinCharacteristic.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_ListIterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_String.hh>
#include <PEL_Timer.hh>
#include <PEL_Vector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_CFootFinder.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_GridFE.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <FE_Galerkin.hh>
#include <FE_TimeIterator.hh>

#include <iostream>
#include <iomanip>
#include <sstream>

using std::ios_base ; using std::setprecision ; using std::setw ;
using std::endl ;

FE_GalerkinCharacteristic const* 
FE_GalerkinCharacteristic:: PROTOTYPE = new FE_GalerkinCharacteristic() ;

struct FE_GalerkinCharacteristic_ERROR
{
   static void n0( PDE_LocalFEcell const* fe ) ;
   static void n1( void ) ;
} ;

//----------------------------------------------------------------------
FE_GalerkinCharacteristic:: FE_GalerkinCharacteristic( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_GalerkinCharacteristic" )
   , CFL_MEAN( 0 )
   , CFL_MIN( 0 )
   , CFL_MAX( 0 )
{
}

//----------------------------------------------------------------------
FE_GalerkinCharacteristic*
FE_GalerkinCharacteristic:: create_replica( 
                                   PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_GalerkinCharacteristic* result =
              new FE_GalerkinCharacteristic( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_GalerkinCharacteristic:: FE_GalerkinCharacteristic( 
                                          PEL_Object* a_owner,
                                          PDE_DomainAndFields const* dom,
                                          FE_SetOfParameters const* prms,
                                          PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , DOM( dom )
   , L_MASKED( PEL::bad_index() )
   , GLIST( PEL_List::create( this ) )
   , G_IT( 0 )
   , SORTED_GLIST( PEL_List::create( this ) )
   , S_G_IT( 0 )
   , TR_TIMER(  PEL_Timer::create( this ) )
   , CV_FIELDS( PEL_Vector::create( this, 0 ) )
   , ADV( dom->set_of_discrete_fields()->item( 
                               exp->string_data( "advective_field" ) ) )
   , L_ADV( exp->int_data( "advective_field_level" ) )
   , cFE( 0 )
   , CFL_MEAN( dom->nb_space_dimensions() )
   , CFL_MIN( dom->nb_space_dimensions() )
   , CFL_MAX( dom->nb_space_dimensions() )
   , NB_ITG_PTS( 0 )
   , ITG_WEIGHT( 0.0 )
   , NB_NO_FOUND_ITG_PTS( 0 )
   , W_NO_FOUND( 0.0 )
   , NB_INFLOWS( 0 )
   , W_INFLOW( 0.0 ) 
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: FE_GalerkinCharacteristic" ) ;

   check_field_storage_depth( ADV, L_ADV ) ;

   cFE = DOM->create_LocalFEcell( this ) ;

   PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "PDE_CFootFinder" ) ;
   cFE->set_foot_finder( e ) ;

   e->destroy() ;
   e = exp->create_subexplorer( 0, "list_of_FE_Galerkin" ) ;
   e->start_module_iterator() ;
   for( ; e->is_valid_module() ; e->go_next_module() )
   {
      PEL_ModuleExplorer* se = e->create_subexplorer( 0 ) ;
      
      FE_Galerkin* gpb = FE_Galerkin::make( GLIST, dom, prms, se ) ;
      insert_Galerkin_problem( gpb ) ;
      se->destroy() ;
   }
   e->destroy() ;

   G_IT = PEL_ListIterator::create( this, GLIST )  ;
   S_G_IT = PEL_ListIterator::create( this, SORTED_GLIST )  ;
}

//----------------------------------------------------------------------
FE_GalerkinCharacteristic:: ~FE_GalerkinCharacteristic( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: do_before_time_stepping( 
                                           FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   start_total_timer( "FE_GalerkinCharacteristic:: do_before_time_stepping" ) ;

   for( G_IT->start() ; G_IT->is_valid() ; G_IT->go_next() ) 
   {
      FE_Galerkin* pb = static_cast<FE_Galerkin*>( G_IT->item() ) ;
      pb->do_before_time_stepping( t_it ) ;
   }

   stop_total_timer() ;
}

//----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: do_before_inner_iterations_stage(
                                            FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: do_before_inner_iterations_stage" );
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   start_total_timer( 
         "FE_GalerkinCharacteristic:: do_before_inner_iterations_stage" ) ;

   for( G_IT->start(); G_IT->is_valid() ; G_IT->go_next() ) 
   {
      FE_Galerkin* pb = static_cast<FE_Galerkin*>( G_IT->item() ) ;
      pb->do_before_inner_iterations_stage( t_it ) ;
   }

   stop_total_timer() ;

   PEL_CHECK_POST( do_before_inner_iterations_stage_POST( t_it ) ) ;
}

//----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: do_one_inner_iteration( 
                                             FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   start_total_timer( "FE_GalerkinCharacteristic:: do_one_inner_iteration" ) ;
   // --------------

   for( G_IT->start(); G_IT->is_valid() ; G_IT->go_next() ) 
   {
      FE_Galerkin* pb = static_cast<FE_Galerkin*>( G_IT->item() ) ;
      pb->reset_discrete_problem( t_it ) ;
   }
   
   start_assembling_timer() ;
   // -------------------

   PDE_CFootFinder* finder = cFE->foot_finder() ;
   reset_statistics() ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      // - find the feet of the characteristics of the integration points
      // - mask the value of the fields with the value at those feet.
      double dtime = t_it->time_step() ;

      GE_QRprovider const* old_qrp = 0 ;
      for( S_G_IT->start() ; S_G_IT->is_valid() ; S_G_IT->go_next() )
      {
         FE_Galerkin* pb = static_cast<FE_Galerkin*>( S_G_IT->item() ) ;
         GE_QRprovider const* qrp = pb->QRprovider_for_material_derivative() ;
         if( qrp != old_qrp )
	 {
            old_qrp = qrp ;
            TR_TIMER->start() ;
            cFE->start_IP_iterator( qrp ) ; 
            for( ; cFE->valid_IP() ; cFE->go_next_IP() )
            {
               finder->search_foot( ADV, L_ADV, dtime,
                                                   cFE->coordinates_of_IP() ) ;
               if( !finder->is_ersatz_foot_investigated() &&
                   !finder->foot_has_been_found() )
               {
                  FE_GalerkinCharacteristic_ERROR::n0( cFE ) ;
               }
               update_statistics( finder ) ;
               for( size_t i=0 ; i<CV_FIELDS->count() ; ++i )
               {
                  PDE_DiscreteField* f =
                      static_cast<PDE_DiscreteField*>( CV_FIELDS->at(i) ) ;
                  for( size_t ic=0 ; ic<f->nb_components() ; ic++ )
                  {
                     double xx = -PEL::max_double() ;
                     if( finder->foot_has_been_found() )
                     {
                        xx = finder->value_at_foot( f, L_MASKED, ic ) ;
                     }
                     else if( finder->is_ersatz_foot_investigated() )
                     {
                        xx = finder->value_at_ersatz_foot( f, L_MASKED, ic ) ;
                     }
                     cFE->mask_value_at_IP( f, L_MASKED, xx, ic ) ;
                   }
                }
             }
             TR_TIMER->stop() ;
	 }

         pb->build_cell_contribution_to_material_derivative( t_it, cFE ) ;
      }
         
   }

   stop_assembling_timer() ;
   // ------------------

   if( verbose_level() >= 2 ) print_transport_statistics( PEL::out(), 5 ) ;

   for( G_IT->start() ; G_IT->is_valid() ; G_IT->go_next() )
   {
      FE_Galerkin* pb = static_cast< FE_Galerkin* >( G_IT->item() ) ;
      pb->terminate_discrete_problem( t_it ) ;
   }

   stop_total_timer() ;
   // -------------
}

//----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: do_after_inner_iterations_stage(
                                        FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: do_after_inner_iterations_stage" );
   PEL_CHECK_PRE( do_after_inner_iterations_stage_PRE( t_it ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( G_IT->start() ; G_IT->is_valid() ; G_IT->go_next() )
   {
      FE_Galerkin* pb = static_cast< FE_Galerkin* >( G_IT->item() ) ;
      pb->do_after_inner_iterations_stage( t_it ) ;
   }
}

//----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: do_after_time_stepping( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: do_after_time_stepping" ) ;

   start_total_timer( "FE_SplitSystem:: do_after_time_stepping" ) ;

   for( G_IT->start() ; G_IT->is_valid() ; G_IT->go_next() )
   {
      FE_Galerkin* pb = static_cast< FE_Galerkin* >( G_IT->item() ) ;
      pb->do_after_time_stepping() ;
   }

   stop_total_timer() ;
}

//----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: print_additional_times( std::ostream& os,
                                                    size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: print_additional_times" ) ;

   os << "FE_GalerkinCharacteristic: transport time = " ;
   TR_TIMER->print( os, 0 ) ;
   os << endl ;
   
   for( G_IT->start() ; G_IT->is_valid() ; G_IT->go_next() ) 
   {
      FE_Galerkin* pb = static_cast< FE_Galerkin* >( G_IT->item() ) ;
      pb->print_additional_times( os, indent_width ) ;
   }
}

//----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: do_additional_savings( FE_TimeIterator const* t_it,
                                                   PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: do_additional_savings" );
   PEL_CHECK_PRE( do_additional_savings_PRE( t_it, rs ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   start_total_timer( "FE_GalerkinCharacteristic:: do_additional_savings" ) ;

   for( G_IT->start() ; G_IT->is_valid() ; G_IT->go_next() ) 
   {
      FE_Galerkin* pb = static_cast< FE_Galerkin* >( G_IT->item() ) ;
      pb->do_additional_savings( t_it, rs ) ;
   }

   stop_total_timer() ;
}

//----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: add_storable_objects( PEL_ListIdentity* list )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: add_storable_objects" ) ;
   PEL_CHECK_PRE( add_storable_objects_PRE( list ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( G_IT->start() ; G_IT->is_valid() ; G_IT->go_next() ) 
   {
      FE_Galerkin* pb = static_cast< FE_Galerkin* >( G_IT->item() ) ;
      pb->add_storable_objects( list ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: print( std::ostream& os,
                                   size_t indent_width ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: print" ) ;

   FE_OneStepIteration:: print( os, indent_width ) ;

   for( G_IT->start() ; G_IT->is_valid() ; G_IT->go_next() )
   {
      FE_Galerkin const* pb = static_cast< FE_Galerkin* >( G_IT->item() ) ;
      pb->print( os, indent_width+3 ) ;
   }
}

//----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: insert_Galerkin_problem( FE_Galerkin* gpb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: insert_Galerkin_problem" ) ;

   if( L_MASKED != PEL::bad_index() )
   {
      if( L_MASKED != gpb->masked_level () )
         FE_GalerkinCharacteristic_ERROR::n1() ;
   }
   else
   {
      L_MASKED = gpb->masked_level() ;
   }

   // Update convected fields :
   PEL_Iterator* itf = gpb->convected_fields()->create_iterator( 0 ) ;
   for( itf->start() ; itf->is_valid() ; itf->go_next() )
   {
      CV_FIELDS->extend( itf->item() ) ;
   }
   itf->destroy() ;

   gpb->transfer_calculation_requirements_for_material_derivative( cFE ) ;

   // Append to solvers list :
   GLIST->append( gpb ) ;
   
   // Solvers sorted for material derivative assembling :
   GE_QRprovider const* qr = gpb->QRprovider_for_material_derivative() ;
   bool found = false ;
   for( size_t i=0 ; !found && i<SORTED_GLIST->count() ; ++i )
   {
      FE_Galerkin* pb = static_cast<FE_Galerkin*>( SORTED_GLIST->at(i) ) ;
      GE_QRprovider const* qr_pb = pb->QRprovider_for_material_derivative() ;
      if( qr==qr_pb )
      {
         found = true ;
         SORTED_GLIST->insert_at( i, gpb ) ;
      }
   }
   if( !found )
   {
      SORTED_GLIST->append( gpb ) ;
   }  
}

//-----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: reset_statistics( void )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: reset_statistics" ) ;
   PEL_CHECK( cFE!=0 && !cFE->is_valid() ) ;
   
   for( size_t dir=0 ; dir<CFL_MEAN.size() ; ++dir )
   {
      CFL_MEAN( dir ) = 0. ;
      CFL_MIN( dir ) =  PEL::max_double() ;
      CFL_MAX( dir ) = -PEL::max_double() ;
   }
   NB_ITG_PTS = 0 ;
   ITG_WEIGHT = 0. ;
   NB_NO_FOUND_ITG_PTS = 0 ;
   W_NO_FOUND = 0. ;
   NB_INFLOWS = 0 ;
   W_INFLOW = 0. ;
}

//-----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: update_statistics( PDE_CFootFinder const* finder )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: update_statistics" ) ;
   PEL_CHECK( finder!=0 ) ;
   PEL_CHECK( cFE!=0 && cFE->is_valid() && cFE->valid_IP() ) ;
   
   NB_ITG_PTS++ ;
   double const w = cFE->weight_of_IP() ;
   ITG_WEIGHT += w ;
   if( finder->foot_is_interior() )
   {
      for( size_t dim=0 ; dim<CFL_MEAN.size() ; ++dim )
      {
         double const cfl = finder->CFL_number( dim ) ;
         CFL_MEAN( dim ) += w*cfl ;
         if( cfl<CFL_MIN( dim ) )
         {
            CFL_MIN( dim ) = cfl ;
         }
         if( cfl>CFL_MAX( dim ) )
         {
            CFL_MAX( dim ) = cfl ;
         }
      }
   }
   else if( finder->foot_is_on_boundary() )
   {
      NB_INFLOWS++ ;
      W_INFLOW += w ;
   }
   else
   {
      NB_NO_FOUND_ITG_PTS++ ;
      W_NO_FOUND += w ;
   }
}

//-----------------------------------------------------------------------
void
FE_GalerkinCharacteristic:: print_transport_statistics(
                                             std::ostream& os,
                                             size_t indent_width ) const
//-----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GalerkinCharacteristic:: print_transport_statistics" ) ;

   std::string s( indent_width, ' ' ) ;
   s += indent() ;
   os << s << "Characteristic feet : "
           << NB_ITG_PTS
           << " integration points" << std::endl ;
   os << s << "Number_of_No_Found_Feet"
           << "  Number_of_InFlow_Feet"
      << std::endl ;
   std::ios_base::fmtflags original_flags = os.flags() ;
   os << std::setiosflags( std::ios::fixed ) ;
   std::streamsize p = os.precision() ;
   os << std::setprecision(2) ;
   os << s << std::setw(14) << NB_NO_FOUND_ITG_PTS
           << " (" << std::setw(5) << 100.*W_NO_FOUND/ITG_WEIGHT<< "%)"
           << std::setw(14) << NB_INFLOWS
           << " (" << std::setw(5) << 100.*W_INFLOW/ITG_WEIGHT << "%)"
           << std::endl ;
   if( CFL_MEAN.size()==1 )
   {
      os << s << "cfl_x_min  cfl_x_max      cfl_x  " << std::endl ;
   }
   else if( CFL_MEAN.size()==2 )
   {
      os << s << "cfl_x_min  cfl_x_max      cfl_x  "
               <<"cfl_y_min  cfl_y_max      cfl_y  " << std::endl ;
   }
   else if( CFL_MEAN.size()==3 )
   {
      os << s << "cfl_x_min  cfl_x_max      cfl_x  "
              << "cfl_y_min  cfl_y_max      cfl_y  "
              << "cfl_z_min  cfl_z_max      cfl_z  " << std::endl ;
   }
   os << s ;
   os << std::setiosflags( std::ios::scientific )
      << std::setprecision(3) ;
   for( size_t i=0 ; i<CFL_MEAN.size() ; ++i )
   {
      os << std::setw(9) << CFL_MIN(i) << "  "
         << std::setw(9) << CFL_MAX(i) << "  "
         << std::setw(9) << CFL_MEAN(i) << "  " ;
   }
   os << std::endl ;
   os << std::setprecision(p) ;
   os.flags( original_flags ) ;
}

//internal-------------------------------------------------------------------
void
FE_GalerkinCharacteristic_ERROR:: n0( PDE_LocalFEcell const* fe )
//internal-------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** FE_GalerkinCharacteristic error:" << endl ;
   mesg << "No characteristic foot found for the current integration point." 
        << std::endl ;
   mesg << "Integration point : " << std::endl ;
   fe->coordinates_of_IP()->print( mesg, 5 ) ;
   mesg << std::endl ;
   mesg << "Mesh : " << std::endl ;
   fe->polyhedron()->print( mesg, 5 ) ;
   mesg << std::endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_GalerkinCharacteristic_ERROR:: n1( void )
//internal-------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** FE_GalerkinCharacteristic error:" << endl ;
   mesg << "    all attached FE_Galerkin instances should" << std::endl ;
   mesg << "    return the same masked level" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

