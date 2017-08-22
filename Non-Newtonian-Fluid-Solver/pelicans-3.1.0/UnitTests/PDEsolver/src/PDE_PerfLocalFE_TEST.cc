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

#include <PDE_PerfLocalFE_TEST.hh>

#include <doubleVector.hh>
#include <stringVector.hh>

#include <PEL_assertions.hh>
#include <PEL_Exec.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Randomizer.hh>
#include <PEL_System.hh>
#include <PEL_Timer.hh>

#include <GE_QRprovider.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFE.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <ios>
#include <iostream>
#include <iomanip>

using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

PDE_PerfLocalFE_TEST const*
PDE_PerfLocalFE_TEST:: REGISTRATOR = new PDE_PerfLocalFE_TEST() ;

//---------------------------------------------------------------------------
PDE_PerfLocalFE_TEST:: PDE_PerfLocalFE_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_LocalFE", 
                     "PDE_PerfLocalFE_TEST" )
{
}

//----------------------------------------------------------------------------
PDE_PerfLocalFE_TEST:: ~PDE_PerfLocalFE_TEST( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
PDE_PerfLocalFE_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PerfLocalFE_TEST:: process_one_test" ) ;

   TEST_NAME = exp->name() ;

   PEL_ModuleExplorer* ee = 
                       exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PEL_Timer* timer = PEL_Timer::create( 0 ) ;

   size_t before = PEL_System::used_memory() ;
   print_memory_result( TEST_NAME+" PDE_DomainAndFields::create before", 
                        before ) ;
   timer->start() ;
   PDE_DomainAndFields* dom = PDE_DomainAndFields:: create( 0, ee ) ;
   timer->stop() ;
   size_t after = PEL_System::used_memory() ;
   print_memory_result( TEST_NAME+" PDE_DomainAndFields::create  after", 
                        after ) ;
   print_memory_result( TEST_NAME+" PDE_DomainAndFields::create   used", 
                        after-before ) ;
   print_time_result( TEST_NAME+" PDE_DomainAndFields::create", 
                      timer->time() ) ;
 
   timer->destroy() ; timer = 0 ;
   ee->destroy() ; ee = 0 ;

   PDE_SetOfDiscreteFields* sdf = dom->set_of_discrete_fields() ;

   if( exp->has_module( "cell_go_next" ) )
   {
      ee = exp->create_subexplorer( 0, "cell_go_next" ) ;
      PDE_LocalFEcell* fe = dom->create_LocalFEcell( 0 ) ;
      for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
      {
         PDE_DiscreteField const* ff = sdf->item() ;
         fe->require_field_calculation( ff, PDE_LocalFE::node ) ;
      }
      do_mesh_go_next( fe, ee, TEST_NAME+" cell_go_next" ) ;
      fe->destroy() ; fe = 0 ;
      ee->destroy() ; ee = 0 ;
   }

   if( exp->has_module( "bound_go_next" ) )
   {
      ee = exp->create_subexplorer( 0, "bound_go_next" ) ;
      PDE_LocalFEbound* fe = dom->create_LocalFEbound( 0 ) ;
      for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
      {
         PDE_DiscreteField const* ff = sdf->item() ;
         fe->require_field_calculation( ff, PDE_LocalFE::node ) ;
      }
      do_mesh_go_next( fe, ee, TEST_NAME+" bound_go_next" ) ;
      fe->destroy() ; fe = 0 ;
      ee->destroy() ; ee = 0 ;
   }

   if( exp->has_module( "cell_go_i_th" ) )
   {
      ee = exp->create_subexplorer( 0, "cell_go_i_th" ) ;
      PDE_LocalFEcell* fe = dom->create_LocalFEcell( 0 ) ;
      for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
      {
         PDE_DiscreteField const* ff = sdf->item() ;
         fe->require_field_calculation( ff, PDE_LocalFE::node ) ;
      }
      do_mesh_go_i_th( fe, ee, TEST_NAME+" cell_go_i_th" ) ;
      fe->destroy() ; fe = 0 ;
      ee->destroy() ; ee = 0 ;
   }

   if( exp->has_module( "bound_go_i_th" ) )
   {
      ee = exp->create_subexplorer( 0, "bound_go_i_th" ) ;
      PDE_LocalFEbound* fe = dom->create_LocalFEbound( 0 ) ;
      for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
      {
         PDE_DiscreteField const* ff = sdf->item() ;
         fe->require_field_calculation( ff, PDE_LocalFE::node ) ;
      }
      do_mesh_go_i_th( fe, ee, TEST_NAME+" bound_go_i_th" ) ;
      fe->destroy() ; fe = 0 ;
      ee->destroy() ; ee = 0 ;
   }

   if( exp->has_module( "cell_iteration_on_IPs" ) )
   {
      ee = exp->create_subexplorer( 0, "cell_iteration_on_IPs" ) ;

      bool ok = ee->bool_data( "test_N" ) ;
      if( ok )
      {
         PDE_LocalFEcell* fe = dom->create_LocalFEcell( 0 ) ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
            fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
         }
         do_mesh_iteration_on_IPs( fe, ee, 
                                   TEST_NAME+" cell_iteration_on_IPs_N" ) ;
         fe->destroy() ; fe = 0 ;
      }

      ok = ee->bool_data( "test_dN" ) ;
      if( ok )
      {
         PDE_LocalFEcell* fe = dom->create_LocalFEcell( 0 ) ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
            fe->require_field_calculation( ff, PDE_LocalFE::dN ) ;
         }
         do_mesh_iteration_on_IPs( fe, ee, 
                                   TEST_NAME+" cell_iteration_on_IPs_dN" ) ;
         fe->destroy() ; fe = 0 ;
      }

      ok = ee->bool_data( "test_N_dN" ) ;
      if( ok )
      {
         PDE_LocalFEcell* fe = dom->create_LocalFEcell( 0 ) ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
            fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
            fe->require_field_calculation( ff, PDE_LocalFE::dN ) ;
         }
         do_mesh_iteration_on_IPs( fe, ee, 
                                   TEST_NAME+" cell_iteration_on_IPs_N_dN" ) ;
         fe->destroy() ; fe = 0 ;
      }

      ok = ee->bool_data( "test_N_dN_d2N" ) ;
      if( ok )
      {
         PDE_LocalFEcell* fe = dom->create_LocalFEcell( 0 ) ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
            fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
            fe->require_field_calculation( ff, PDE_LocalFE::dN ) ;
            fe->require_field_calculation( ff, PDE_LocalFE::d2N ) ;
         }
         do_mesh_iteration_on_IPs( fe, ee, 
                               TEST_NAME+" cell_iteration_on_IPs_N_dN_d2N" ) ;
         fe->destroy() ; fe = 0 ;
      }

      ee->destroy() ; ee = 0 ;
   }

   if( exp->has_module( "bound_iteration_on_IPs" ) )
   {
      ee = exp->create_subexplorer( 0, "bound_iteration_on_IPs" ) ;

      bool ok = ee->bool_data( "test_N" ) ;
      if( ok )
      {
         PDE_LocalFEbound* fe = dom->create_LocalFEbound( 0 ) ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
               fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
         }
         do_mesh_iteration_on_IPs( fe, ee, 
                                   TEST_NAME+" bound_iteration_on_IPs_N" ) ;
         fe->destroy() ; fe = 0 ;
      }

      ok = ee->bool_data( "test_dN" ) ;
      if( ok )
      {
         PDE_LocalFEbound* fe = dom->create_LocalFEbound( 0 ) ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
            if( !dom->is_defined_on_boundary( ff->name() ) )
            {
               fe->require_field_calculation( ff, PDE_LocalFE::dN ) ;
            }
         }
         do_mesh_iteration_on_IPs( fe, ee,
                                   TEST_NAME+" bound_iteration_on_IPs_dN" ) ;
         fe->destroy() ; fe = 0 ;
      }

      ok = ee->bool_data( "test_N_dN" ) ;
      if( ok )
      {
         PDE_LocalFEbound* fe = dom->create_LocalFEbound( 0 ) ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
            fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
            if( !dom->is_defined_on_boundary( ff->name() ) )
            {
               fe->require_field_calculation( ff, PDE_LocalFE::dN ) ;
            }
         }
         do_mesh_iteration_on_IPs( fe, ee,
                                   TEST_NAME+" bound_iteration_on_IPs_N_dN" ) ;
         fe->destroy() ; fe = 0 ;
      }

      ok = ee->bool_data( "test_N_dN_d2N" ) ;
      if( ok )
      {
         PDE_LocalFEbound* fe = dom->create_LocalFEbound( 0 ) ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
            fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
            if( !dom->is_defined_on_boundary( ff->name() ) )
            {
               fe->require_field_calculation( ff, PDE_LocalFE::dN ) ;
               fe->require_field_calculation( ff, PDE_LocalFE::d2N ) ;
            }
         }
         do_mesh_iteration_on_IPs( fe, ee,
                              TEST_NAME+" bound_iteration_on_IPs_N_dN_d2N" ) ;
         fe->destroy() ; fe = 0 ;
      }

      ee->destroy() ; ee = 0 ;
   }

   dom->destroy() ;
}

//-----------------------------------------------------------------------------
void
PDE_PerfLocalFE_TEST:: do_mesh_go_next( PDE_LocalFE* fe,
                                        PEL_ModuleExplorer const* exp,
                                        std::string const& mesg )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PerfLocalFE_TEST:: do_mesh_go_next" ) ;

   size_t nb_reps = exp->int_data( "nb_repetitions" ) ;
   
   PEL_Timer* timer = PEL_Timer::create( 0 ) ;
   timer->start() ;   
   for( size_t i=0 ; i<nb_reps ; ++i )
   {
      for( fe->start() ; fe->is_valid() ; fe->go_next() )
      {
      }
   }
   timer->stop() ;

   print_time_result( mesg, timer->time() ) ;

   timer->destroy() ;
}

//-----------------------------------------------------------------------------
void
PDE_PerfLocalFE_TEST:: do_mesh_go_i_th( PDE_LocalFE* fe,
                                        PEL_ModuleExplorer const* exp,
                                        std::string const& mesg )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PerfLocalFE_TEST:: do_mesh_go_i_th" ) ;

   PEL_Randomizer* bill = PEL_Randomizer::create( 0, 1 ) ;
   bill->start() ;
   size_t nbm = fe->nb_meshes() ;
   
   size_t nb2 = exp->int_data( "nb_access" ) ;
   PEL_Timer* timer = PEL_Timer::create( 0 ) ;
   timer->start() ;   
   for( size_t i=0 ; i<nb2 ; ++i )
   {
      bill->go_next() ;
      fe->go_i_th( (size_t)(nbm*bill->item()) ) ;
   }
   timer->stop() ;

   print_time_result( mesg, timer->time() ) ;

   timer->destroy() ;
   bill->destroy() ;
}

//-----------------------------------------------------------------------------
void
PDE_PerfLocalFE_TEST:: do_mesh_iteration_on_IPs( PDE_LocalFE* fe,
                                                 PEL_ModuleExplorer const* exp,
                                                 std::string const& mesg )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PerfLocalFE_TEST:: do_mesh_iteration_on_IPs" ) ;

   GE_QRprovider const* qrp = 
      GE_QRprovider::object( exp->string_data( "quadrature_rule_provider" ) ) ;
   size_t nb_reps = exp->int_data( "nb_repetitions" ) ;
   
   PEL_Timer* timer = PEL_Timer::create( 0 ) ;
   for( size_t i=0 ; i<nb_reps ; ++i )
   {
      for( fe->start() ; fe->is_valid() ; fe->go_next() )
      {
         timer->start() ;   
         fe->start_IP_iterator( qrp ) ;
         for( ; fe->valid_IP() ; fe->go_next_IP() )
         {
         }
         timer->stop() ;
      }
   }

   print_time_result( mesg, timer->time() ) ;

   timer->destroy() ;
}

