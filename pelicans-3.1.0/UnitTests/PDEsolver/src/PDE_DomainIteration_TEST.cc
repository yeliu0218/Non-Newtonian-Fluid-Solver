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

#include <PDE_DomainIteration_TEST.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>
#include <boolVector.hh>
#include <stringVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Transform.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_InterfaceAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalFEmortarSide.hh>
#include <PDE_CursorFEside.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfDomains.hh>

#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>

using std::cout ;
using std::endl ;
using std::ios_base ;
using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

//---------------------------------------------------------------------------
PDE_DomainIteration_TEST*
PDE_DomainIteration_TEST:: registered_test = new PDE_DomainIteration_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_DomainIteration_TEST:: PDE_DomainIteration_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_LocalFE", "PDE_DomainIteration_TEST" )
   , QRP( 0 )
{
}

//---------------------------------------------------------------------------
PDE_DomainIteration_TEST:: ~PDE_DomainIteration_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: process_one_test" ) ;

   out() << "| ... "  <<  exp->name() << endl ;
   
   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "iterations_trace" ) ;
   std::ofstream ofs( se->string_data( "output_file" ).c_str() ) ;
   se->destroy() ;

   bool do_print_IPs = false ;
   bool do_print_vals_at_IPs = false ;

   if( exp->has_module( "calculations_consistency" ) )
   {
      se = exp->create_subexplorer( 0, "calculations_consistency" ) ;
      D_EPS = se->double_data( "dbl_epsilon" ) ;
      D_MIN = se->double_data( "dbl_minimum" ) ;
      QRP = GE_QRprovider::object( 
                se->string_data( "quadrature_rule_provider" ) ) ;
      if( se->has_entry( "do_print_IPs" ) )
         do_print_IPs = se->bool_data( "do_print_IPs" ) ;
      if( se->has_entry( "do_print_values_at_IPs" ) )
         do_print_vals_at_IPs = se->bool_data( "do_print_values_at_IPs" ) ;
      se->destroy() ;
   }

   if( exp->has_module( "PDE_SetOfDomains" ) )
   {
      se = exp->create_subexplorer( 0, "PDE_SetOfDomains" ) ;
      PDE_SetOfDomains* sdoms = PDE_SetOfDomains::create( 0, se ) ;
      for( size_t i=0 ; i<sdoms->nb_domains() ; ++i )
      {
         PDE_DomainAndFields* dom = sdoms->domain( i ) ;
         process_domain( dom, ofs, do_print_IPs, do_print_vals_at_IPs ) ;
      }
      for( size_t i=0 ; i<sdoms->nb_interfaces() ; ++i )
      {
         PDE_InterfaceAndFields* interf = sdoms->interface( i ) ;
         PDE_DomainAndFields* dom_0 = 
           sdoms->domain( sdoms->index_of_interface_adjacent_domain( i, 0 ) ) ;
         PDE_DomainAndFields* dom_1 = 
           sdoms->domain( sdoms->index_of_interface_adjacent_domain( i, 1 ) ) ;
         process_interface( interf, dom_0, dom_1, 
                            ofs, do_print_IPs, do_print_vals_at_IPs ) ;
      }
      sdoms->destroy() ;
      se->destroy() ;
   }

   if( exp->has_module( "PDE_DomainAndFields" ) )
   {
      se = exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
      PDE_DomainAndFields* dom = PDE_DomainAndFields::create( 0, se ) ;
      process_domain( dom, ofs, do_print_IPs, do_print_vals_at_IPs ) ;
      dom->destroy() ;
      se->destroy() ;
   }

   ofs.close() ;
}

//---------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: process_domain( PDE_DomainAndFields const* dom,
                                           std::ostream& os,
                                           bool do_print_IPs,
                                           bool do_print_vals_at_IPs )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: process_domain" ) ;

   PDE_SetOfDiscreteFields* sdf = dom->set_of_discrete_fields() ;

   PEL_ModuleExplorer const* exp = dom->decohesion_explorer() ;
   if( exp != 0 ) exp->print( os, 3 ) ;

   PDE_LocalFEcell*  cfe_1 = dom->create_LocalFEcell( 0 ) ;
   cfe_1->include_color( GE_Color::halo_color() ) ;
   PDE_LocalFEcell*  cfe_2 = dom->create_LocalFEcell( 0 ) ;
   cfe_2->include_color( GE_Color::halo_color() ) ;
   PDE_LocalFEbound* bfe_1 = dom->create_LocalFEbound( 0 ) ;
   bfe_1->include_color( GE_Color::halo_color() ) ;
   PDE_LocalFEbound* bfe_2 = dom->create_LocalFEbound( 0 ) ;
   bfe_2->include_color( GE_Color::halo_color() ) ;
   PDE_CursorFEside* sfe_1 = dom->create_CursorFEside( 0 ) ;
   sfe_1->include_color( GE_Color::halo_color() ) ;
   PDE_CursorFEside* sfe_2 = dom->create_CursorFEside( 0 ) ;
   sfe_2->include_color( GE_Color::halo_color() ) ;
   
   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      bfe_1->require_field_calculation( ff, PDE_LocalFE::N ) ;
      bfe_2->require_field_calculation( ff, PDE_LocalFE::N ) ;
      if( !dom->is_defined_on_boundary( ff->name() ) )
      {
         cfe_1->require_field_calculation( ff, PDE_LocalFE::N ) ;
         cfe_1->require_field_calculation( ff, PDE_LocalFE::dN ) ;
         cfe_1->require_field_calculation( ff, PDE_LocalFE::d2N ) ;
         cfe_2->require_field_calculation( ff, PDE_LocalFE::N ) ;
         bfe_1->require_field_calculation( ff, PDE_LocalFE::dN ) ;
         bfe_1->require_field_calculation( ff, PDE_LocalFE::d2N ) ;
         sfe_1->require_field_calculation( ff, PDE_LocalFE::node ) ;
         sfe_2->require_field_calculation( ff, PDE_LocalFE::node ) ;
      }
   }
   for( size_t i=0 ; i<dom->nb_conformal_adjacent_domains() ; ++i )
   {
      PDE_DomainAndFields const* adjdom = dom->conformal_adjacent_domain( i ) ;
      PDE_SetOfDiscreteFields* adjsdf = adjdom->set_of_discrete_fields() ;
      for( adjsdf->start() ; adjsdf->is_valid() ;adjsdf->go_next() )
      {
         PDE_DiscreteField const* ff = adjsdf->item() ;
         bfe_1->require_field_calculation( ff, PDE_LocalFE::N ) ;
         bfe_2->require_field_calculation( ff, PDE_LocalFE::N ) ;
         if( !adjdom->is_defined_on_boundary( ff->name() ) )
         {
            bfe_1->require_field_calculation( ff, PDE_LocalFE::dN ) ;
            bfe_1->require_field_calculation( ff, PDE_LocalFE::d2N ) ;
         }
      }
   }
   
   iterate_over_meshes( cfe_1, cfe_2, 
                        os, do_print_IPs, do_print_vals_at_IPs ) ;
   
   iterate_over_cells( cfe_1, sfe_1, bfe_1 ) ;

   iterate_over_meshes( bfe_1, bfe_2, 
                        os, do_print_IPs, do_print_vals_at_IPs ) ;

   iterate_over_sides( sfe_1, sfe_2, dom, os ) ;
   
   cfe_1->destroy() ; cfe_1 = 0 ;
   cfe_2->destroy() ; cfe_2 = 0 ;
   bfe_1->destroy() ; bfe_1 = 0 ;
   bfe_2->destroy() ; bfe_2 = 0 ;
   sfe_1->destroy() ; sfe_1 = 0 ;
   sfe_2->destroy() ; sfe_2 = 0 ;
}

//---------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: process_interface( 
                                   PDE_InterfaceAndFields const* interf,
                                   PDE_DomainAndFields const* dom_0,
                                   PDE_DomainAndFields const* dom_1,
                                   std::ostream& os,
                                   bool do_print_IPs,
                                   bool do_print_vals_at_IPs )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: process_interface" ) ;

   PDE_LocalFEmortarSide* fe_1 = interf->create_LocalFEmortarSide( 0 ) ;
   PDE_LocalFEmortarSide* fe_2 = interf->create_LocalFEmortarSide( 0 ) ;

   PDE_SetOfDiscreteFields* sdf = interf->set_of_discrete_fields() ;
   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      fe_1->require_field_calculation( sdf->item(), PDE_LocalFE::N ) ;
      fe_2->require_field_calculation( sdf->item(), PDE_LocalFE::N ) ;
   }

   sdf = dom_0->set_of_discrete_fields() ;
   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      fe_1->require_field_calculation( ff, PDE_LocalFE::N ) ;
      if( !dom_0->is_defined_on_boundary( ff->name() ) )
      {
         fe_1->require_field_calculation( ff, PDE_LocalFE::dN ) ;
         fe_1->require_field_calculation( ff, PDE_LocalFE::d2N ) ;
      }
      fe_2->require_field_calculation( ff, PDE_LocalFE::N ) ;
   }

   sdf = dom_1->set_of_discrete_fields() ;
   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      fe_1->require_field_calculation( sdf->item(), PDE_LocalFE::N ) ;
      fe_2->require_field_calculation( sdf->item(), PDE_LocalFE::N ) ;
   }

   iterate_over_meshes( fe_1, fe_2, os, 
                        do_print_IPs, do_print_vals_at_IPs ) ;

   fe_1->destroy() ;
   fe_2->destroy() ;
}

//---------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: iterate_over_meshes( PDE_LocalFE* fe_1,
                                                PDE_LocalFE* fe_2,
                                                std::ostream& os,
                                                bool do_print_IPs,
                                                bool do_print_vals_at_IPs )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: iterate_over_meshes"  ) ;

   bool ok_go_i_th = true ;
   bool ok_node_loc = true ;
   bool ok_val = true ;
   bool done_der1 = false ; bool ok_der1 = true ;
   bool done_der2 = false ; bool ok_der2 = true ;

   os << "***************************************" << endl ;
   os << fe_1->type_name() <<  endl ;
   os << "   " << fe_1->nb_meshes() << " meshes" << endl ;
   os << "***************************************" << endl ;
   os << "---------------------------------------" << endl ;
   os << "handled fields" << endl ;
   fe_1->print_handled_fields( os, 3 ) ;
   for( fe_1->start() ; fe_1->is_valid() ; fe_1->go_next() )
   {
      os << "---------------------------------------" << endl ;
      os << "current mesh" << endl ;
      fe_1->print_current_mesh( os, 3 ) ;
      os << "Local discretization of fields" << endl ;
      fe_1->print_local_discretization_of_fields( os, 3 ) ;

      fe_2->go_i_th( fe_1->mesh_id() ) ;
      ok_go_i_th = ok_go_i_th && ( fe_2->is_valid() ) ;
      ok_go_i_th = ok_go_i_th && ( fe_2->mesh_id() == fe_1->mesh_id() ) ;

      check_node_location( fe_1, ok_node_loc ) ;

      if( QRP != 0 ) iterate_over_IPs( fe_1, ok_val, 
                                       done_der1, ok_der1, 
                                       done_der2, ok_der2,
                                       os, do_print_IPs,
                                       do_print_vals_at_IPs ) ;
   }
   os << "---------------------------------------" << endl ;

   notify_one_test_result( fe_1->type_name()+" go_i_th/go_next", ok_go_i_th ) ;
   notify_one_test_result( fe_1->type_name()+" node_location", ok_node_loc ) ;
   notify_one_test_result( fe_1->type_name()+" values", ok_val ) ;
   if( done_der1 )
      notify_one_test_result( fe_1->type_name()+" derivatives_1", ok_der1 ) ;
   if( done_der2 )
      notify_one_test_result( fe_1->type_name()+" derivatives_2", ok_der2 ) ;
}

//----------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: iterate_over_cells( PDE_LocalFEcell*  cfe,
                                               PDE_CursorFEside* sfe,
                                               PDE_LocalFEbound* bfe )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: iterate_over_cells" ) ;

   bool ok_sides = true ;
   bool ok_bounds = true ;
   
   for( cfe->start() ; cfe->is_valid() ; cfe->go_next() )
   {
      size_t_vector const& side_ids = cfe->adjacent_side_ids() ;
      for( size_t i=0 ; i<side_ids.size() ; ++i )
      {
         sfe->go_i_th( side_ids( i ) ) ;
         ok_sides &=
               sfe->adjacent_localFEcell( 0 )->mesh_id() == cfe->mesh_id()
            || sfe->adjacent_localFEcell( 1 )->mesh_id() == cfe->mesh_id() ;
      }
      size_t_vector const& bound_ids = cfe->adjacent_bound_ids() ;
      for( size_t i=0 ; i<bound_ids.size() ; ++i )
      {
         bfe->go_i_th( bound_ids( i ) ) ;
         ok_bounds &= bfe->adjacent_cell_id() == cfe->mesh_id() ;
      }
   }
   notify_one_test_result( "PDE_LocalFEcell::adjacent_side_ids", ok_sides ) ;
   notify_one_test_result( "PDE_LocalFEcell::adjacent_bound_ids", ok_bounds ) ;
}

//----------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: check_node_location( PDE_LocalFE* fe,
                                                bool& ok_node_loc ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: check_node_location" ) ;

   GE_Mpolyhedron const* poly = fe->polyhedron() ;
   for( size_t iF=0 ; iF<fe->nb_handled_fields() ; ++iF )
   {
      PDE_DiscreteField const* ff = fe->handled_field( iF ) ;
      if( fe->nb_local_nodes( ff ) != 0 )
      {
         for( size_t i=0 ; i<fe->nb_local_nodes( ff ) ; ++i )
         {
            GE_Point const* pt = fe->local_node_location( ff, i ) ;
            bool ok = 
              ( poly->contains( pt ) &&  fe->local_node_is_in_mesh( ff, i )) ||
              (!poly->contains( pt ) && !fe->local_node_is_in_mesh( ff, i ) ) ;
            if( !ok )
            {
               out() << "------------------------------------" << endl ;
               out() << "current mesh" << endl ;
               fe->print_current_mesh( out(), 3 ) ;
               out() << "Local discretization of fields" << endl ;
               fe->print_local_discretization_of_fields( out(), 3 ) ;
               out() << "------" << endl ;
               out() << "node : " ; 
               pt->print( out(), 3 ) ; out() << endl ;
               out() << "   local_node_is_in_mesh --> "
                     << fe->local_node_is_in_mesh( ff, i ) << endl ;
               out() << "                contains --> "
                     << poly->contains( pt ) << endl ;
            }
            ok_node_loc = ok_node_loc && ok ;
         }
      }
   }
}

//----------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: iterate_over_IPs( 
                                PDE_LocalFE* fe,
                                bool& ok_val,
                                bool& done_der1, bool& ok_der1,
                                bool& done_der2, bool& ok_der2,
                                std::ostream& os,
                                bool do_print_IPs,
                                bool do_print_vals_at_IPs ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: iterate_over_IPs" ) ;

   if( do_print_IPs )
   {
      os << "Integration points" << endl ;
      fe->start_IP_iterator( QRP ) ;
      for( ; fe->valid_IP() ; fe->go_next_IP() )
      {
         fe->print_current_IP( os, 3 ) ;
      }
   }

   for( size_t iF=0 ; iF<fe->nb_handled_fields() ; ++iF )
   {
      PDE_DiscreteField const* ff = fe->handled_field( iF ) ;
      if( fe->nb_local_nodes( ff ) != 0 )
      {
         fe->set_row_and_col_fields( ff, ff ) ;
         if( do_print_vals_at_IPs )
            os << "Integration points" << endl ;
         fe->start_IP_iterator( QRP ) ;
         for( ; fe->valid_IP() ; fe->go_next_IP() )
         {
            if( do_print_vals_at_IPs )
               fe->print_values_at_current_IP( os, 3 ) ;

            check_IP_vals( fe, ff, ok_val ) ;

            if( fe->field_calculation_is_handled( ff, PDE_LocalFE::dN ) )
            {
               done_der1 = true ;
               check_IP_der1( fe, ff, ok_der1 ) ;
            }

            if( fe->field_calculation_is_handled( ff, PDE_LocalFE::d2N ) )
            {
               done_der2 = true ;
               check_IP_der2( fe, ff, ok_der2 ) ;
            }
         }
      }
   }
}

//----------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: check_IP_vals( PDE_LocalFE* fe,
                                          PDE_DiscreteField const* ff,
                                          bool& ok ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: check_IP_vals" ) ;

   bool eq ;
   GE_Point const* pt = fe->coordinates_of_IP() ;
   fe->set_calculation_point( pt ) ;
   for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
   {
      double ref1 = fe->value_at_IP( ff, 0, ic ) ;

      double ref2 = fe->value_at_pt( ff, 0, ic ) ;

      eq = PEL::double_equality( ref1, ref2, D_EPS, D_MIN ) ;
      if( !eq ) display_error( fe->type_name()+","+ff->name(), 
                               " value_at_pt", "value_at_IP", 
                               ref1, ref2 ) ;
      ok = ok && eq ;

      double val1 = 0.0 ;
      for( size_t il=0 ; il<fe->nb_local_nodes( ff ) ; ++il )
      {
         val1 += ff->DOF_value( 0, fe->global_node( ff, il ), ic ) *
                 fe->N_at_IP( row, il ) ;
      }
      eq = PEL::double_equality( val1, ref1, D_EPS, D_MIN ) ;
      if( !eq ) display_error( fe->type_name()+","+ff->name(), 
                               "sum_i u_i N_i", "value_at_IP", 
                               val1, ref1 ) ;
      ok = ok && eq ;

      double val2 = 0.0 ;
      for( size_t il=0 ; il<fe->nb_local_nodes( ff ) ; ++il )
      {
         val2 += ff->DOF_value( 0, fe->global_node( ff, il ), ic ) *
                 fe->N_at_pt( ff, il ) ;
      }
      eq = PEL::double_equality( val2, ref2, D_EPS, D_MIN ) ;
      if( !eq ) display_error( fe->type_name()+","+ff->name(), 
                               "sum_i u_i N_i", "value_at_pt", 
                               val2, ref2 ) ;
      ok = ok && eq ;
   }
}

//----------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: check_IP_der1( PDE_LocalFE* fe,
                                          PDE_DiscreteField const* ff,
                                          bool& ok ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: check_IP_der1" ) ;

   bool eq ;
   GE_Point const* pt = fe->coordinates_of_IP() ;
   fe->set_calculation_point( pt ) ;
   for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
   {
      for( size_t a=0 ; a<fe->nb_space_dimensions() ; ++a )
      {
         double ref1 = fe->gradient_at_IP( ff, 0, a, ic ) ;

         double ref2 = fe->gradient_at_pt( ff, 0, a, ic ) ;

         eq = PEL::double_equality( ref1, ref2, D_EPS, D_MIN ) ;
         if( !eq ) display_error( fe->type_name()+","+ff->name(), 
                                  " gradient_at_pt", "gradient_at_IP", 
                                  ref1, ref2 ) ;

         ok = ok && eq ;

         double val1 = 0.0 ;
         for( size_t il=0 ; il<fe->nb_local_nodes( ff ) ; ++il )
         {
            val1 += ff->DOF_value( 0, fe->global_node( ff, il ), ic ) *
                    fe->dN_at_IP( row, il, a ) ;
         }
         eq = PEL::double_equality( val1, ref1, D_EPS, D_MIN ) ;
         if( !eq ) display_error( fe->type_name()+","+ff->name(), 
                                  "sum_i u_i dN_i", "gradient_at_IP", 
                                  val1, ref1 ) ;
         ok = ok && eq ;

         double val2 = 0.0 ;
         for( size_t il=0 ; il<fe->nb_local_nodes( ff ) ; ++il )
         {
            val2 += ff->DOF_value( 0, fe->global_node( ff, il ), ic ) *
                    fe->dN_at_pt( ff, il, a ) ;
         }
         eq = PEL::double_equality( val2, ref2, D_EPS, D_MIN ) ;
         if( !eq ) display_error( fe->type_name()+","+ff->name(), 
                                  "sum_i u_i dN_i", "gradient_at_pt", 
                                  val2, ref2 ) ;
         ok = ok && eq ;
      }
   }
}

//----------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: check_IP_der2( PDE_LocalFE* fe,
                                          PDE_DiscreteField const* ff,
                                          bool& ok ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: check_IP_der2" ) ;

   bool eq ;
   GE_Point const* pt = fe->coordinates_of_IP() ;
   fe->set_calculation_point( pt ) ;
   for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
   {
      for( size_t a=0 ; a<fe->nb_space_dimensions() ; ++a )
      {
         for( size_t b=0 ; b<fe->nb_space_dimensions() ; ++b )
         {
            double ref1 = fe->hessian_at_IP( ff, 0, a, b, ic ) ;

            double ref2 = fe->hessian_at_pt( ff, 0, a, b, ic ) ;

            eq = PEL::double_equality( ref1, ref2, D_EPS, D_MIN ) ;
            if( !eq ) display_error( fe->type_name()+","+ff->name(), 
                                     " hessian_at_pt", "hessian_at_IP", 
                                     ref1, ref2 ) ;

            ok = ok && eq ;

            double val1 = 0.0 ;
            for( size_t il=0 ; il<fe->nb_local_nodes( ff ) ; ++il )
            {
               val1 += ff->DOF_value( 0, fe->global_node( ff, il ), ic ) *
                       fe->d2N_at_IP( row, il, a, b ) ;
            }
            eq = PEL::double_equality( val1, ref1, D_EPS, D_MIN ) ;
            if( !eq ) display_error( fe->type_name()+","+ff->name(), 
                                     "sum_i u_i d2N_i", "hessian_at_IP", 
                                     val1, ref1 ) ;
            ok = ok && eq ;

            double val2 = 0.0 ;
            for( size_t il=0 ; il<fe->nb_local_nodes( ff ) ; ++il )
            {
               val2 += ff->DOF_value( 0, fe->global_node( ff, il ), ic ) *
                       fe->d2N_at_pt( ff, il, a, b ) ;
            }
            eq = PEL::double_equality( val2, ref2, D_EPS, D_MIN ) ;
            if( !eq ) display_error( fe->type_name()+","+ff->name(), 
                                     "sum_i u_i d2N_i", "hessian_at_pt", 
                                     val2, ref2 ) ;
            ok = ok && eq ;
         }
      }
   }
}

//----------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: iterate_over_sides( PDE_CursorFEside* sfe_1,
                                               PDE_CursorFEside* sfe_2,
                                               PDE_DomainAndFields const* dom,
                                               std::ostream& os )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: iterate_over_sides" ) ;

   PDE_SetOfDiscreteFields* sdf = dom->set_of_discrete_fields() ;
   bool ok_go = true ;
   bool ok_perio = true ;

   os << "***************************************" << endl ;
   os << sfe_1->type_name() <<  endl ;
   os << "   " << sfe_1->nb_meshes() << " sides" << endl ;
   os << "***************************************" << endl ;
   for( sfe_1->start() ; sfe_1->is_valid() ; sfe_1->go_next() )
   {
      os << "---------------------------------------" << endl ;
      os << "current mesh" << endl ;
      sfe_1->print_current_mesh( os, 3 ) ;

      sfe_2->go_i_th( sfe_1->mesh_id() ) ;
      ok_go &= ( sfe_2->is_valid() ) ;
      ok_go &= ( sfe_2->mesh_id() == sfe_1->mesh_id() ) ;

      os << "Local discretization of fields on adjacent cells" << endl ;
      for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
      {
         PDE_DiscreteField const* ff = sdf->item() ;
         if( !dom->is_defined_on_boundary( ff->name() ) )
         {
            print_adjacent_cell_nodes( 0, sfe_1, sfe_2, ff, os, ok_go ) ;
            print_adjacent_cell_nodes( 1, sfe_1, sfe_2, ff, os, ok_go ) ;
         }
      }
      
      if( sfe_1->is_periodic() ) check_periodicity( sfe_1, ok_perio ) ;
   }
   os << "---------------------------------------" << endl ;

   notify_one_test_result( "PDE_CursorFEside go_i_th/go_next", ok_go ) ;
   notify_one_test_result( "PDE_CursorFEside periodicity", ok_perio ) ;
}

//----------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: check_periodicity( PDE_CursorFEside* sfe,
                                              bool& ok_perio )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: check_periodicity" ) ;

   size_t nb_dims = sfe->nb_space_dimensions() ;
   GE_Point* pt = GE_Point::create( 0, nb_dims ) ;

   GE_Mpolyhedron const* poly_0 = sfe->polyhedron( 0 ) ;
   GE_Mpolyhedron const* poly_1 = sfe->polyhedron( 1 ) ;
   GE_SetOfPoints* verts_0 = GE_SetOfPoints::create( 0, nb_dims ) ;
   GE_SetOfPoints* verts_1 = GE_SetOfPoints::create( 0, nb_dims ) ;

   size_t nb_verts = poly_0->nb_vertices() ;
   ok_perio &= ( poly_1->nb_vertices() == nb_verts ) ;
   for( size_t i=0 ; i<nb_verts ; ++i )
   {
      verts_0->append( poly_0->vertex( i ) ) ;
      verts_1->append( poly_1->vertex( i ) ) ;
   }

   GE_Transform const* tr_01 = sfe->periodic_transform( 0 ) ;
   GE_Transform const* tr_10 = sfe->periodic_transform( 1 ) ;
   ok_perio &= ( tr_01->inverse() == tr_10 ) ;
   ok_perio &= ( tr_10->inverse() == tr_01 ) ;

   boolVector found_0( nb_verts ) ; found_0.set( false ) ;
   boolVector found_1( nb_verts ) ; found_1.set( false ) ;
   for( size_t i=0 ; i<nb_verts ; ++i )
   {
      pt->set( poly_0->vertex( i ) ) ;
      tr_01->apply( pt ) ;
      if( verts_1->has( pt ) )
      {
         size_t ii = verts_1->index( pt ) ;
         ok_perio &= !found_0( ii ) ;
         found_0( ii ) = true ;
      }
      else
      {
         ok_perio = false ;
      }

      pt->set( poly_1->vertex( i ) ) ;
      tr_10->apply( pt ) ;
      if( verts_0->has( pt ) )
      {
         size_t ii = verts_0->index( pt ) ;
         ok_perio &= !found_1( ii ) ;
         found_1( ii ) = true ;
      }
      else
      {
         ok_perio = false ;
      }
   }

   for( size_t i=0 ; i<nb_verts ; ++i )
   {
      if( !found_0( i ) || !found_1( i ) ) 
      {
         out() << "--- periodicity error ----------------" << endl ;
         out() << "source polyhedron:" << endl ;
         poly_0->print( out(), 3 ) ;
         out() << "target polyhedron:" << endl ;
         poly_1->print( out(), 3 ) ;
         out() << "direct transformation:" << endl ;
         tr_01->print( out(), 3 ) ;
         out() << "inverse transformation:" << endl ;
         tr_10->print( out(), 3 ) ;
      }
      ok_perio &= found_0( i ) ;
      ok_perio &= found_1( i ) ;
   }

   pt->destroy() ;
   verts_0->destroy() ; verts_0 = 0 ;
   verts_1->destroy() ; verts_1 = 0 ;
}

//----------------------------------------------------------------------------
void
PDE_DomainIteration_TEST:: print_adjacent_cell_nodes( 
                                   size_t i_adj,
                                   PDE_CursorFEside* sfe_1,
                                   PDE_CursorFEside* sfe_2,
                                   PDE_DiscreteField const* ff,
                                   std::ostream& os,
                                   bool& ok ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainIteration_TEST:: print_adjacent_cell_nodes" ) ;

   os << "   \"" << ff->name() << "\"  " 
      << sfe_1->adjacent_localFEcell(i_adj)->nb_local_nodes( ff ) 
      << " local nodes in adjacent cell " << i_adj << endl ;

   ok = ok && ( sfe_1->adjacent_localFEcell(i_adj)->nb_local_nodes( ff ) ==
                sfe_2->adjacent_localFEcell(i_adj)->nb_local_nodes( ff ) ) ;

   for( size_t i=0 ; i<sfe_1->adjacent_localFEcell(i_adj)->nb_local_nodes( ff ) ; ++i )
   {
      os << "      " ;
      os << sfe_1->adjacent_localFEcell(i_adj)->global_node( ff, i ) << "   " ;
      sfe_1->adjacent_localFEcell(i_adj)->local_node_location( ff, i )->print( os, 0 ) ;
      os << endl ;

      ok = ok && ( 
           sfe_1->adjacent_localFEcell(i_adj)->local_node_location( ff, i )->is_equal( 
                   sfe_2->adjacent_localFEcell(i_adj)->local_node_location( ff, i ) ) ) ;
   }  
}

//-----------------------------------------------------------------------
void 
PDE_DomainIteration_TEST:: display_error( std::string const& mesg,
                                          std::string const& title_1,
                                          std::string const& title_2,
                                          double xx_1, double xx_2 ) const
//-----------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = out().flags() ;
   out().setf( ios_base::uppercase | ios_base::scientific ) ;
   out() << std::setprecision( 10 ) ;

   string bl( mesg.length(), ' ' ) ;

   out() << mesg
         << std::setw( 20 ) << title_1
         << std::setw( 20 ) << title_2
         << endl ;
   out() << bl
         << std::setw( 20 ) << xx_1
         << std::setw( 20 ) << xx_2
         << endl ;

   out().flags( original_flags ) ;
}

