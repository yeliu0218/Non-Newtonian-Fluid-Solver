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

#include <PDE_LocalFEboundAdjacencies_TEST.hh>

#include <PEL.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>
#include <doubleVector.hh>
#include <stringVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_GridMover.hh>
#include <PDE_InterfaceAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfDomains.hh>

#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>

using std::endl ;
using std::ios_base ;
using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

//---------------------------------------------------------------------------
PDE_LocalFEboundAdjacencies_TEST*
PDE_LocalFEboundAdjacencies_TEST:: registered_test = 
                                   new PDE_LocalFEboundAdjacencies_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_LocalFEboundAdjacencies_TEST:: PDE_LocalFEboundAdjacencies_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_LocalFEbound", "PDE_LocalFEboundAdjacencies_TEST" )
{
   CONTEXT = PEL_ContextSimple::create( this ) ;
   COORDS  = PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) ;
   CONTEXT->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;   
}

//---------------------------------------------------------------------------
PDE_LocalFEboundAdjacencies_TEST:: ~PDE_LocalFEboundAdjacencies_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_LocalFEboundAdjacencies_TEST:: process_one_test( 
                                              PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEboundAdjacencies_TEST:: process_one_test" ) ;

   D_EPS = exp->double_data( "dbl_epsilon" ) ;
   D_MIN = exp->double_data( "dbl_minimum" ) ;
   QRP = GE_QRprovider::object( 
                        exp->string_data( "quadrature_rule_provider" ) ) ;

   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "PDE_SetOfDomains" ) ;
   PDE_SetOfDomains* sdoms = PDE_SetOfDomains::create( 0, se ) ;
   se->destroy() ; se = 0 ;

   std::queue< double, std::list<double> > valeurs ;

   do_all_checks( sdoms, exp, valeurs, true ) ;

   if( exp->has_module( "grid_displacement" ) )
   {
      PEL_ModuleExplorer* ssee = 
                           exp->create_subexplorer( 0, "grid_displacement" ) ;
      PDE_DomainAndFields const* dom = 
                           sdoms->domain( ssee->string_data( "domain" ) ) ;
      PDE_GridMover* gm = PDE_GridMover::create( 0, dom ) ;

      out() << "*** deplacement du maillage" << endl ;
      gm->move_grid( dom->set_of_discrete_fields()->item( 
                              ssee->string_data( "field" ) ),
                     0, 1.0 ) ;

      do_all_checks( sdoms, exp, valeurs, false ) ;

      gm->destroy() ; gm = 0 ;
      ssee->destroy() ; ssee = 0 ;
   }

   sdoms->destroy() ;
}

//---------------------------------------------------------------------------
void
PDE_LocalFEboundAdjacencies_TEST:: do_all_checks( 
                            PDE_SetOfDomains const* sdoms,
                            PEL_ModuleExplorer const* exp,
                            std::queue< double, std::list<double> >& valeurs,
                            bool learn )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEboundAdjacencies_TEST:: do_all_checks" ) ;

   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "list_of_checks" ) ;
   se->start_module_iterator() ;
   for( ; se->is_valid_module() ; se->go_next_module() )
   {
      PEL_ModuleExplorer* ee = se->create_subexplorer( 0 ) ;

      PEL_DataWithContext* val = 0 ;
      PEL_DataWithContext* jac = 0 ;
      PEL_DataWithContext* hes = 0 ;

      PDE_DomainAndFields* dom = 
         sdoms->domain( ee->string_data( "domain_for_iteration" ) ) ;
      PDE_DomainAndFields* dom_ext =
         sdoms->domain( ee->string_data( "adjacent_domain" ) ) ;

      std::string const& fn = ee->string_data( "field_of_adjacent_domain" ) ;

      PEL_ASSERT( !dom->set_of_discrete_fields()->has( fn ) ) ;

      PDE_SetOfDiscreteFields* sdf_ext = dom_ext->set_of_discrete_fields() ;
      PDE_DiscreteField const* ff = sdf_ext->item( fn ) ;

      GE_Color const* interf_color = GE_Color::object( 
         ee->string_data( "color_of_bounds_with_adjacency" ) ) ;

      PDE_LocalFEbound* fe = dom->create_LocalFEbound( 0 ) ;
      fe->require_field_calculation( ff, PDE_LocalFE::node ) ;
      if( ee->has_entry( "value" ) )
      {
         val = ee->abstract_data( ee, "value", CONTEXT ) ;
         if( val->data_type()!=PEL_Data::DoubleVector )
         {
            PEL_Error::object()->raise_bad_data_type(
               ee, "value", PEL_Data::DoubleVector  ) ;
         }
         if( !val->value_can_be_evaluated(0) )
         {
            PEL_Error::object()->raise_not_evaluable(
               ee, "value", val->undefined_variables(0) ) ;
         }
         fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
         if( ee->has_entry( "jacobian" ) )
         {
            jac = ee->abstract_data( ee, "jacobian", CONTEXT ) ;
            if( jac->data_type()!=PEL_Data::DoubleArray2D )
            {
               PEL_Error::object()->raise_bad_data_type(
                  ee, "jacobian", PEL_Data::DoubleArray2D  ) ;
            }
            if( !jac->value_can_be_evaluated(0) )
            {
               PEL_Error::object()->raise_not_evaluable(
                  ee, "jacobian", jac->undefined_variables(0) ) ;
            }
            fe->require_field_calculation( ff, PDE_LocalFE::dN ) ;
            if( ee->has_entry( "hessian" ) )
            {
               hes = ee->abstract_data( ee, "hessian", CONTEXT ) ;
               if( hes->data_type()!=PEL_Data::DoubleArray3D )
               {
                  PEL_Error::object()->raise_bad_data_type(
                     ee, "hessian", PEL_Data::DoubleArray3D  ) ;
               }
               if( !hes->value_can_be_evaluated(0) )
               {
                  PEL_Error::object()->raise_not_evaluable(
                     ee, "hessian", hes->undefined_variables(0) ) ;
               }
               fe->require_field_calculation( ff, PDE_LocalFE::d2N ) ;
            }
         }
      }

      iterate_and_check_CP( fe, interf_color, ff, val, jac, hes, 
                            valeurs, learn ) ;

      iterate_and_check_IP( fe, interf_color, ff, val, jac, hes,
                            valeurs, learn ) ;

      fe->destroy() ;
      ee->destroy() ;
   }
   se->destroy() ;
}

//---------------------------------------------------------------------------
void
PDE_LocalFEboundAdjacencies_TEST:: iterate_and_check_IP(
                            PDE_LocalFEbound* fe,
                            GE_Color const* interf_color,
                            PDE_DiscreteField const* ff,
                            PEL_DataWithContext const* val,
                            PEL_DataWithContext const* jac,
                            PEL_DataWithContext const* hes,
                            std::queue< double, std::list<double> >& valeurs,
                            bool learn )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEboundAdjacencies_TEST:: iterate_and_check_IP" ) ;

   bool ok = true ;
   bool ok_val = true ;
   bool ok_grad = true ;
   bool ok_hes = true ;
   bool eq ;

   double theo, xx ;

   size_t nb_dims = fe->nb_space_dimensions() ;

   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      GE_Color const* color = fe->color() ;
      if( color->is_matching( interf_color ) )
      {
         ok = ok && ( fe->nb_local_nodes( ff ) != 0 ) ;
         if( ok && val!=0 )
         {
            fe->start_IP_iterator( QRP ) ;
            for( ; fe->valid_IP() ; fe->go_next_IP() )
            {
               fill_context_with_coordinates( fe->coordinates_of_IP() ) ;
               for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
               {
                  get_theo( val->to_double_vector()( ic ), theo,
                            valeurs, learn ) ;
                  xx = fe->value_at_IP( ff, 0, ic ) ;
                  eq = PEL::double_equality( theo, xx, D_EPS, D_MIN ) ;
                  ok_val = ok_val && eq ;
                  if( !eq ) display_error( ff->name(),
                                           "data deck", "value_at_IP",
                                           theo, xx ) ;
                  for( size_t a=0 ; (jac != 0) && a<nb_dims ; ++a )
                  {
                     get_theo( jac->to_double_array2D()( ic, a ),
                               theo, valeurs, learn ) ;
                     xx = fe->gradient_at_IP( ff, 0, a, ic ) ;
                     eq = PEL::double_equality( theo, xx, D_EPS, D_MIN ) ;
                     ok_grad = ok_grad && eq ;
                     if( !eq ) display_error( ff->name(),  
                                              "data deck", "gradient_at_IP",
                                              theo, xx ) ;

                     for( size_t b=0 ; (hes!=0) && b<nb_dims ; ++b )
                     {
                        get_theo( hes->to_double_array3D()( ic, a, b ),
                                  theo, valeurs, learn ) ;
                        xx = fe->hessian_at_IP( ff, 0, a, b, ic ) ;
                        eq = PEL::double_equality( theo, xx, D_EPS, D_MIN ) ;
                        ok_hes = ok_hes && eq ;
                        if( !eq ) display_error( ff->name(), 
                                                 "data deck", "hessian_at_IP",
                                                 theo, xx ) ;
                     }
                  }
               }
            }
         }
      }
      else
      {
         ok = ok && ( fe->nb_local_nodes( ff ) == 0 ) ;
      }
   }
   notify_one_test_result( ff->name()+" localization", ok ) ;
   if( ok )
   {
      if( val ) notify_one_test_result( ff->name() + " value_at_IP", 
                                        ok_val ) ; 
      if( jac ) notify_one_test_result( ff->name() + " gradient_at_IP", 
                                        ok_grad ) ;
      if( hes ) notify_one_test_result( ff->name() + " hessian_at_IP", 
                                        ok_hes ) ;
   }
}

//---------------------------------------------------------------------------
void
PDE_LocalFEboundAdjacencies_TEST:: iterate_and_check_CP(
                            PDE_LocalFEbound* fe,
                            GE_Color const* interf_color,
                            PDE_DiscreteField const* ff,
                            PEL_DataWithContext const* val,
                            PEL_DataWithContext const* jac,
                            PEL_DataWithContext const* hes,
                            std::queue< double, std::list<double> >& valeurs,
                            bool learn )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFEboundAdjacencies_TEST:: iterate_and_check_CP" ) ;

   bool ok = true ;
   bool ok_val = true ;
   bool ok_grad = true ;
   bool ok_hes = true ;
   bool eq ;

   double theo, xx ;

   size_t nb_dims = fe->nb_space_dimensions() ;

   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      GE_Color const* color = fe->color() ;
      if( color->is_matching( interf_color ) )
      {
         ok = ok && ( fe->nb_local_nodes( ff ) != 0 ) ;
         if( ok && val!=0 )
         {
            GE_Point const* pt = fe->polyhedron()->center() ;
            fe->set_calculation_point( pt ) ;
            fill_context_with_coordinates( pt ) ;
            for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
            {
               get_theo( val->to_double_vector()( ic ), theo,
                         valeurs, learn ) ; ;
               xx = fe->value_at_pt( ff, 0, ic ) ;
               eq = PEL::double_equality( theo, xx, D_EPS, D_MIN ) ;
               ok_val = ok_val && eq ;
               if( !eq ) display_error( ff->name(),
                                        "data deck", "value_at_pt",
                                        theo, xx ) ;
               for( size_t a=0 ; (jac != 0) && a<nb_dims ; ++a )
               {
                  get_theo( jac->to_double_array2D( CONTEXT)( ic, a ),
                            theo, valeurs, learn ) ;
                  xx = fe->gradient_at_pt( ff, 0, a, ic ) ;
                  eq = PEL::double_equality( theo, xx, D_EPS, D_MIN ) ;
                  ok_grad = ok_grad && eq ;
                  if( !eq ) display_error( ff->name(),  
                                           "data deck", "gradient_at_pt",
                                           theo, xx ) ;

                  for( size_t b=0 ; (hes!=0) && b<nb_dims ; ++b )
                  {
                     get_theo( hes->to_double_array3D()( ic, a, b ),
                               theo, valeurs, learn ) ;
                     xx = fe->hessian_at_pt( ff, 0, a, b, ic ) ;
                     eq = PEL::double_equality( theo, xx, D_EPS, D_MIN ) ;
                     ok_hes = ok_hes && eq ;
                     if( !eq ) display_error( ff->name(), 
                                              "data deck", "hessian_at_pt",
                                              theo, xx ) ;
                  }
               }
            }
         }
      }
      else
      {
         ok = ok && ( fe->nb_local_nodes( ff ) == 0 ) ;
      }
   }
   notify_one_test_result( ff->name() + " localization", ok ) ;
   if( ok )
   {
      if( val ) notify_one_test_result( ff->name() + " value_at_pt", 
                                        ok_val ) ;
      if( jac ) notify_one_test_result( ff->name() + " gradient_at_pt", 
                                        ok_grad ) ;
      if( hes ) notify_one_test_result( ff->name() + " hessian_at_pt", 
                                        ok_hes ) ;
   }
}

//-----------------------------------------------------------------------
void
PDE_LocalFEboundAdjacencies_TEST:: get_theo( 
                            double val, double& theo,
                            std::queue< double, std::list<double> >& valeurs, 
                            bool learn ) const
//-----------------------------------------------------------------------
{
   if( learn )
   {
      theo = val ;
      valeurs.push( theo ) ;
   }
   else
   {
      theo = valeurs.front() ;
      valeurs.pop() ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFEboundAdjacencies_TEST:: fill_context_with_coordinates( 
                                                     GE_Point const* pt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE_TEST:: fill_context_with_coordinates" ) ;

   size_t n = pt->nb_coordinates() ;
   doubleVector dv( n ) ;
   for( size_t i=0 ; i<n ; ++i )
   {
      dv(i) = pt->coordinate( i ) ;
   }
   COORDS->set( dv ) ;
}

//-----------------------------------------------------------------------
void PDE_LocalFEboundAdjacencies_TEST::  display_error( 
                                          std::string const& mesg,
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

