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

#include <FE_Parameter_TEST.hh>

#include <PEL.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_assertions.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>

#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <iostream>

FE_Parameter_TEST const*
FE_Parameter_TEST::PROTOTYPE = new FE_Parameter_TEST() ;

//----------------------------------------------------------------------------
FE_Parameter_TEST:: FE_Parameter_TEST( void ) 
//----------------------------------------------------------------------------
   : PEL_ObjectTest( "FE_Parameter", "FE_Parameter_TEST" )
   , CTX( 0 )
   , TIME( 0 )
   , COORDS( 0 )
{
   PEL_ContextSimple* ctx_dummy = PEL_ContextSimple::create( this ) ;
   TIME = PEL_Double::create( ctx_dummy, 0. ) ;
   ctx_dummy->extend( PEL_Variable::object( "DS_T" ), TIME ) ;
   COORDS = PEL_DoubleVector::create( ctx_dummy, doubleVector(0) ) ;
   ctx_dummy->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
   CTX = ctx_dummy ;
}

//----------------------------------------------------------------------------
FE_Parameter_TEST:: ~FE_Parameter_TEST( void )
//----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: process_one_test" ) ;

   PEL_ModuleExplorer* dom_exp =
                          exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields const* dom = PDE_DomainAndFields::create( 0, dom_exp ) ;
   dom_exp->destroy() ; dom_exp = 0 ;

   PEL_ModuleExplorer* param_exp =
                           exp->create_subexplorer( 0, "FE_SetOfParameters" ) ;
   FE_SetOfParameters const* params =
                              FE_SetOfParameters::create( 0, dom, param_exp ) ;
   param_exp->destroy() ; param_exp = 0 ;

   PEL_ModuleExplorer* time_exp =
                              exp->create_subexplorer( 0, "FE_TimeIterator" ) ;
   FE_TimeIterator* t_it = FE_TimeIterator::create( 0, time_exp ) ;
   time_exp->destroy() ; time_exp = 0 ;

   PEL_ModuleExplorer* tests_exp = exp->create_subexplorer( 0, "tests" ) ;
   for( tests_exp->start_module_iterator() ;
        tests_exp->is_valid_module() ;
        tests_exp->go_next_module() )
   {
      PEL_ModuleExplorer* t_exp = tests_exp->create_subexplorer( 0 ) ;
      
      out() << "| " <<  exp->name() << "/" << t_exp->name() << std::endl ;

      FE_Parameter* param = 0 ;

      std::string const& p = t_exp->string_data( "parameter_name" ) ;
      bool ok = params->has( p ) ;
      notify_one_test_result( "FE_SetOfParameters::has", ok ) ;
      if( ok )
      {
         param = params->item( p ) ;
         int const nb_comps = t_exp->int_data( "parameter_nb_components" ) ;
         ok = ( param->nb_components() == (size_t) nb_comps ) ;
         notify_one_test_result( "nb_components", ok ) ;
      }

      if( ok )
      {
         std::string const& type = t_exp->string_data( "type") ;
         D_EPS = t_exp->double_data( "dbl_epsilon" ) ;
         D_MIN = t_exp->double_data( "dbl_minimum" ) ;
         for( t_it->start() ; !t_it->is_finished() ; t_it->go_next_time() )
         {
            if( type == "cell_value" )
            {
               cell_test( dom, param, t_it, t_exp ) ;
            }
            else if(  type == "cell_value_at_IPs" )
            {
               cell_at_IPs_test( dom, param, t_it, t_exp ) ;
            }
            else if(  type == "cell_value_at_centers" )
            {
               cell_at_centers_test( dom, param, t_it, t_exp ) ;
            }
            else if( type == "side_value" )
            {
               side_test( dom, param, t_it, t_exp ) ;
            }
            else if(  type == "side_value_at_IPs" )
            {
               side_at_IPs_test( dom, param, t_it, t_exp ) ;
            }
            else if(  type == "side_value_at_centers" )
            {
               side_at_centers_test( dom, param, t_it, t_exp ) ;
            }
            else if( type == "bound_value" )
            {
               bound_test( dom, param, t_it, t_exp ) ;
            }
            else if(  type == "bound_value_at_IPs" )
            {
               bound_at_IPs_test( dom, param, t_it, t_exp ) ;
            }
            else if(  type == "bound_value_at_centers" )
            {
               bound_at_centers_test( dom, param, t_it, t_exp ) ;
            }
            else
            {
               PEL_Error::object()->raise_bad_data_value(
                  tests_exp, "type",
                  "expected values are: \n"
                  "   - \"cell_value\"\n"
                  "   - \"cell_value_at_IPs\"\n"
                  "   - \"cell_value_at_centers\"\n"
                  "   - \"side_value\"\n"
                  "   - \"side_value_at_IPs\"\n"
                  "   - \"side_value_at_centers\"\n"
                  "   - \"bound_value\"\n"
                  "   - \"bound_value_at_IPs\"\n"
                  "   - \"bound_value_at_centers\"" ) ;
            }
         }
      }
      t_exp->destroy() ; t_exp = 0 ;
   }

   dom->destroy() ; dom = 0 ;
   params->destroy() ; params = 0 ;
   t_it->destroy() ; t_it = 0 ;
   tests_exp->destroy() ; tests_exp = 0 ;
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: cell_test( PDE_DomainAndFields const* dom,
                               FE_Parameter* param,
                               FE_TimeIterator* t_it,
                               PEL_ModuleExplorer* t_exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: cell_test" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( param != 0 ) ;
   PEL_CHECK( t_it != 0 && t_it->is_started() ) ;
   PEL_CHECK( t_exp != 0 ) ;
   
   PDE_LocalFEcell* cFE = dom->create_LocalFEcell( 0 ) ;

   if( t_exp->has_entry( "value" ) )
   {
      doubleVector val( param->nb_components() ) ;
      param->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
      bool ok = param->ok_for_cell_calculations( cFE, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_cell_calculations(Val)", ok ) ;
      if( ok )
      {
         for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
         {
            PEL::out() << cFE->mesh_id() << std::endl ;
            PEL_Context const* ct = context( cFE->polyhedron()->center(),
                                             t_it->time() ) ;
            doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
            if( v.size() != param->nb_components() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "value", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               PEL::out() << "calcul de value : " << std::endl ;
               double value = param->cell_value( t_it, cFE, ic ) ;
               PEL::out() << "   value=" << value << std::endl ;
               ok &= test_equality( value, v(ic) ) ;
            }
         }
         notify_one_test_result( "cell_value", ok ) ;
      }
   }
   
   if( t_exp->has_entry( "gradient" ) )
   {
      doubleArray2D grad( param->nb_components(), dom->nb_space_dimensions() ) ;
      param->transfer_cell_calculation_requirements( cFE, FE_Parameter::Grad ) ;
      bool ok = param->ok_for_cell_calculations( cFE, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_cell_calculations(Grad)", ok ) ;
      if( ok )
      {
         for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
         {
            PEL_Context const* ct = context( cFE->polyhedron()->center(),
                                             t_it->time() ) ;
            doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
            if( g.index_bound(0) != param->nb_components() ||
                g.index_bound(1) != dom->nb_space_dimensions() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "gradient", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
               {
                  double gradient = param->cell_gradient( t_it, cFE, d, ic ) ;
                  ok &= test_equality( gradient, g(ic,d) ) ;
               }
            }
         }
         notify_one_test_result( "cell_gradient", ok ) ;
      }
   }

   cFE->destroy() ; cFE = 0 ;
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: cell_at_centers_test( PDE_DomainAndFields const* dom,
                                          FE_Parameter* param,
                                          FE_TimeIterator* t_it,
                                          PEL_ModuleExplorer* t_exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: cell_at_centers_test" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( param != 0 ) ;
   PEL_CHECK( t_it != 0 && t_it->is_started() ) ;
   PEL_CHECK( t_exp != 0 ) ;
   
   PDE_LocalFEcell* cFE = dom->create_LocalFEcell( 0 ) ;

   if( t_exp->has_entry( "value" ) )
   {
      doubleVector val( param->nb_components() ) ;
      param->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
      bool ok = param->ok_for_cell_calculations( cFE, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_cell_calculations(Val)", ok ) ;
      if( ok )
      {
         for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
         {
            cFE->set_calculation_point( cFE->polyhedron()->center() ) ;
            PEL_Context const* ct = context( cFE->calculation_point(),
                                             t_it->time() ) ;
            doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
            if( v.size() != param->nb_components() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "value", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               double value = param->cell_value_at_pt( t_it, cFE, ic ) ;
               ok &= test_equality( value, v(ic) ) ;
            }
         }
         notify_one_test_result( "cell_value_at_pt", ok ) ;
      }
   }
   
   if( t_exp->has_entry( "gradient" ) )
   {
      doubleArray2D grad( param->nb_components(), dom->nb_space_dimensions() ) ;
      param->transfer_cell_calculation_requirements( cFE, FE_Parameter::Grad ) ;
      bool ok = param->ok_for_cell_calculations( cFE, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_cell_calculations(Grad)", ok ) ;
      if( ok )
      {
         for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
         {
            cFE->set_calculation_point( cFE->polyhedron()->center() ) ;
            PEL_Context const* ct = context( cFE->calculation_point(),
                                             t_it->time() ) ;
            doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
            if( g.index_bound(0) != param->nb_components() ||
                g.index_bound(1) != dom->nb_space_dimensions() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "gradient", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
               {
                  double gradient = param->cell_gradient_at_pt( t_it, cFE,
                                                                d, ic ) ;
                  ok &= test_equality( gradient, g(ic,d) ) ;
               }
            }
         }
         notify_one_test_result( "cell_gradient_at_pt", ok ) ;
      }
   }

   cFE->destroy() ; cFE = 0 ;
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: cell_at_IPs_test( PDE_DomainAndFields const* dom,
                                      FE_Parameter* param,
                                      FE_TimeIterator* t_it,
                                      PEL_ModuleExplorer* t_exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: cell_at_IPs_test" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( param != 0 ) ;
   PEL_CHECK( t_it != 0 && t_it->is_started() ) ;
   PEL_CHECK( t_exp != 0 ) ;
   
   double xxx = PEL::bad_double() ;
   
   PDE_LocalFEcell*  cFE = dom->create_LocalFEcell( 0 ) ;
   PDE_CursorFEside* sFE = dom->create_CursorFEside( 0 ) ;
   PDE_LocalFEcell const* cfe0 = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* cfe1 = sFE->adjacent_localFEcell( 1 ) ;
   
   GE_QRprovider const* qr =
      GE_QRprovider::object( t_exp->string_data( "QRprovider_name" ) ) ;

   if( t_exp->has_entry( "value" ) )
   {
      doubleVector val( param->nb_components() ) ;
      param->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
      bool ok = param->ok_for_cell_calculations( cFE, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_cell_calculations(Val)", ok ) ;
      if( !ok ) return ; 
      
      for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
      {
         cFE->start_IP_iterator( qr ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            PEL_Context const* ct = context( cFE->coordinates_of_IP(),
                                             t_it->time() ) ;
            doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
            if( v.size() != param->nb_components() )
            {
               PEL_Error::object()->raise_bad_data_value(
                                                         t_exp, "value", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               xxx = param->cell_value_at_IP( t_it, cFE, ic ) ;
               ok &= test_equality( xxx, v(ic) ) ;
            }
         }
      }
      notify_one_test_result( "cell_value_at_IP", ok ) ;
      
      param->transfer_side_calculation_requirements( sFE, FE_Parameter::Val ) ;
      ok = param->ok_for_side_calculations( sFE, FE_Parameter::Val ) ;
      ok &= param->ok_for_cell_calculations( cfe0, FE_Parameter::Val ) ;
      ok &= param->ok_for_cell_calculations( cfe1, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_side_calculations(Val)", ok ) ;
      if( !ok ) return ;
      
      for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
      {
         sFE->start_IP_iterator( qr ) ;
         for( ; sFE->valid_IP() ; sFE->go_next_IP() )
         {
            PEL_Context const* ct = context( sFE->coordinates_of_IP(),
                                             t_it->time() ) ;
            doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               xxx = param->cell_value_at_IP( t_it, cfe0, ic ) ;
               ok &= test_equality( xxx, v(ic) ) ;
               xxx = param->cell_value_at_IP( t_it, cfe1, ic ) ;
               ok &= test_equality( xxx, v(ic) ) ;
            }
         }
      }      
      notify_one_test_result( "cell_value_at_IP (from sides)", ok ) ;
   }
            
   if( t_exp->has_entry( "gradient" ) )
   {
      doubleArray2D grad( param->nb_components(), dom->nb_space_dimensions() ) ;
      param->transfer_cell_calculation_requirements( cFE, FE_Parameter::Grad ) ;
      bool ok = param->ok_for_cell_calculations( cFE, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_cell_calculations(Grad)", ok ) ;
      if( !ok ) return ;

      for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
      {
         cFE->start_IP_iterator( qr ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            PEL_Context const* ct = context( cFE->coordinates_of_IP(),
                                             t_it->time() ) ;
            doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
               {
                  xxx = param->cell_gradient_at_IP( t_it, cFE, d, ic ) ;
                  ok &= test_equality( xxx, g(ic,d) ) ;
               }
            }
         }
      }
      notify_one_test_result( "cell_gradient_at_IP", ok ) ;
      
      param->transfer_side_calculation_requirements( sFE, FE_Parameter::Grad ) ;
      ok = param->ok_for_side_calculations( sFE, FE_Parameter::Grad ) ;
      ok &= param->ok_for_cell_calculations( cfe0, FE_Parameter::Grad ) ;
      ok &= param->ok_for_cell_calculations( cfe1, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_side_calculations(Val)", ok ) ;
      if( !ok ) return ;
      
      for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
      {
         sFE->start_IP_iterator( qr ) ;
         for( ; sFE->valid_IP() ; sFE->go_next_IP() )
         {
            PEL_Context const* ct = context( sFE->coordinates_of_IP(),
                                             t_it->time() ) ;
            doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
               {
                  xxx = param->cell_gradient_at_IP( t_it, cfe0, d, ic ) ;
                  ok &= test_equality( xxx, g(ic,d) ) ;
                  xxx = param->cell_gradient_at_IP( t_it, cfe1, d, ic ) ;
                  ok &= test_equality( xxx, g(ic,d) ) ;
               }
            }
         }
      }
      notify_one_test_result( "cell_gradient_at_IP (from sides)", ok ) ;
   }

   cFE->destroy() ; cFE = 0 ;
   sFE->destroy() ; sFE = 0 ;
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: bound_test( PDE_DomainAndFields const* dom,
                               FE_Parameter* param,
                               FE_TimeIterator* t_it,
                               PEL_ModuleExplorer* t_exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: bound_test" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( param != 0 ) ;
   PEL_CHECK( t_it != 0 && t_it->is_started() ) ;
   PEL_CHECK( t_exp != 0 ) ;
   
   PDE_LocalFEbound* bFE = dom->create_LocalFEbound( 0 ) ;

   if( t_exp->has_entry( "value" ) )
   {
      doubleVector val( param->nb_components() ) ;
      param->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
      bool ok = param->ok_for_bound_calculations( bFE, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_bound_calculations(Val)", ok ) ;
      if( ok )
      {
         for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
         {
            PEL_Context const* ct = context( bFE->polyhedron()->center(),
                                             t_it->time() ) ;
            doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
            if( v.size() != param->nb_components() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "value", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               double value = param->bound_value( t_it, bFE, ic ) ;
               ok &= test_equality( value, v(ic) ) ;
            }
         }
         notify_one_test_result( "bound_value", ok ) ;
      }
   }
   
   if( t_exp->has_entry( "gradient" ) )
   {
      doubleArray2D grad( param->nb_components(), dom->nb_space_dimensions() ) ;
      param->transfer_bound_calculation_requirements( bFE, FE_Parameter::Grad ) ;
      bool ok = param->ok_for_bound_calculations( bFE, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_bound_calculations(Grad)", ok ) ;
      if( ok )
      {
         for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
         {
            PEL_Context const* ct = context( bFE->polyhedron()->center(),
                                             t_it->time() ) ;
            doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
            if( g.index_bound(0) != param->nb_components() ||
                g.index_bound(1) != dom->nb_space_dimensions() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "gradient", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
               {
                  double gradient = param->bound_gradient( t_it, bFE, d, ic ) ;
                  ok &= test_equality( gradient, g(ic,d) ) ;
               }
            }
         }
         notify_one_test_result( "bound_gradient", ok ) ;
      }
   }

   bFE->destroy() ; bFE = 0 ;
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: bound_at_centers_test( PDE_DomainAndFields const* dom,
                                          FE_Parameter* param,
                                          FE_TimeIterator* t_it,
                                          PEL_ModuleExplorer* t_exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: bound_at_centers_test" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( param != 0 ) ;
   PEL_CHECK( t_it != 0 && t_it->is_started() ) ;
   PEL_CHECK( t_exp != 0 ) ;
   
   PDE_LocalFEbound* bFE = dom->create_LocalFEbound( 0 ) ;

   if( t_exp->has_entry( "value" ) )
   {
      doubleVector val( param->nb_components() ) ;
      param->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
      bool ok = param->ok_for_bound_calculations( bFE, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_bound_calculations(Val)", ok ) ;
      if( ok )
      {
         for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
         {
            bFE->set_calculation_point( bFE->polyhedron()->center() ) ;
            PEL_Context const* ct = context( bFE->calculation_point(),
                                             t_it->time() ) ;
            doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
            if( v.size() != param->nb_components() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "value", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               double value = param->bound_value_at_pt( t_it, bFE, ic ) ;
               ok &= test_equality( value, v(ic) ) ;
            }
         }
         notify_one_test_result( "bound_value_at_pt", ok ) ;
      }
   }
   
   if( t_exp->has_entry( "gradient" ) )
   {
      doubleArray2D grad( param->nb_components(), dom->nb_space_dimensions() ) ;
      param->transfer_bound_calculation_requirements( bFE, FE_Parameter::Grad ) ;
      bool ok = param->ok_for_bound_calculations( bFE, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_bound_calculations(Grad)", ok ) ;
      if( ok )
      {
         for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
         {
            bFE->set_calculation_point( bFE->polyhedron()->center() ) ;
            PEL_Context const* ct = context( bFE->calculation_point(),
                                             t_it->time() ) ;
            doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
            if( g.index_bound(0) != param->nb_components() ||
                g.index_bound(1) != dom->nb_space_dimensions() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "gradient", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
               {
                  double gradient = param->bound_gradient_at_pt( t_it, bFE,
                                                                d, ic ) ;
                  ok &= test_equality( gradient, g(ic,d) ) ;
               }
            }
         }
         notify_one_test_result( "bound_gradient_at_pt", ok ) ;
      }
   }

   bFE->destroy() ; bFE = 0 ;
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: bound_at_IPs_test( PDE_DomainAndFields const* dom,
                                      FE_Parameter* param,
                                      FE_TimeIterator* t_it,
                                      PEL_ModuleExplorer* t_exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: bound_at_IPs_test" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( param != 0 ) ;
   PEL_CHECK( t_it != 0 && t_it->is_started() ) ;
   PEL_CHECK( t_exp != 0 ) ;
   
   PDE_LocalFEbound* bFE = dom->create_LocalFEbound( 0 ) ;
   
   GE_QRprovider const* qr =
      GE_QRprovider::object( t_exp->string_data( "QRprovider_name" ) ) ;

   if( t_exp->has_entry( "value" ) )
   {
      doubleVector val( param->nb_components() ) ;
      param->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
      bool ok = param->ok_for_bound_calculations( bFE, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_bound_calculations(Val)", ok ) ;
      if( ok )
      {
         for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
         {
            bFE->start_IP_iterator( qr ) ;
            for( ; bFE->valid_IP() ; bFE->go_next_IP() )
            {
               PEL_Context const* ct = context( bFE->coordinates_of_IP(),
                                                t_it->time() ) ;
               doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
               if( v.size() != param->nb_components() )
               {
                  PEL_Error::object()->raise_bad_data_value(
                     t_exp, "value", "bad number of components" ) ;
               }
               for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
               {
                  double value = param->bound_value_at_IP( t_it, bFE, ic ) ;
                  ok &= test_equality( value, v(ic) ) ;
               }
            }
         }
         notify_one_test_result( "bound_value_at_IP", ok ) ;
      }
   }
   
   if( t_exp->has_entry( "gradient" ) )
   {
      doubleArray2D grad( param->nb_components(), dom->nb_space_dimensions() ) ;
      param->transfer_bound_calculation_requirements( bFE, FE_Parameter::Grad ) ;
      bool ok = param->ok_for_bound_calculations( bFE, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_bound_calculations(Grad)", ok ) ;
      if( ok )
      {
         for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
         {
            bFE->start_IP_iterator( qr ) ;
            for( ; bFE->valid_IP() ; bFE->go_next_IP() )
            {
               PEL_Context const* ct = context( bFE->coordinates_of_IP(),
                                                t_it->time() ) ;
               doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
               if( g.index_bound(0) != param->nb_components() ||
                   g.index_bound(1) != dom->nb_space_dimensions() )
               {
                  PEL_Error::object()->raise_bad_data_value(
                     t_exp, "gradient", "bad number of components" ) ;
               }
               for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
               {
                  for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
                  {
                     double gradient = param->bound_gradient_at_IP( t_it, bFE,
                                                                   d, ic ) ;
                     ok &= test_equality( gradient, g(ic,d) ) ;
                  }
               }
            }
         }
         notify_one_test_result( "bound_gradient_at_IP", ok ) ;
      }
   }

   bFE->destroy() ; bFE = 0 ;
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: side_test( PDE_DomainAndFields const* dom,
                               FE_Parameter* param,
                               FE_TimeIterator* t_it,
                               PEL_ModuleExplorer* t_exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: side_test" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( param != 0 ) ;
   PEL_CHECK( t_it != 0 && t_it->is_started() ) ;
   PEL_CHECK( t_exp != 0 ) ;
   
   PDE_CursorFEside* sFE = dom->create_CursorFEside( 0 ) ;

   if( t_exp->has_entry( "value" ) )
   {
      doubleVector val( param->nb_components() ) ;
      param->transfer_side_calculation_requirements( sFE, FE_Parameter::Val ) ;
      bool ok = param->ok_for_side_calculations( sFE, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_side_calculations(Val)", ok ) ;
      if( ok )
      {
         for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
         {
            PEL_Context const* ct = context( sFE->polyhedron()->center(),
                                             t_it->time() ) ;
            doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
            if( v.size() != param->nb_components() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "value", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               double value = param->side_value( t_it, sFE, ic ) ;
               ok &= test_equality( value, v(ic) ) ;
            }
         }
         notify_one_test_result( "side_value", ok ) ;
      }
   }
   
   if( t_exp->has_entry( "gradient" ) )
   {
      doubleArray2D grad( param->nb_components(), dom->nb_space_dimensions() ) ;
      param->transfer_side_calculation_requirements( sFE, FE_Parameter::Grad ) ;
      bool ok = param->ok_for_side_calculations( sFE, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_side_calculations(Grad)", ok ) ;
      if( ok )
      {
         for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
         {
            PEL_Context const* ct = context( sFE->polyhedron()->center(),
                                             t_it->time() ) ;
            doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
            if( g.index_bound(0) != param->nb_components() ||
                g.index_bound(1) != dom->nb_space_dimensions() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "gradient", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
               {
                  double gradient = param->side_gradient( t_it, sFE, d, ic ) ;
                  ok &= test_equality( gradient, g(ic,d) ) ;
               }
            }
         }
         notify_one_test_result( "side_gradient", ok ) ;
      }
   }

   sFE->destroy() ; sFE = 0 ;
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: side_at_centers_test( PDE_DomainAndFields const* dom,
                                          FE_Parameter* param,
                                          FE_TimeIterator* t_it,
                                          PEL_ModuleExplorer* t_exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: side_at_centers_test" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( param != 0 ) ;
   PEL_CHECK( t_it != 0 && t_it->is_started() ) ;
   PEL_CHECK( t_exp != 0 ) ;
   
   PDE_CursorFEside* sFE = dom->create_CursorFEside( 0 ) ;

   if( t_exp->has_entry( "value" ) )
   {
      doubleVector val( param->nb_components() ) ;
      param->transfer_side_calculation_requirements( sFE, FE_Parameter::Val ) ;
      bool ok = param->ok_for_side_calculations( sFE, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_side_calculations(Val)", ok ) ;
      if( ok )
      {
         for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
         {
            sFE->set_calculation_point( sFE->polyhedron()->center() ) ;
            PEL_Context const* ct = context( sFE->calculation_point(),
                                             t_it->time() ) ;
            doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
            if( v.size() != param->nb_components() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "value", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               double value = param->side_value_at_pt( t_it, sFE, ic ) ;
               ok &= test_equality( value, v(ic) ) ;
            }
         }
         notify_one_test_result( "side_value_at_pt", ok ) ;
      }
   }
   
   if( t_exp->has_entry( "gradient" ) )
   {
      doubleArray2D grad( param->nb_components(), dom->nb_space_dimensions() ) ;
      param->transfer_side_calculation_requirements( sFE, FE_Parameter::Grad ) ;
      bool ok = param->ok_for_side_calculations( sFE, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_side_calculations(Grad)", ok ) ;
      if( ok )
      {
         for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
         {
            sFE->set_calculation_point( sFE->polyhedron()->center() ) ;
            PEL_Context const* ct = context( sFE->calculation_point(),
                                             t_it->time() ) ;
            doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
            if( g.index_bound(0) != param->nb_components() ||
                g.index_bound(1) != dom->nb_space_dimensions() )
            {
               PEL_Error::object()->raise_bad_data_value(
                  t_exp, "gradient", "bad number of components" ) ;
            }
            for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
            {
               for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
               {
                  double gradient = param->side_gradient_at_pt( t_it, sFE,
                                                                d, ic ) ;
                  ok &= test_equality( gradient, g(ic,d) ) ;
               }
            }
         }
         notify_one_test_result( "side_gradient_at_pt", ok ) ;
      }
   }

   sFE->destroy() ; sFE = 0 ;
}

//-----------------------------------------------------------------------------
void
FE_Parameter_TEST:: side_at_IPs_test( PDE_DomainAndFields const* dom,
                                      FE_Parameter* param,
                                      FE_TimeIterator* t_it,
                                      PEL_ModuleExplorer* t_exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: side_at_IPs_test" ) ;
   PEL_CHECK( dom != 0 ) ;
   PEL_CHECK( param != 0 ) ;
   PEL_CHECK( t_it != 0 && t_it->is_started() ) ;
   PEL_CHECK( t_exp != 0 ) ;
   
   PDE_CursorFEside* sFE = dom->create_CursorFEside( 0 ) ;
   
   GE_QRprovider const* qr =
      GE_QRprovider::object( t_exp->string_data( "QRprovider_name" ) ) ;

   if( t_exp->has_entry( "value" ) )
   {
      doubleVector val( param->nb_components() ) ;
      param->transfer_side_calculation_requirements( sFE, FE_Parameter::Val ) ;
      bool ok = param->ok_for_side_calculations( sFE, FE_Parameter::Val ) ;
      notify_one_test_result( "ok_for_side_calculations(Val)", ok ) ;
      if( ok )
      {
         for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
         {
            sFE->start_IP_iterator( qr ) ;
            for( ; sFE->valid_IP() ; sFE->go_next_IP() )
            {
               PEL_Context const* ct = context( sFE->coordinates_of_IP(),
                                                t_it->time() ) ;
               doubleVector const& v = t_exp->doubleVector_data( "value", ct ) ;
               if( v.size() != param->nb_components() )
               {
                  PEL_Error::object()->raise_bad_data_value(
                     t_exp, "value", "bad number of components" ) ;
               }
               for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
               {
                  double value = param->side_value_at_IP( t_it, sFE, ic ) ;
                  ok &= test_equality( value, v(ic) ) ;
               }
            }
         }
         notify_one_test_result( "side_value_at_IP", ok ) ;
      }
   }
   
   if( t_exp->has_entry( "gradient" ) )
   {
      doubleArray2D grad( param->nb_components(), dom->nb_space_dimensions() ) ;
      param->transfer_side_calculation_requirements( sFE, FE_Parameter::Grad ) ;
      bool ok = param->ok_for_side_calculations( sFE, FE_Parameter::Grad ) ;
      notify_one_test_result( "ok_for_side_calculations(Grad)", ok ) ;
      if( ok )
      {
         for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
         {
            sFE->start_IP_iterator( qr ) ;
            for( ; sFE->valid_IP() ; sFE->go_next_IP() )
            {
               PEL_Context const* ct = context( sFE->coordinates_of_IP(),
                                                t_it->time() ) ;
               doubleArray2D const& g = t_exp->doubleArray2D_data( "gradient", ct ) ;
               if( g.index_bound(0) != param->nb_components() ||
                   g.index_bound(1) != dom->nb_space_dimensions() )
               {
                  PEL_Error::object()->raise_bad_data_value(
                     t_exp, "gradient", "bad number of components" ) ;
               }
               for( size_t ic = 0 ; ic<param->nb_components() ; ++ic )
               {
                  for( size_t d=0 ; d<dom->nb_space_dimensions() ; ++d )
                  {
                     double gradient = param->side_gradient_at_IP( t_it, sFE,
                                                                   d, ic ) ;
                     ok &= test_equality( gradient, g(ic,d) ) ;
                  }
               }
            }
         }
         notify_one_test_result( "side_gradient_at_IP", ok ) ;
      }
   }

   sFE->destroy() ; sFE = 0 ;
}

//-----------------------------------------------------------------------------
PEL_Context const*
FE_Parameter_TEST:: context( GE_Point const* pt, double time )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: context" ) ;
   PEL_CHECK( pt != 0 ) ;

   TIME->set( time ) ;
   COORDS->set( pt->coordinate_vector() ) ;
   return( CTX ) ;
}

//-----------------------------------------------------------------------------
bool
FE_Parameter_TEST:: test_equality( double value, double expected_value )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter_TEST:: test_equality" ) ;

   bool result = PEL::double_equality( value, expected_value,
                                       D_EPS, D_MIN ) ;
   if( !result )
   {
      out() << "   expected : " << expected_value << std::endl ;
      out() << "   computed : " << value << std::endl ;
      out() << std::endl ;
   }
   return( result ) ;
}
