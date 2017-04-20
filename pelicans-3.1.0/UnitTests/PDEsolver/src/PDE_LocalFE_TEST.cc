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

#include <PDE_LocalFE_TEST.hh>

#include <PEL.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <doubleArray2D.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>

#include <PDE_AdapterCHARMS.hh>
#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalFEinterface.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <iostream>

using std::cout ;
using std::endl ;
using std::string ;

//---------------------------------------------------------------------------
PDE_LocalFE_TEST*
PDE_LocalFE_TEST:: registered_test = new PDE_LocalFE_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_LocalFE_TEST:: PDE_LocalFE_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_LocalFE", "PDE_LocalFE_TEST" )
   , MY_EPS( PEL::bad_double() )
   , MY_MIN( PEL::bad_double() )
   , FES( 0 )
   , I_FIELDS( 0 )
   , I_VALUES( 0 )
   , JACOBIANS( 0 )
   , HESSIANS( 0 )
   , BD_FIELDS( 0 )
   , BD_VALUES( 0 )
   , CONTEXT( 0 )
   , COORDS( 0 )
   , CFE( 0 )
   , SFE( 0 )
   , BFE( 0 )
{
}

//---------------------------------------------------------------------------
PDE_LocalFE_TEST:: ~PDE_LocalFE_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_LocalFE_TEST:: initialize_internals( PEL_Object* a_owner )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE_TEST:: initialize_internals" ) ;

   FES       = PEL_Vector::create( a_owner, 0 ) ;
   I_FIELDS  = PEL_Vector::create( a_owner, 0 ) ;
   I_VALUES  = PEL_Vector::create( a_owner, 0 ) ;
   JACOBIANS = PEL_Vector::create( a_owner, 0 ) ;
   HESSIANS  = PEL_Vector::create( a_owner, 0 ) ;
   BD_FIELDS = PEL_Vector::create( a_owner, 0 ) ;
   BD_VALUES = PEL_Vector::create( a_owner, 0 ) ;
   CONTEXT   = PEL_ContextSimple::create( a_owner ) ;
   COORDS    = PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) ;
   CONTEXT->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
}

//---------------------------------------------------------------------------
void
PDE_LocalFE_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE_TEST:: process_one_test" ) ;

   MY_EPS = exp->double_data( "dbl_epsilon" ) ;
   MY_MIN = exp->double_data( "dbl_minimum" ) ;
   
   PDE_DomainAndFields* dom = PDE_DomainAndFields::create( 0, exp/*, PEL_Exec::communicator()*/ ) ;

   initialize_internals( dom ) ;

   PDE_SetOfDiscreteFields const* dfs = dom->set_of_discrete_fields() ;

   PDE_AdapterCHARMS* da = dom->adapter_CHARMS() ;
   if( da != 0 )
   {
      da->reset() ;
      bool keep_adapting = false ;
      do
      {
         da->adapt() ;
         keep_adapting = da->something_changed() ;
         if( keep_adapting ) 
         {
            dom->apply_requests_of_DOFs_values_modules( true ) ;
         }
         da->print_statistics( PEL::out(), 3 ) ;
      }
      while( keep_adapting ) ;
   }
   
   PEL_ModuleExplorer* e_sol = exp->create_subexplorer( 0, "solution" ) ;
   PEL_ModuleExplorer* ee_ff =
                             exp->create_subexplorer( 0, "interior_fields" ) ;
   ee_ff->start_module_iterator() ;
   for( ; ee_ff->is_valid_module() ; ee_ff->go_next_module() )
   {
      PEL_ModuleExplorer const* ee = ee_ff->create_subexplorer( 0 ) ;
      std::string nn = ee->string_data( "name" ) ;
      std::string nn2 = nn + "_dupl" ;
      dom->duplicate_field( nn, nn2 ) ;
      PDE_DiscreteField const* ff = dfs->item( nn ) ;
      PDE_DiscreteField const* ff2 = dfs->item( nn2 ) ;
      I_FIELDS->append( const_cast<PDE_DiscreteField*>( ff ) ) ;
      I_FIELDS->append( const_cast<PDE_DiscreteField*>( ff2 ) ) ;
      PEL_Data const* val = ee->abstract_data( I_VALUES,
                                               "DOFs_values/value", CONTEXT ) ;
      if( val->data_type()!=PEL_Data::DoubleVector )
      {
         PEL_Error::object()->raise_bad_data_type(
            ee, "DOFs_values/value", PEL_Data::DoubleVector  ) ;
      }
      if( !val->value_can_be_evaluated(0) )
      {
         PEL_Error::object()->raise_not_evaluable(
            ee, "DOFs_values/value", val->undefined_variables(0) ) ;
      }
      I_VALUES->append( const_cast<PEL_Data*>( val ) ) ;
      I_VALUES->append( const_cast<PEL_Data*>( val ) ) ;
      PEL_Data const* jac = e_sol->abstract_data( JACOBIANS, 
                                                  nn+"/jacobian", CONTEXT ) ;
      if( jac->data_type()!=PEL_Data::DoubleArray2D )
      {
         PEL_Error::object()->raise_bad_data_type(
            e_sol, nn+"/jacobian", PEL_Data::DoubleArray2D  ) ;
      }
      if( !jac->value_can_be_evaluated(0) )
      {
         PEL_Error::object()->raise_not_evaluable(
            e_sol, nn+"/jacobian", jac->undefined_variables(0) ) ;
      }
      JACOBIANS->append( const_cast<PEL_Data*>( jac ) ) ;
      JACOBIANS->append( const_cast<PEL_Data*>( jac ) ) ;
      PEL_Data const* hes = e_sol->abstract_data( HESSIANS,
                                                  nn+"/hessian", CONTEXT ) ;
      if( hes->data_type()!=PEL_Data::DoubleArray3D )
      {
         PEL_Error::object()->raise_bad_data_type(
            e_sol, nn+"/hessian", PEL_Data::DoubleArray3D  ) ;
      }
      if( !hes->value_can_be_evaluated(0) )
      {
         PEL_Error::object()->raise_not_evaluable(
            e_sol, nn+"/hessian", hes->undefined_variables(0) ) ;
      }
      HESSIANS->append( const_cast<PEL_Data*>( hes ) ) ;
      HESSIANS->append( const_cast<PEL_Data*>( hes ) ) ;
      ee->destroy() ;
   }
   ee_ff->destroy() ;
   
   if( exp->has_module( "boundary_fields" ) )
   {
      ee_ff = exp->create_subexplorer( 0, "boundary_fields" ) ;
      ee_ff->start_module_iterator() ;
      for( ; ee_ff->is_valid_module() ; ee_ff->go_next_module() )
      {
         PEL_ModuleExplorer const* ee = ee_ff->create_subexplorer( 0 ) ;
         std::string nn = ee->string_data( "name" ) ;
         std::string nn2 = nn + "_dupl" ;
         dom->duplicate_field( nn, nn2 ) ;
         PDE_DiscreteField const* ff = dfs->item( nn ) ;
         PDE_DiscreteField const* ff2 = dfs->item( nn2 ) ;
         BD_FIELDS->append( const_cast<PDE_DiscreteField*>( ff ) ) ;
         BD_FIELDS->append( const_cast<PDE_DiscreteField*>( ff2 ) ) ;
         PEL_Data const* val = ee->abstract_data( BD_VALUES,
                                               "DOFs_values/value", CONTEXT ) ;
         if( val->data_type()!=PEL_Data::DoubleVector )
         {
            PEL_Error::object()->raise_bad_data_type(
               ee, "DOFs_values/value", PEL_Data::DoubleVector  ) ;
         }
         if( !val->value_can_be_evaluated(0) )
         {
            PEL_Error::object()->raise_not_evaluable(
               ee, "DOFs_values/value", val->undefined_variables(0) ) ;
         }
         BD_VALUES->append( const_cast<PEL_Data*>( val ) ) ;
         BD_VALUES->append( const_cast<PEL_Data*>( val ) ) ;
         ee->destroy() ;
      }
      ee_ff->destroy() ;
   }
   
   CFE = dom->create_LocalFEcell( dom ) ;
   SFE = dom->create_CursorFEside( dom ) ;
   BFE = dom->create_LocalFEbound( dom ) ;
   PDE_LocalFEbound* bFE = dom->create_LocalFEbound( FES ) ;
   FES->append( bFE ) ;
   PDE_LocalFEcell* cFE = dom->create_LocalFEcell( FES ) ;
   FES->append( cFE ) ;
   if( e_sol->has_module( "interface" ) )
   {
      PEL_ModuleExplorer* ee = e_sol->create_subexplorer( 0, "interface" ) ;
      stringVector const& ic = ee->stringVector_data( "inward_colors" ) ;
      stringVector const& oc = ee->stringVector_data( "outward_colors" ) ;
      for( size_t i=0 ; i<ic.size() ; ++i )
      {
         GE_Color const* c1 = GE_Color::object( ic( i ) ) ;
         GE_Color const* c2 = GE_Color::object( oc( i ) ) ;
         FES->append( dom->create_LocalFEinterface( FES, c2, c1 ) ) ;
         FES->append( dom->create_LocalFEinterface( FES, c1, c2 ) ) ;
      }
      ee->destroy() ;
   }
   e_sol->destroy() ;

   if( I_FIELDS->count() > 0 )
   {
      for( size_t i=0 ; i<FES->count() ; ++i )
      {
         PDE_LocalFE* fe = static_cast<PDE_LocalFE*>( FES->at(i) ) ;
         for( size_t j=0 ; j<I_FIELDS->count() ; ++j )
         {
            PDE_DiscreteField const* ff = 
                        static_cast<PDE_DiscreteField*>( I_FIELDS->at(j) ) ;
            fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
            fe->require_field_calculation( ff, PDE_LocalFE::dN ) ;
            fe->require_field_calculation( ff, PDE_LocalFE::d2N ) ;
         }
         iterate_over_meshes( exp->name(), fe, dom->nb_space_dimensions() ) ;
      }
   }
   if( BD_FIELDS->count() > 0 )
   {
      for( size_t j=0 ; j<BD_FIELDS->count() ; ++j )
      {
         PDE_DiscreteField const* ff = 
                     static_cast<PDE_DiscreteField*>( BD_FIELDS->at(j) ) ;
         bFE->require_field_calculation( ff, PDE_LocalFE::N ) ;
         cFE->require_field_calculation( ff, PDE_LocalFE::node ) ;
      }
      iterate_over_bounds( exp->name(), bFE, dom->nb_space_dimensions() ) ;
   }

   iterate_over_cells( exp->name(), cFE ) ;

   dom->destroy() ;
   CFE = 0 ; SFE = 0 ; BFE = 0 ;
}

//----------------------------------------------------------------------
void
PDE_LocalFE_TEST:: iterate_over_meshes( std::string const& test_name,
                                        PDE_LocalFE* fe,
                                        size_t nb_sp_dims )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE_TEST:: iterate_over_meshes" ) ;

   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_3" ) ;

   string fe_name = test_name+"::"+fe->type_name() ;

   double theo, xx ;
   bool eq ;

   bool ok_val  = true ;
   bool ok_grad = true ;
   bool ok_hes  = true ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      fe->start_IP_iterator( qrp ) ;
      for( ; fe->valid_IP() ; fe->go_next_IP() )
      {
         fill_context_with_coordinates( fe->coordinates_of_IP() ) ;
         for( size_t i=0 ; i<I_FIELDS->count() ; ++i )
         {
            PDE_DiscreteField const* ff = 
                           static_cast<PDE_DiscreteField*>( I_FIELDS->at(i) ) ;
            PEL_Data const* val = static_cast<PEL_Data*>( I_VALUES->at(i) ) ;
            PEL_Data const* jac = static_cast<PEL_Data*>( JACOBIANS->at(i) ) ;
            PEL_Data const* hes = static_cast<PEL_Data*>( HESSIANS->at(i) ) ;
            for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
            {
               theo = val->to_double_vector()( ic ) ;
               xx = fe->value_at_IP( ff, 0, ic ) ;
               eq = PEL::double_equality( theo, xx, MY_EPS, MY_MIN ) ;
               ok_val = ok_val && eq ;
               if( !eq ) display_error( theo, xx ) ;

               for( size_t a=0 ; a<nb_sp_dims ; ++a )
               {
                  theo = jac->to_double_array2D()( ic, a ) ;
                  xx = fe->gradient_at_IP( ff, 0, a, ic ) ;
                  eq = PEL::double_equality( theo, xx, MY_EPS, MY_MIN ) ;
                  ok_grad = ok_grad && eq ;
                  if( !eq ) display_error( theo, xx ) ;

                  for( size_t b=0 ; b<nb_sp_dims ; ++b )
                  {
                     theo = hes->to_double_array3D()( ic, a, b ) ;
                     xx = fe->hessian_at_IP( ff, 0, a, b, ic ) ;
                     eq = PEL::double_equality( theo, xx, MY_EPS, MY_MIN ) ;
                     ok_hes = ok_hes && eq ;
                     if( !eq ) display_error( theo, xx ) ;
                  }
               }
            }
         }
      }
   }
   notify_one_test_result( fe_name+"::value_at_IP", ok_val ) ;
   notify_one_test_result( fe_name+"::gradient_at_IP", ok_grad ) ;
   notify_one_test_result( fe_name+"::hessian_at_IP", ok_hes ) ;

   ok_val  = true ;
   ok_grad = true ;
   ok_hes  = true ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      GE_Point const* pt = fe->polyhedron()->center() ;
      fe->set_calculation_point( pt ) ;
      fill_context_with_coordinates( pt ) ;
      for( size_t i=0 ; i<I_FIELDS->count() ; ++i )
      {
         PDE_DiscreteField const* ff = 
                           static_cast<PDE_DiscreteField*>( I_FIELDS->at(i) ) ;
         PEL_Data const* val = static_cast<PEL_Data*>( I_VALUES->at(i) ) ;
         PEL_Data const* jac = static_cast<PEL_Data*>( JACOBIANS->at(i) ) ;
         PEL_Data const* hes = static_cast<PEL_Data*>( HESSIANS->at(i) ) ;
         for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
         {
            theo = val->to_double_vector()( ic ) ;
            xx = fe->value_at_pt( ff, 0, ic ) ;
            eq = PEL::double_equality( theo, xx, MY_EPS, MY_MIN ) ;
            ok_val = ok_val && eq ;
            if( !eq ) display_error( theo, xx ) ;
               
            for( size_t a=0 ; a<nb_sp_dims ; ++a )
            {
               theo = jac->to_double_array2D()( ic, a ) ;
               xx = fe->gradient_at_pt( ff, 0, a, ic ) ;
               eq = PEL::double_equality( theo, xx, MY_EPS, MY_MIN ) ;
               ok_grad = ok_grad && eq ;
               if( !eq ) display_error( theo, xx ) ;

               for( size_t b=0 ; b<nb_sp_dims ; ++b )
               {
                  theo = hes->to_double_array3D()( ic, a, b ) ;
                  xx = fe->hessian_at_pt( ff, 0, a, b, ic ) ;
                  eq = PEL::double_equality( theo, xx, MY_EPS, MY_MIN ) ;
                  ok_hes = ok_hes && eq ;
                  if( !eq ) display_error( theo, xx ) ;
               }
            }
         }
      }
   }
   notify_one_test_result( fe_name+"::value_at_pt", ok_val ) ;
   notify_one_test_result( fe_name+"::gradient_at_pt", ok_grad ) ;
   notify_one_test_result( fe_name+"::hessian_at_pt", ok_hes ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFE_TEST:: iterate_over_bounds( std::string const& test_name,
                                        PDE_LocalFEbound* fe, 
                                        size_t nb_sp_dims )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE_TEST:: iterate_over_bounds" ) ;

   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_3" ) ;

   double theo, xx ;
   bool eq ;

   bool ok_val  = true ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      fe->start_IP_iterator( qrp ) ;
      for( ; fe->valid_IP() ; fe->go_next_IP() )
      {
         fill_context_with_coordinates( fe->coordinates_of_IP() ) ;
         for( size_t i=0 ; i<BD_FIELDS->count() ; ++i )
         {
            PDE_DiscreteField const* ff = 
                          static_cast<PDE_DiscreteField*>( BD_FIELDS->at(i) ) ;
            PEL_Data const* val = static_cast<PEL_Data*>( BD_VALUES->at(i) ) ;
            for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
            {
               theo = val->to_double_vector()( ic ) ;
               xx = fe->value_at_IP( ff, 0, ic ) ;
               eq = PEL::double_equality( theo, xx, MY_EPS, MY_MIN ) ;
               ok_val = ok_val && eq ;
               if( !eq ) display_error( theo, xx ) ;
            }
         }
      }
   }
   notify_one_test_result( test_name+"::PDE_LocalFEbound::value_at_IP",
                           ok_val ) ;

   ok_val = true ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      GE_Point const* pt = fe->polyhedron()->center() ;
      fe->set_calculation_point( pt ) ;
      fill_context_with_coordinates( pt ) ;
      for( size_t i=0 ; i<BD_FIELDS->count() ; ++i )
      {
         PDE_DiscreteField const* ff = 
                          static_cast<PDE_DiscreteField*>( BD_FIELDS->at(i) ) ;
         PEL_Data const* val = static_cast<PEL_Data*>( BD_VALUES->at(i) ) ;
         for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
         {
            theo = val->to_double_vector()( ic ) ;
            xx = fe->value_at_pt( ff, 0, ic ) ;
            eq = PEL::double_equality( theo, xx, MY_EPS, MY_MIN ) ;
            ok_val = ok_val && eq ;
            if( !eq ) display_error( theo, xx ) ;
         }
      }
   }
   notify_one_test_result( test_name+"::PDE_LocalFEbound::value_at_pt",
                           ok_val ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFE_TEST:: iterate_over_cells( std::string const& test_name,
                                       PDE_LocalFEcell* fe )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE_TEST:: iterate_over_cells" ) ;

   string fe_name = test_name+"::"+fe->type_name() ;

   bool ok_cells = true ;
   bool ok_sides = true ;
   bool ok_bounds = true ;
   bool ok_bd_fields = true ;
   
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      size_t_vector const& cell_ids = fe->adjacent_cell_ids() ;
      for( size_t i=0 ; i<cell_ids.size() ; ++i )
      {
         CFE->go_i_th( cell_ids( i ) ) ;
         ok_cells &= CFE->adjacent_cell_ids().has( fe->mesh_id() ) ;
      }
      size_t_vector const& side_ids = fe->adjacent_side_ids() ;
      for( size_t i=0 ; i<side_ids.size() ; ++i )
      {
         SFE->go_i_th( side_ids( i ) ) ;
         size_t const cell0_id = SFE->adjacent_localFEcell( 0 )->mesh_id() ;
         size_t const cell1_id = SFE->adjacent_localFEcell( 1 )->mesh_id() ;
         ok_sides &= cell0_id == fe->mesh_id() || cell1_id == fe->mesh_id() ;
         if( cell0_id == fe->mesh_id() )
         {
            ok_cells &= cell_ids.has( cell1_id ) ;
         }
         if( cell1_id == fe->mesh_id() )
         {
            ok_cells &= cell_ids.has( cell0_id ) ;
         }
      }
      size_t_vector const& bound_ids = fe->adjacent_bound_ids() ;
      for( size_t i=0 ; i<bound_ids.size() ; ++i )
      {
         BFE->go_i_th( bound_ids( i ) ) ;
         ok_bounds &= BFE->adjacent_cell_id() == fe->mesh_id() ;
         ok_bounds &= BFE->adjacent_cell_color() == fe->color() ;
         ok_bounds &= BFE->adjacent_cell_polyhedron() == fe->polyhedron() ;
      }
      for( size_t i=0 ; ok_bd_fields && i<BD_FIELDS->count() ; ++i )
      {
         PDE_DiscreteField const* ff =
                   static_cast<PDE_DiscreteField const*>( BD_FIELDS->at(i) ) ;
         ok_bd_fields = ( fe->nb_local_nodes( ff ) == 0 ) ;
      }
   }
   notify_one_test_result( fe_name+"::adjacent_cell_ids", ok_cells ) ;
   notify_one_test_result( fe_name+"::adjacent_side_ids", ok_sides ) ;
   notify_one_test_result( fe_name+"::adjacent_bound_ids", ok_bounds ) ;
   if( BD_FIELDS->count() != 0 )
   {
      notify_one_test_result(
            fe_name+"::boundary fields discretization", ok_bd_fields ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFE_TEST:: fill_context_with_coordinates( GE_Point const* pt )
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
void PDE_LocalFE_TEST:: display_error( double theo, double xx ) const
//-----------------------------------------------------------------------
{
   std::cout << " theo " << theo << " xx " << xx
             << " relative " << PEL::relative( theo, xx )
             << std::endl ;
}



