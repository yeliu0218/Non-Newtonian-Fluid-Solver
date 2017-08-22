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

#include <PDE_FieldInitialization_TEST.hh>

#include <PEL.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <doubleVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>

#include <PDE_AdapterCHARMS.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalFEinterface.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <vector>
#include <iostream>

using std::cout ;
using std::endl ;
using std::string ;

PDE_FieldInitialization_TEST const*
PDE_FieldInitialization_TEST:: REGISTRATOR = 
                               new PDE_FieldInitialization_TEST() ;

//---------------------------------------------------------------------------
PDE_FieldInitialization_TEST:: PDE_FieldInitialization_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_LocalFE", "PDE_FieldInitialization_TEST" ),
     CONTEXT( 0 ),
     COORDS( 0 )
{
   CONTEXT   = PEL_ContextSimple::create( this ) ;
   COORDS    = PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) ;
   CONTEXT->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
}

//---------------------------------------------------------------------------
PDE_FieldInitialization_TEST:: ~PDE_FieldInitialization_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_FieldInitialization_TEST:: process_one_test( 
                                               PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldInitialization_TEST:: process_one_test" ) ;

   PEL_ModuleExplorer* se = 
             exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom = PDE_DomainAndFields::create( 0, se ) ;
   se->destroy() ; se = 0 ;
   
   PDE_SetOfDiscreteFields const* dfs =  dom->set_of_discrete_fields() ;

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
   
   std::vector< std::string > names ;
   for( dfs->start() ; dfs->is_valid() ; dfs->go_next() )
   {
      names.push_back( dfs->item()->name() ) ;
   }
   std::vector< std::string >::const_iterator it = names.begin() ;
   for( ; it != names.end() ; ++it )
   {
      dom->duplicate_field( (*it), (*it)+"_dupl" ) ;
   }

   PDE_LocalFEbound* bFE = dom->create_LocalFEbound( 0 ) ;
   PDE_LocalFEcell* cFE = dom->create_LocalFEcell( 0 ) ;

   if( exp->has_module( "solution" ) )
   {
      PEL_ModuleExplorer const* eee = exp->create_subexplorer( 0, "solution" ) ;
      if( eee->has_module( "interior_fields" ) )
      {
         se = eee->create_subexplorer( 0, "interior_fields" ) ;
         se->start_module_iterator() ;
         for( ; se->is_valid_module() ; se->go_next_module() )
         {
            PEL_ModuleExplorer* e = se->create_subexplorer( 0 ) ;
            std::string nn = e->string_data( "name" ) ;
            std::string nn2 = nn + "_dupl" ;
            PDE_DiscreteField const* ff  = dfs->item( nn  ) ;
            PDE_DiscreteField const* ff2 = dfs->item( nn2 ) ;
            PEL_Data const* val = e->abstract_data( e, "value", CONTEXT ) ;
            if( val->data_type()!=PEL_Data::DoubleVector )
            {
               PEL_Error::object()->raise_bad_data_type( 
                                    e, "value", PEL_Data::DoubleVector  ) ;
            }
            if( !val->value_can_be_evaluated(0) )
            {
               PEL_Error::object()->raise_not_evaluable(
                                    e, "value", val->undefined_variables(0) ) ;
            }
            notify_one_test_result( "field \""+nn+"\" on cells",
                                    check_from_meshes( cFE, ff, val ) ) ;
            notify_one_test_result( "field \""+nn2+"\" on cells",
                                    check_from_meshes( cFE, ff2, val ) ) ;
            notify_one_test_result( "field \""+nn+"\" on bounds",
                                    check_from_meshes( bFE, ff, val ) ) ;
            notify_one_test_result( "field \""+nn2+"\" on bounds",
                                    check_from_meshes( bFE, ff2, val ) ) ;
            e->destroy() ; e = 0 ;
         }
         se->destroy() ; se = 0 ;
      }

      if( eee->has_module( "boundary_fields" ) )
      {
         se = eee->create_subexplorer( 0, "boundary_fields" ) ;
         se->start_module_iterator() ;
         for( ; se->is_valid_module() ; se->go_next_module() )
         {
            PEL_ModuleExplorer* e = se->create_subexplorer( 0 ) ;
            std::string nn = e->string_data( "name" ) ;
            std::string nn2 = nn + "_dupl" ;
            PDE_DiscreteField const* ff = dfs->item( nn ) ;
            PDE_DiscreteField const* ff2 = dfs->item( nn2 ) ;
            PEL_Data const* val = e->abstract_data( e, "value", CONTEXT ) ;
            if( val->data_type()!=PEL_Data::DoubleVector )
            {
               PEL_Error::object()->raise_bad_data_type(
                                    e, "value", PEL_Data::DoubleVector  ) ;
            }
            if( !val->value_can_be_evaluated(0) )
            {
               PEL_Error::object()->raise_not_evaluable(
                                    e, "value", val->undefined_variables(0) ) ;
            }
            notify_one_test_result( "field \""+nn+"\" on bounds",
                                    check_from_meshes( bFE, ff, val ) ) ;
            notify_one_test_result( "field \""+nn2+"\" on bounds",
                                    check_from_meshes( bFE, ff2, val ) ) ;
            e->destroy() ;
         }
         se->destroy() ; se = 0 ;
      }
      eee->destroy() ; eee = 0 ;
   }

   bFE->destroy() ;
   cFE->destroy() ;
   dom->destroy() ;
}

//----------------------------------------------------------------------
bool
PDE_FieldInitialization_TEST:: check_from_meshes( PDE_LocalFE* fe,
                                                  PDE_DiscreteField const* ff,
                                                  PEL_Data const* val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldInitialization_TEST:: check_from_meshes" ) ;

   bool result  = true ;

   fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      for( size_t il=0 ; il<fe->nb_local_nodes(ff) ; ++il )
      {
         size_t nn = fe->global_node( ff, il ) ;
         GE_Point const* pt = fe->local_node_location( ff, il ) ;
         size_t ref_level = fe->node_refinement_level( ff, il ) ;
         fill_context_with_coordinates( pt ) ;
         doubleVector const& v_theo = val->to_double_vector() ;
         if( ref_level == 0 )
         {
            for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
            {
               double theo = v_theo( ic ) ;
               for( size_t level=0 ; level<ff->storage_depth() ; ++level )
               {
                  double xx = ff->DOF_value( level, nn, ic ) ;
                  bool eq = PEL::equal( theo, xx ) ;
                  result = result && eq ;
                  if( !eq ) display_error( pt, ic, theo, xx ) ;
               }
            }
         }
         if( fe->local_node_is_in_mesh( ff, il ) )
         {
            fe->set_calculation_point( pt ) ;
            for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
            {
               double theo = v_theo( ic ) ;
               for( size_t level=0 ; level<ff->storage_depth() ; ++level )
               {
                  double xx = fe->value_at_pt( ff, level, ic ) ;
                  bool eq = PEL::equal( theo, xx ) ;
                  result = result && eq ;
                  if( !eq ) display_error( pt, ic, theo, xx ) ;
               }
            }
         }
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_FieldInitialization_TEST:: fill_context_with_coordinates( 
                                                          GE_Point const* pt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldInitialization_TEST:: fill_context_with_coordinates" );
   size_t n = pt->nb_coordinates() ;
   doubleVector dv( n ) ;
   for( size_t i=0 ; i<n ; ++i )
   {
      dv(i) = pt->coordinate( i ) ;
   }
   COORDS->set( dv ) ;
}

//-----------------------------------------------------------------------
void
PDE_FieldInitialization_TEST:: display_error( GE_Point const* pt,
                                              size_t ic,
                                              double xx_data_deck,
                                              double xx_calculated ) const
//------------------------------------------------------------------------
{
   std::cout << "node : " ;
   pt->print( cout, 0 ) ;
   std::cout << " component : " << ic << endl ; 
   std::cout << "   from data deck : " << xx_data_deck ;
   std::cout << " computed : " <<  xx_calculated << std::endl ;
}

