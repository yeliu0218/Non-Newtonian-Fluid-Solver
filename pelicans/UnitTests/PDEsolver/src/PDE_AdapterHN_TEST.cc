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

#include <PDE_AdapterHN_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>

#include <PDE_AdapterHN.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_CursorFEside.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::ios_base ;
using std::ostringstream ;
using std::setprecision ; using std::setw ;
using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

PDE_AdapterHN_TEST*
PDE_AdapterHN_TEST:: REGISTRATOR = new PDE_AdapterHN_TEST() ;

//---------------------------------------------------------------------------
PDE_AdapterHN_TEST:: PDE_AdapterHN_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_AdapterHN", "PDE_AdapterHN_TEST" )
   , QRP( 0 )
{
}

//---------------------------------------------------------------------------
PDE_AdapterHN_TEST:: ~PDE_AdapterHN_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_AdapterHN_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN_TEST:: process_one_test" ) ;

   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "iterations_trace" ) ;
   std::ofstream ofs( se->string_data( "output_file" ).c_str() ) ;
   se->destroy() ;

   se = exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom =
       PDE_DomainAndFields::create( 0, se, PEL_Exec::communicator() ) ;
   se->destroy() ; se = 0 ;
   PDE_AdapterHN* da = dom->adapter_HN() ;

   TEST_NAME = exp->name() ;
   
   if( exp->has_entry( "quadrature_rule_provider" ) )
   {
      QRP = GE_QRprovider::object(
                          exp->string_data( "quadrature_rule_provider" ) ) ;
   }

   da->reset() ;

   da->adapt() ;
   process_domain( 1, dom, ofs ) ;

   dom->destroy() ;
   ofs.close() ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterHN_TEST:: process_domain( size_t ref_step,
                                     PDE_DomainAndFields const* dom,
                                     std::ostream& os )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN_TEST:: process_domain" ) ;

   os << "********* GRID *********************" << endl ;
   dom->print_grid( os, 0 ) ;

   PDE_SetOfDiscreteFields* sdf = dom->set_of_discrete_fields() ;

   os << "********* CELLS ********************" << endl ;
   PDE_LocalFEcell* cfe = dom->create_LocalFEcell( 0 ) ;
   cfe->include_color( GE_Color::halo_color() ) ;
   iterate( ref_step, sdf, cfe, os ) ;
   cfe->destroy() ; cfe = 0 ;

   os << "********* BOUNDS ********************" << endl ;
   PDE_LocalFEbound* bfe = dom->create_LocalFEbound( 0 ) ;
   iterate( ref_step, sdf, bfe, os ) ;
   bfe->destroy() ; bfe = 0 ;

   os << "********* SIDES ********************" << endl ;
   PDE_CursorFEside* sfe = dom->create_CursorFEside( 0 ) ;
   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      sfe->require_field_calculation( ff, PDE_LocalFE::N ) ;
   }
   for( sfe->start() ; sfe->is_valid() ; sfe->go_next() )
   {
      os << "------------------------------" << endl ;
      sfe->print_current_mesh( os, 0 ) ;
   }
   sfe->destroy() ; sfe=0 ;
   
   iterate_fields( ref_step, sdf, os ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterHN_TEST:: iterate( size_t ref_step,
                              PDE_SetOfDiscreteFields* sdf,
                              PDE_LocalFE* fe,
                              std::ostream& os )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN_TEST:: iterate" ) ;

   bool ok_loc_node = true ;

   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
   }
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      os << "------------------------------" << endl ;
      check_node_location( ref_step, sdf, fe, ok_loc_node ) ;
      if( QRP == 0 )
      {
         fe->print( os, 0 ) ;
      }
      else
      {
         os << "current mesh" << endl ;
         fe->print_current_mesh( os, 3 ) ;
         os << "Local discretization of fields" << endl ;
         fe->print_local_discretization_of_fields( os, 3 ) ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
            fe->set_row_and_col_fields( ff, ff ) ;
            os << "Integration points" << endl ;
            fe->start_IP_iterator( QRP ) ;
            for( ; fe->valid_IP() ; fe->go_next_IP() )
            {
               fe->print_values_at_current_IP( os, 3 ) ;
            }
         }
      }
   }

   ostringstream mesg ;
   mesg << TEST_NAME << "(" << ref_step << ") "
        << fe->type_name() << " node location";
   notify_one_test_result( mesg.str(), ok_loc_node ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterHN_TEST:: check_node_location( size_t ref_step,
                                          PDE_SetOfDiscreteFields* sdf,
                                          PDE_LocalFE* fe,
                                          bool& ok_node_loc ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN_TEST:: check_node_location" ) ;

   GE_Mpolyhedron const* poly = fe->polyhedron() ;

   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      for( size_t i=0 ; i<fe->nb_local_nodes( ff ) ; ++i )
      {
         GE_Point const* pt = fe->local_node_location( ff, i ) ;
         bool ok =
              ( poly->contains( pt ) &&  fe->local_node_is_in_mesh( ff, i )) ||
              (!poly->contains( pt ) && !fe->local_node_is_in_mesh( ff, i ) ) ;
         if( !ok )
         {
            PEL::out() << "------------------------------------" << endl ;
            PEL::out() << "current mesh" << endl ;
            fe->print_current_mesh( PEL::out(), 3 ) ;
            PEL::out() << "Local discretization of fields" << endl ;
            fe->print_local_discretization_of_fields( PEL::out(), 3 ) ;
            PEL::out() << "------" << endl ;
            PEL::out() << "node : " ;
            pt->print( PEL::out(), 3 ) ; PEL::out() << endl ;
            PEL::out() << "   local_node_is_in_mesh --> "
                       << fe->local_node_is_in_mesh( ff, i ) << endl ;
            PEL::out() << "                contains --> "
                       << poly->contains( pt ) << endl ;
         }
         ok_node_loc = ok_node_loc && ok ;
      }
   }
}

//---------------------------------------------------------------------------
void
PDE_AdapterHN_TEST:: iterate_fields( size_t ref_step,
                                     PDE_SetOfDiscreteFields* sdf,
                                     std::ostream& os ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN_TEST:: iterate_fields" ) ;
   
   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      os << "********* \"" << ff->name() << "\" ********************" << endl ;
      print_field( ff, os ) ;
      PDE_LinkDOF2Unknown* ff_link = PDE_LinkDOF2Unknown::create( 0, ff,
                                               "sequence_of_the_components",
                                               true ) ;
      PDE_SystemNumbering* nmb = PDE_SystemNumbering::create( 0, ff_link ) ;
      os << "-------------" << std::endl ;
      nmb->print( os, 3 ) ;
      nmb->destroy() ;
   }   
}

//-----------------------------------------------------------------------
void
PDE_AdapterHN_TEST:: display_error( PDE_LocalFE const* fe,
                                    double theo, double xx ) const
//-----------------------------------------------------------------------
{
   std::cout << fe->type_name() << " mesh : " << fe->mesh_id()
             << " theo " << theo << " xx " << xx << std::endl ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterHN_TEST:: print_field( PDE_DiscreteField const* ff,
                                  std::ostream& os ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterHN_TEST:: print_field" ) ;
   
   PEL_ASSERT( ff->nb_components() == 1 ) ;
   size_t ic = 0 ;
   
   for( size_t n=0 ; n<ff->nb_nodes() ; ++n )
   {
      os << "node " << setw( 3 ) << n ;
      os << " " ;
      if( !ff->node_is_active( n ) )
      {
         os << "inactive" ;
      }
      else
      {
         os << "  active" ;
      }
      os << " " ;
      if( ff->DOF_has_imposed_value( n, ic ) )
      {
         os << "imposed" ;
      }
      else
      {
         os << "       " ;
      }
      os << " " ;
      if( ff->DOF_is_constrained( n, ic ) )
      {
         os << "constrained" ;
      }
      else
      {
         os << "           " ;         
      }
      os << std::endl ;
   }
}
