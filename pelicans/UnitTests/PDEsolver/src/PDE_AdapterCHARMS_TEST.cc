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

#include <PDE_AdapterCHARMS_TEST.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_System.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>
#include <doubleArray2D.hh>
#include <stringVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>
#include <GE_SegmentPolyhedron_INT.hh>
#include <GE_SetOfPoints.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_AdapterCHARMS.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_MeshingCoarsener.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfDomains.hh>

#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::ios_base ;
using std::ostringstream ;
using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

PDE_AdapterCHARMS_TEST const*
PDE_AdapterCHARMS_TEST:: REGISTRATOR = new PDE_AdapterCHARMS_TEST() ;

//---------------------------------------------------------------------------
PDE_AdapterCHARMS_TEST:: PDE_AdapterCHARMS_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_AdapterCHARMS", "PDE_AdapterCHARMS_TEST" )
   , QRP( 0 )
   , VAL_CHECK( 0 )
   , UU_CHECK( 0 )
   , EPS_INTERP( PEL::bad_double() )
   , MIN_INTERP( PEL::bad_double() )
   , FFS_VERTS( 0 )
   , EPS_VERTS( PEL::bad_double() )
   , MIN_VERTS( PEL::bad_double() )
   , CONTEXT( 0 )
   , COORDS( 0 )
   , VAL_BDIMP( 0 )
   , UU_BDIMP( 0 )
   , COLOR_BDIMP( 0 )
   , EPS_BDIMP( PEL::bad_double() )
   , MIN_BDIMP( PEL::bad_double() )
{
}

//---------------------------------------------------------------------------
PDE_AdapterCHARMS_TEST:: ~PDE_AdapterCHARMS_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: process_one_test" ) ;

   CHECK_BD_CELLS = exp->has_entry( "check_covering_of_cells_boundary" ) ?
                     exp->bool_data( "check_covering_of_cells_boundary" ) :
                    true ;
   QRP = 0 ;
   VAL_CHECK=0 ; UU_CHECK=0 ;
   EPS_INTERP=PEL::bad_index() ; MIN_INTERP=PEL::bad_index() ;

   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "iterations_trace" ) ;
   std::ofstream ofs( se->string_data( "output_file" ).c_str() ) ;
   se->destroy() ;

   se = exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom =
       PDE_DomainAndFields::create( 0, se, PEL_Exec::communicator() ) ;
   se->destroy() ; se = 0 ;
   PDE_AdapterCHARMS* da = dom->adapter_CHARMS() ;

   da->reset() ;

   CONTEXT = PEL_ContextSimple::create( dom ) ;
   COORDS  = PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) ;
   CONTEXT->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;

   TEST_NAME = exp->name() ;

   if( exp->has_entry( "quadrature_rule_provider" ) )
   {
      QRP = GE_QRprovider::object(
                          exp->string_data( "quadrature_rule_provider" ) ) ;
   }

   if( exp->has_module( "interpolation_check" ) )
   {
      se = exp->create_subexplorer( 0, "interpolation_check" ) ;
      UU_CHECK =
         dom->set_of_discrete_fields()->item( se->string_data( "field" ) ) ;
      VAL_CHECK = se->abstract_data( dom, "value", CONTEXT ) ;
      if( VAL_CHECK->data_type()!=PEL_Data::DoubleVector )
      {
         PEL_Error::object()->raise_bad_data_type(
                                    se, "value", PEL_Data::DoubleVector  ) ;
      }
      if( !VAL_CHECK->value_can_be_evaluated() )
      {
         PEL_Error::object()->raise_not_evaluable(
                           se, "value", VAL_CHECK->undefined_variables() ) ;
      }
      EPS_INTERP = se->double_data( "dbl_epsilon" ) ;
      MIN_INTERP = se->double_data( "dbl_minimum" ) ;
      se->destroy() ; se = 0 ;

      da->add_excluded_field( UU_CHECK ) ;
   }

   if( exp->has_module( "values_at_vertices_check" ) )
   {
      se = exp->create_subexplorer( 0, "values_at_vertices_check" ) ;
      FFS_VERTS = se->stringVector_data( "fields" ) ;
      EPS_VERTS = se->double_data( "dbl_epsilon" ) ;
      MIN_VERTS = se->double_data( "dbl_minimum" ) ;
      se->destroy() ; se = 0 ;
   }
   
   if( exp->has_module( "imposed_values_check" ) )
   {
      se = exp->create_subexplorer( 0, "imposed_values_check" ) ;
      UU_BDIMP =
         dom->set_of_discrete_fields()->item( se->string_data( "field" ) ) ;
      VAL_BDIMP = se->abstract_data( dom, "value", CONTEXT ) ;
      if( VAL_BDIMP->data_type()!=PEL_Data::DoubleVector )
      {
         PEL_Error::object()->raise_bad_data_type(
                                    se, "value", PEL_Data::DoubleVector  ) ;
      }
      if( !VAL_BDIMP->value_can_be_evaluated() )
      {
         PEL_Error::object()->raise_not_evaluable(
                           se, "value", VAL_CHECK->undefined_variables() ) ;
      }
      COLOR_BDIMP = GE_Color::object( se->string_data( "color" ) ) ;
      EPS_BDIMP = se->double_data( "dbl_epsilon" ) ;
      MIN_BDIMP = se->double_data( "dbl_minimum" ) ;
      se->destroy() ; se = 0 ;
   }

   size_t i=0 ;

   bool keep_adapting = false ;
   do
   {
      ofs << "======================================================" << endl ;
      ofs << "            REFINEMEMENT : " << i << endl ;
      ofs << "======================================================" << endl ;

      da->adapt() ;
      keep_adapting = da->something_changed() ;
      if( keep_adapting )
      {
         process_domain( i, dom, ofs ) ;
      }

      ++i ;
   }
   while( keep_adapting ) ;

   PDE_MeshingCoarsener* coar = da->meshing_coarsener() ;
   check_meshing_coarsening( dom, coar ) ;

   dom->destroy() ;

   FFS_VERTS.re_initialize( 0 ) ; EPS_VERTS = MIN_VERTS = PEL::bad_double() ;
   UU_CHECK = 0 ; VAL_CHECK = 0 ; EPS_INTERP = MIN_INTERP = PEL::bad_double() ;
   UU_BDIMP = 0 ; VAL_BDIMP = 0 ; EPS_BDIMP = MIN_BDIMP = PEL::bad_double() ;
   ofs.close() ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: check_meshing_coarsening(
                                               PDE_DomainAndFields const* dom,
                                               PDE_MeshingCoarsener* coar )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: check_meshing_coarsening" ) ;

   size_t ll = 0 ;
   PDE_SetOfDiscreteFields* sdf = dom->set_of_discrete_fields() ;
   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField* ff = sdf->item() ;
      for( size_t n=0 ; n<ff->nb_nodes() ; ++n )
      {
         for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
         {
            // on affecte une valeur quelconque, qui ne doit pas
            // etre modifiee par les operation de coarsening/uncoarsening
            double xx = (double)n ;
            ff->set_DOF_value( ll, n, xx, ic ) ;
         }
      }
   }

   bool ok_coar = true ;
   coar->prepare_for_coarsening() ;
   size_t level_max = coar->nb_levels()-1 ;

   if( level_max != 0 )
   {
      size_t level = level_max ;
      print_discretization( dom, "a", level ) ;
      for( ; level != 0 ; --level )
      {
         coar->do_one_coarsening() ;
         size_t current_level = level-1 ;
         if( current_level != 0 )
            print_discretization( dom, "a", current_level ) ;
      }
      for( ; level != level_max ; ++level )
      {
         coar->do_one_uncoarsening() ;
         size_t current_level = level+1 ;
         print_discretization( dom, "b", current_level ) ;
         bool ok = true ;
         comp_discretizations( "a", "b", current_level, ok ) ;
         ok_coar = ok_coar && ok ;
      }
   }

   ostringstream mesg ;
   mesg << TEST_NAME << " coarsening/uncoarsening " ;
   notify_one_test_result( mesg.str(), ok_coar ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: print_discretization( PDE_DomainAndFields const* dom,
                                               std::string const& fname,
                                               size_t level )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: print_discretizations" ) ;

   PEL_ASSERT( level != 0 ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;

   std::ostringstream mm ;
   mm << fname << "_" << com->rank() << "_" << level ;
   std::ofstream ofs( mm.str().c_str() ) ;

   ofs << "********* FIELDS *******************" << endl ;
   PDE_SetOfDiscreteFields* sdf = dom->set_of_discrete_fields() ;
   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      ff->print( ofs, 0 ) ;
   }

   ofs << "********* CELLS ********************" << endl ;
   PDE_LocalFEcell* cfe = dom->create_LocalFEcell( 0 ) ;
   iterate_1( level, sdf, cfe, ofs ) ;
   cfe->destroy() ; cfe=0 ;

   ofs << "********* BOUNDS ********************" << endl ;
   PDE_LocalFEbound* bfe = dom->create_LocalFEbound( 0 ) ;
   iterate_1( level, sdf, bfe, ofs ) ;
   bfe->destroy() ; bfe=0 ;

   ofs << "********* SIDES ********************" << endl ;
   PDE_CursorFEside* sfe = dom->create_CursorFEside( 0 ) ;
   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      sfe->require_field_calculation( ff, PDE_LocalFE::N ) ;
   }
   for( sfe->start() ; sfe->is_valid() ; sfe->go_next() )
   {
      ofs << "------------------------------" << endl ;
      sfe->print_current_mesh( ofs, 0 ) ;
   }
   sfe->destroy() ; sfe=0 ;

   if( CHECK_BD_CELLS )
   {
      check_cell_and_faces( level, dom ) ;
      check_covering_of_cells_boundary( level, dom ) ;
   }
   ofs.close() ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: comp_discretizations( std::string const& fname_1,
                                               std::string const& fname_2,
                                               size_t level,
                                               bool& ok )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: comp_discretizations" ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;

   std::ostringstream m1 ;
   m1 << fname_1 << "_" << com->rank() << "_" << level ;
   std::ifstream ifs1( m1.str().c_str() ) ;

   std::ostringstream m2 ;
   m2 << fname_2 << "_" << com->rank() << "_" << level ;
   std::ifstream ifs2( m2.str().c_str() ) ;

   ok = ok && ifs1.good() && ifs2.good() ;
   std::string fline ;
   std::string line1, line2 ;
   size_t iline = 0 ;

   while( ok && ifs1.good() && ifs2.good() )
   {
      iline++ ;
      getline( ifs1, line1 ) ;
      getline( ifs2, line2 ) ;
      if( line1 != line2 )
      {
         ok = false ;
      }
   }
   if( ok && ( ifs1.good() || ifs2.good() ) )
   {
      ok = false ;
   }
   ifs1.close() ;
   ifs2.close() ;

   if( !ok )
   {
      out() << "*** files " <<  m1.str() << " and " << m2.str()
            << " differ ***" << endl ;
   }
   else
   {
      out() << "|     " << TEST_NAME
            << "(" << level << ") recovering discretization : OK " << endl ;
      PEL_System::erase( m1.str() ) ;
      PEL_System::erase( m2.str() ) ;
   }
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: iterate_1( size_t ref_step,
                                    PDE_SetOfDiscreteFields* sdf,
                                    PDE_LocalFE* fe,
                                    std::ostream& os )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: iterate" ) ;

   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      fe->require_field_calculation( ff, PDE_LocalFE::N ) ;
   }
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      os << "------------------------------" << endl ;
      os << "current mesh" << endl ;
      fe->print_current_mesh( os, 3 ) ;
      os << "Local discretization of fields" << endl ;
      fe->print_local_discretization_of_fields( os, 3 ) ;
   }
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: process_domain( size_t ref_step,
                                         PDE_DomainAndFields const* dom,
                                         std::ostream& os )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: process_domain" ) ;

   os << "********* GRID *********************" << endl ;
   dom->print_grid( os, 0 ) ;

   PDE_SetOfDiscreteFields* sdf = dom->set_of_discrete_fields() ;

   os << "********* CELLS ********************" << endl ;
   PDE_LocalFEcell* cfe = dom->create_LocalFEcell( 0 ) ;
   cfe->include_color( GE_Color::halo_color() ) ;
   iterate( ref_step, sdf, cfe, os ) ;
   for( size_t i=0 ; i<FFS_VERTS.size() ; ++i )
   {
      check_values_at_vertices( ref_step,
                                sdf->item( FFS_VERTS( i ) ),
                                dom->set_of_vertices(),
                                cfe ) ;
   }
   if( VAL_CHECK != 0 ) check_interpolation( ref_step, cfe ) ;
   cfe->destroy() ; cfe=0 ;

   os << "********* BOUNDS ********************" << endl ;
   PDE_LocalFEbound* bfe = dom->create_LocalFEbound( 0 ) ;
   iterate( ref_step, sdf, bfe, os ) ;
   for( size_t i=0 ; i<FFS_VERTS.size() ; ++i )
   {
      check_values_at_vertices( ref_step,
                                sdf->item( FFS_VERTS( i ) ),
                                dom->set_of_vertices(),
                                bfe ) ;
   }
   if( VAL_BDIMP != 0 ) check_imposed_values( ref_step, bfe ) ;
   bfe->destroy() ; bfe=0 ;

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

   if( CHECK_BD_CELLS )
   {
      check_cell_and_faces( ref_step, dom ) ;
      check_covering_of_cells_boundary( ref_step, dom ) ;
   }
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: iterate( size_t ref_step,
                                  PDE_SetOfDiscreteFields* sdf,
                                  PDE_LocalFE* fe,
                                  std::ostream& os )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: iterate" ) ;

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
PDE_AdapterCHARMS_TEST:: check_node_location( size_t ref_step,
                                              PDE_SetOfDiscreteFields* sdf,
                                              PDE_LocalFE* fe,
                                              bool& ok_node_loc ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: check_node_location" ) ;

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
PDE_AdapterCHARMS_TEST:: check_interpolation( size_t ref_step,
                                              PDE_LocalFE* fe )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: check_interpolation" ) ;

   PEL_ASSERT( QRP != 0 ) ;

   bool ok_val = true ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      fe->start_IP_iterator( QRP ) ;
      for( ; fe->valid_IP() ; fe->go_next_IP() )
      {
         fill_context_with_coordinates( fe->coordinates_of_IP() ) ;
         for( size_t ic=0 ; ic<UU_CHECK->nb_components() ; ++ic )
         {
            double yy = VAL_CHECK->to_double_vector( CONTEXT )( ic ) ;
            double xx = fe->value_at_IP( UU_CHECK, 0 ) ;
            bool eq = PEL::double_equality( yy, xx, EPS_INTERP, MIN_INTERP ) ;
            if( !eq ) display_error( fe, yy, xx ) ;
            ok_val = ok_val && eq ;
         }
      }
   }

   ostringstream mesg ;
   mesg << TEST_NAME << "(" << ref_step << ") \""
        << UU_CHECK->name() << "\" interpolation" ;
   notify_one_test_result( mesg.str(), ok_val ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: check_imposed_values( size_t ref_step,
                                               PDE_LocalFEbound* bfe )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: check_interpolation" ) ;

   PEL_ASSERT( QRP != 0 ) ;

   bool ok_val = true ;
   bool done = false ;
   for( bfe->start() ; bfe->is_valid() ; bfe->go_next() )
   {
      if( COLOR_BDIMP->is_matching( bfe->color() ) )
      {
         done = true ;
         bfe->start_IP_iterator( QRP ) ;
         for( ; bfe->valid_IP() ; bfe->go_next_IP() )
         {
            fill_context_with_coordinates( bfe->coordinates_of_IP() ) ;
            for( size_t ic=0 ; ic<UU_BDIMP->nb_components() ; ++ic )
            {
               double yy = VAL_BDIMP->to_double_vector( CONTEXT )( ic ) ;
               double xx = bfe->value_at_IP( UU_BDIMP, 0 ) ;
               bool eq = PEL::double_equality( yy, xx, EPS_BDIMP, MIN_BDIMP ) ;
               if( !eq ) display_error( bfe, yy, xx ) ;
               ok_val = ok_val && eq ;
            }
         }
      }
   }

   ostringstream mesg ;
   mesg << TEST_NAME << "(" << ref_step << ") \""
        << UU_BDIMP->name() << "\" imposed values" ;
   notify_one_test_result( mesg.str(), ok_val && done ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: check_values_at_vertices(
                                              size_t ref_step,
                                              PDE_DiscreteField* ff,
                                              GE_SetOfPoints const* vertices,
                                              PDE_LocalFE* fe )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: check_values_at_vertices" ) ;

   doubleArray2D XX( vertices->nb_points(), ff->nb_components() ) ;
   XX.set( PEL::bad_double() ) ;

   size_t level = 0 ;

   // *******
   // etape 1 : on donne une valeur quelconque aux DOFs
   // *******
   if( ff != UU_CHECK )
   {
      for( size_t n=0 ; n<ff->nb_nodes() ; ++n )
      {
         for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
         {
            if( ff->node_is_active( n ) &&
                ff->DOF_value( level, n, ic ) == PEL::bad_double() )
            {
               // on affecte une valeur quelconque
               double xx = (double)n ;
               ff->set_DOF_value( level, n, xx, ic ) ;
            }
         }
      }
   }

   // *******
   // etape 2 : on vérifie que la valeur au sommet est bien la meme quelle que
   //           soit la maille dans laquelle on se place
   // *******
   bool ok = true ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      GE_Mpolyhedron const* poly = fe->polyhedron() ;
      for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
      {
         GE_Point const* pt = poly->vertex( i ) ;
         size_t iv = vertices->index( pt ) ;
         fe->set_calculation_point( pt ) ;
         for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
         {
            double val = fe->value_at_pt( ff, level, ic ) ;
            if( XX( iv, ic ) != PEL::bad_double() )
            {
               bool eq = PEL::double_equality( val, XX(iv,ic),
                                               EPS_VERTS, MIN_VERTS ) ;
               if( !eq )
               {
                  pt->print( cout, 3 ) ;
                  cout << " from " << fe->mesh_id() << " : " << val
                       << "  from previous : " << XX( iv, ic ) << endl ;
               }
               ok = ok && eq ;
            }
            XX( iv, ic ) = val ;
         }
      }
   }

   ostringstream mesg ;
   mesg << TEST_NAME << "(" << ref_step << ") \""
        << ff->name() << "\" at vertices" ;
   notify_one_test_result( mesg.str(), ok ) ;
}

//----------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: check_cell_and_faces( size_t ref_step,
                                               PDE_DomainAndFields const* dom )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: check_cell_and_faces" ) ;

   GE_SetOfPoints const* vertices = dom->set_of_vertices() ;
   
   // first way to determine the active vertices :
   //    vertices of active cells and vertices of active periodic faces
   size_t nb_active_1 = 0 ;
   size_t_vector active_idx_1( vertices->nb_points() ) ;
   active_idx_1.set( PEL::bad_index() ) ;
   
   bool ok_ref_level = true ;
   
   PDE_LocalFEcell* cfe = dom->create_LocalFEcell( 0 ) ;
   for( cfe->start() ; cfe->is_valid() ; cfe->go_next() )
   {
      GE_Mpolyhedron const* poly = cfe->polyhedron() ;
      PDE_ResultSaver::update_active_vertices( vertices, cfe->polyhedron(),
                                               active_idx_1, nb_active_1 ) ;
      size_t_vector const& c_ids = cfe->adjacent_cell_ids() ;
      size_t_vector const& s_ids = cfe->adjacent_side_ids() ;
      size_t_vector const& b_ids = cfe->adjacent_bound_ids() ;
      ok_ref_level &= ( c_ids.size() == s_ids.size() ) ;
      ok_ref_level &= ( s_ids.size() + b_ids.size() >= poly->nb_faces() ) ;
   }
   PDE_CursorFEside* sfe = dom->create_CursorFEside( 0 ) ;
   for( sfe->start() ; sfe->is_valid() ; sfe->go_next() )
   {
      if( sfe->is_periodic() )
      {
         PDE_ResultSaver::update_active_vertices( vertices, sfe->polyhedron(0),
                                                  active_idx_1, nb_active_1 ) ;
         PDE_ResultSaver::update_active_vertices( vertices, sfe->polyhedron(1),
                                                  active_idx_1, nb_active_1 ) ;
      }
   }

   // second way to determine the active vertices :
   //    vertices of active faces
   size_t nb_active_2 = 0 ;
   size_t_vector active_idx_2( vertices->nb_points() ) ;
   active_idx_2.set( PEL::bad_index() ) ;

   PDE_LocalFEbound* bfe = dom->create_LocalFEbound( 0 ) ;
   for( bfe->start() ; bfe->is_valid() ; bfe->go_next() )
   {
      PDE_ResultSaver::update_active_vertices( vertices, bfe->polyhedron(),
                                               active_idx_2, nb_active_2 ) ;
      cfe->go_i_th( bfe->adjacent_cell_id() ) ;
      bool ook = ( bfe->refinement_level() == cfe->refinement_level() ) ;
      ok_ref_level &= ook ;
      if( !ook )
      {
         PEL::out() << "********" << endl ;
         bfe->print_current_mesh( PEL::out(), 3 ) ;
         PEL::out() << "********" << endl ;
         cfe->print_current_mesh( PEL::out(), 3 ) ;
         PEL_Error::object()->raise_plain( "error in refinement levels" ) ;
      }
   }
   bfe->destroy() ; bfe = 0 ;
   cfe->destroy() ; cfe = 0 ;
   for( sfe->start() ; sfe->is_valid() ; sfe->go_next() )
   {
      if( sfe->is_periodic() )
      {
         PDE_ResultSaver::update_active_vertices( vertices, sfe->polyhedron(0),
                                                  active_idx_2, nb_active_2 ) ;
         PDE_ResultSaver::update_active_vertices( vertices, sfe->polyhedron(1),
                                                  active_idx_2, nb_active_2 ) ;
      }
      else
      {
         PDE_ResultSaver::update_active_vertices( vertices, sfe->polyhedron(),
                                                  active_idx_2, nb_active_2 ) ;
         
      }

      size_t s_level = sfe->refinement_level() ;
      PDE_LocalFEcell const* cfe_0 = sfe->adjacent_localFEcell( 0 ) ;
      size_t c0_level = cfe_0->refinement_level() ;
      PDE_LocalFEcell const* cfe_1 = sfe->adjacent_localFEcell( 1 ) ;
      size_t c1_level = cfe_1->refinement_level() ;
      bool ook = ( ( (s_level == c0_level) && (s_level >= c1_level) ) ||
                   ( (s_level == c1_level) && (s_level >= c0_level) ) ) ;
      ok_ref_level &= ook ;
      if( !ook )
      {
         PEL::out() << "********" << endl ;
         sfe->print_current_mesh( PEL::out(), 3 ) ;
         PEL::out() << "********" << endl ;
         cfe_0->print_current_mesh( PEL::out(), 3 ) ;
         PEL::out() << "********" << endl ;
         cfe_1->print_current_mesh( PEL::out(), 3 ) ;
         PEL_Error::object()->raise_plain( "error in refinement levels" ) ;
      }
   }
   sfe->destroy() ; sfe = 0 ;

   ostringstream mesg1 ;
   mesg1 << TEST_NAME 
         << "(" << ref_step << ") consistency between cells and faces" ;
   notify_one_test_result( mesg1.str(), ok_ref_level ) ;
   
   bool ok = ( nb_active_1 == nb_active_2 ) ;
   for( size_t i=0 ; i<vertices->nb_points() ; ++i )
   {
      bool ook = ( ( active_idx_1( i ) == PEL::bad_index() ) && 
                   ( active_idx_2( i ) == PEL::bad_index() ) )
                 ||
                 ( ( active_idx_1( i ) != PEL::bad_index() ) && 
                   ( active_idx_2( i ) != PEL::bad_index() ) ) ;
      ok &= ook ;
      if( !ook )
      {
         GE_Point const* pt = vertices->point( i ) ;
         pt->print( PEL::out(), 0 ) ; PEL::out() << endl ;
         PEL::out() << "   active_idx_1 : " << active_idx_1( i )
                    << endl ;
         PEL::out() << "   active_idx_2 : " << active_idx_2( i )
                    << endl ;
      }
   }

   ostringstream mesg2 ;
   mesg2 << TEST_NAME << "(" << ref_step << ") vertices from cells and faces" ;
   notify_one_test_result( mesg2.str(), ok ) ;
}

//---------------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: check_covering_of_cells_boundary(
                                            size_t ref_step,
                                            PDE_DomainAndFields const* dom  )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: check_covering_of_cells_boundary" ) ;

   if( dom->nb_space_dimensions() != 2 ) return ;

   bool ok = true ;

   size_t dim = 2 ;
   size_t nb_points = 30 ;

   GE_SegmentPolyhedron_INT* sp =
    GE_SegmentPolyhedron_INT::make( 0, "GE_SegmentPolyhedron1D_INT", dim, 0 ) ;

   GE_Point* pt = GE_Point::create( 0, dim ) ;
   PDE_LocalFEcell*  c_fe = dom->create_LocalFEcell( 0 ) ;
   PDE_LocalFEbound* b_fe = dom->create_LocalFEbound( 0 ) ;
   PDE_CursorFEside* s_fe = dom->create_CursorFEside( 0 ) ;
   for( c_fe->start() ; c_fe->is_valid() ; c_fe->go_next() )
   {
      size_t_vector const& sid = c_fe->adjacent_side_ids() ;
      size_t_vector const& bid = c_fe->adjacent_bound_ids() ;

      GE_Mpolyhedron const* poly = c_fe->polyhedron() ;
      GE_Point const* center = poly->center() ;
      double x_c = center->coordinate( 0 ) ;
      double y_c = center->coordinate( 1 ) ;
      double rr = poly->inter_vertices_maximum_distance() ;
      for( size_t i=0 ; i<nb_points ; ++i )
      {
         double theta = 2.0*PEL::pi()*((double)i)/((double)nb_points) ;
         pt->set_coordinate( 0, x_c + rr*PEL::cos( theta ) ) ;
         pt->set_coordinate( 1, y_c + rr*PEL::sin( theta ) ) ;
         size_t nb_inter = 0 ;
         for( size_t j=0 ; j<sid.size() ; ++j )
         {
            s_fe->go_i_th( sid( j ) ) ;
            GE_Mpolyhedron const* spoly = 0 ;
            if( s_fe->is_periodic() )
            {
               if( s_fe->adjacent_localFEcell( 0 )->mesh_id() == 
                   c_fe->mesh_id() )
               {
                  spoly = s_fe->polyhedron( 0 ) ;
                  PEL_ASSERT( s_fe->adjacent_localFEcell( 1 )->mesh_id() != 
                              c_fe->mesh_id() ) ;
               }
               else
               {
                  PEL_ASSERT( s_fe->adjacent_localFEcell( 1 )->mesh_id() == 
                              c_fe->mesh_id() ) ;
                  spoly = s_fe->polyhedron( 1 ) ;
               }
            }
            else
            {
               spoly = s_fe->polyhedron() ;
            }
            sp->check_intersection( center, pt, spoly ) ;
            if( sp->one_single_intersection() )
            {
               nb_inter++ ;
            }
         }
         for( size_t j=0 ; j<bid.size() ; ++j )
         {
            b_fe->go_i_th( bid( j ) ) ;
            sp->check_intersection( center, pt, b_fe->polyhedron() ) ;
            if( sp->one_single_intersection() )
            {
               nb_inter++ ;
            }
         }
         if( nb_inter!=1 && nb_inter!=2 )
         {
            PEL::out() << "******** nb_inter = " << nb_inter << endl ;
            c_fe->print( PEL::out(), 3 ) ;
            // poly->print( PEL::out(), 3 ) ;
            PEL::out() << "segment : " ;
            center->print( PEL::out(), 0 ) ; PEL::out() << " -> " ;
            pt->print( PEL::out(), 0 ) ; PEL::out() << endl ;
            for( size_t j=0 ; j<sid.size() ; ++j )
            {
               s_fe->go_i_th( sid( j ) ) ;
               GE_Mpolyhedron const* spoly = ( s_fe->is_periodic() ? 
                              s_fe->polyhedron( 0 ) : s_fe->polyhedron() ) ;
               spoly->print( PEL::out(), 6 ) ;
            }
            for( size_t j=0 ; j<bid.size() ; ++j )
            {
               b_fe->go_i_th( bid( j ) ) ;
               b_fe->polyhedron()->print( PEL::out(), 6 ) ;
            }
            ok = false ;
            break ;
         }
      }
   }
   c_fe->destroy() ;
   s_fe->destroy() ;
   b_fe->destroy() ;
   pt->destroy() ;
   sp->destroy() ;

   ostringstream mesg ;
   mesg << TEST_NAME << "(" << ref_step << ") closure of cells" ;
   notify_one_test_result( mesg.str(), ok ) ;
}

//----------------------------------------------------------------------
void
PDE_AdapterCHARMS_TEST:: fill_context_with_coordinates( GE_Point const* pt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdapterCHARMS_TEST:: fill_context_with_coordinates" ) ;

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
PDE_AdapterCHARMS_TEST:: display_error( PDE_LocalFE const* fe,
                                        double theo, double xx ) const
//-----------------------------------------------------------------------
{
   std::cout << fe->type_name() << " mesh : " << fe->mesh_id()
             << " theo " << theo << " xx " << xx << std::endl ;
}
