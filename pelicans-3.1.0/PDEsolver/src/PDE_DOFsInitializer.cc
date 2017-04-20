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

#include <PDE_DOFsInitializer.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <boolArray2D.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>
#include <GE_SetOfPoints.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_InterfaceAndFields.hh>
#include <PDE_LocalFEmortarSide.hh>
#include <PDE_ValuesOnMeshingProjector.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultReader.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_ValuesOnMeshing.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;

struct PDE_DOFsInitializer_ERROR
{
   static void n0( std::string const& name, size_t nbc, size_t given_nbc ) ;
   static void n1( std::string const& name, size_t ic, GE_Point const* pt,
                   double val1, double val2 ) ;
   static void n3( std::string const& name, size_t ic ) ;
   static void n4( std::string const& name ) ;
   static void n5( PEL_ModuleExplorer const* exp, size_t nb_comps ) ;
   static void n6( std::string const& name, size_t nbc, size_t given_nbc ) ;
   static void n7( std::string const& name ) ;
   static void n8( std::string const& name ) ;
} ;

//-------------------------------------------------------------------------
PDE_DOFsInitializer const*
PDE_DOFsInitializer:: object( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: object" ) ;

   static PDE_DOFsInitializer const* result = new PDE_DOFsInitializer() ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==PEL_Root::object() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_DOFsInitializer:: PDE_DOFsInitializer( void )
//-------------------------------------------------------------------------
   : PEL_Object( PEL_Root::object() )
   , CTX( 0 )
   , COORDS( 0 )
{
   // Context for indicator :
   {
      PEL_ContextSimple* c = PEL_ContextSimple::create( this ) ;
      COORDS = PEL_DoubleVector::create( c, doubleVector(0) ) ;
      c->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
      CTX = c ;
   }
}

//-------------------------------------------------------------------------
PDE_DOFsInitializer:: ~PDE_DOFsInitializer( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: initialize_discrete_fields(
                                PDE_DomainAndFields const* df,
                                PEL_ModuleExplorer const* df_exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: initialize_discrete_fields" ) ;
   PEL_CHECK_PRE( df!=0 ) ;
   PEL_CHECK_PRE( df_exp!=0 ) ;
   PEL_CHECK_PRE( df_exp->has_entry( "verbose_level" ) ) ;
   PEL_CHECK_PRE( df_exp->has_module( "interior_fields" ) ) ;

   bool verbose = df_exp->int_data( "verbose_level" ) ;

   set_discrete_fields_imposed_DOFs( df, df_exp, verbose ) ;

   set_discrete_fields_DOFs_values( df, df_exp, verbose ) ;
}
   
//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: initialize_discrete_fields(
                                PDE_InterfaceAndFields const* interf,
                                PEL_ModuleExplorer const* interf_exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: initialize_discrete_fields" ) ;
   PEL_CHECK_PRE( interf!=0 ) ;
   PEL_CHECK_PRE( interf_exp!=0 ) ;
   PEL_CHECK_PRE( interf_exp->has_entry( "verbose_level" ) ) ;
   PEL_CHECK_PRE( interf_exp->has_module( "fields" ) ) ;

   bool verbose = interf_exp->int_data( "verbose_level" ) ;

   set_discrete_fields_imposed_DOFs( interf, interf_exp, verbose ) ;
}
   
//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: set_discrete_fields_DOFs_values(
                                PDE_DomainAndFields const* df,
                                PEL_ModuleExplorer const* df_exp,
                                bool verbose ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: initialize_discrete_fields" ) ;
   PEL_CHECK_PRE( df!=0 ) ;
   PEL_CHECK_PRE( df_exp!=0 ) ;
   PEL_CHECK_PRE( df_exp->has_module( "interior_fields" ) ) ;
   
   PEL_Vector* readers = PEL_Vector::create( 0, 0 ) ;
   if( df_exp->has_module( "PDE_ResultReader" ) &&
       df_exp->has_module( "set_of_PDE_ResultReaders" ) )
   {
      PEL_Error::object()->raise_plain(
         "*** Definition of PDE_ResultReader error :\n"
         "    modules PDE_ResultReader and set_of_PDE_ResultReaders" ) ;
   }
   if( df_exp->has_module( "PDE_ResultReader" ) )
   {
      PEL_ModuleExplorer const* reader_exp =
                  df_exp->create_subexplorer( 0, "PDE_ResultReader" ) ;
      PDE_ResultReader* reader =
                      PDE_ResultReader::create( readers, reader_exp ) ;
      readers->append( reader ) ;
      reader_exp->destroy() ; reader_exp = 0 ;
   }
   if( df_exp->has_module( "set_of_PDE_ResultReaders" ) )
   {
      PEL_ModuleExplorer* ee =
           df_exp->create_subexplorer( 0, "set_of_PDE_ResultReaders" ) ;
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer const* reader_exp = ee->create_subexplorer(0) ;
         PDE_ResultReader* reader =
                      PDE_ResultReader::create( readers, reader_exp ) ;
         readers->append( reader ) ;
         reader_exp->destroy() ; reader_exp = 0 ;
      }
      ee->destroy() ; ee = 0 ;
   }

   if( df_exp->has_module( "L2_projection" ) )
   {
      // For data deck checker...
   }

   if( verbose )
   {
      PEL::out() << "*** Initializing fields" << std::endl ;
   }

   PEL_ModuleExplorer* iexp =
                           df_exp->create_subexplorer( 0, "interior_fields" ) ;
   iexp->start_module_iterator() ;
   for( ; iexp->is_valid_module() ;  iexp->go_next_module() )
   {
      PEL_ModuleExplorer* fexp = iexp->create_subexplorer( iexp ) ;
      std::string const& fname = fexp->string_data( "name" ) ;
      PDE_DiscreteField* f = df->set_of_discrete_fields()->item( fname ) ;
      apply_DOFs_values( f, df, df_exp, fexp, readers, false, verbose ) ;   
   }
   iexp->destroy() ; iexp=0 ;

   if( df_exp->has_module( "boundary_fields" ) )
   {
      iexp = df_exp->create_subexplorer( 0, "boundary_fields" ) ;
      iexp->start_module_iterator() ;
      for( ; iexp->is_valid_module() ; iexp->go_next_module() )
      {
         PEL_ModuleExplorer* fexp = iexp->create_subexplorer( iexp ) ;
         std::string const& fname = fexp->string_data( "name" ) ;
         PDE_DiscreteField* f = df->set_of_discrete_fields()->item( fname ) ;
         apply_DOFs_values( f, df, df_exp, fexp, readers, true, verbose ) ;
      }
      iexp->destroy() ; iexp=0 ;
   }
   
   if( verbose )
   {
      PEL::out() << std::endl ;
   }

   readers->destroy() ; readers = 0 ;
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: set_discrete_fields_imposed_DOFs(
                                PDE_DomainAndFields const* df,
                                PEL_ModuleExplorer const* df_exp,
                                bool verbose ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: set_discrete_fields_imposed_DOFs" ) ;
   PEL_CHECK_PRE( df!=0 ) ;
   PEL_CHECK_PRE( df_exp!=0 ) ;
   PEL_CHECK_PRE( df_exp->has_module( "interior_fields" ) ) ;

   if( verbose )
   {
      PEL::out() << "*** Taking into account requests of "
                 << "\"DOFs_imposed_value\" modules" << std::endl ;
   }

   PEL_ModuleExplorer* iexp =
                           df_exp->create_subexplorer( 0, "interior_fields" ) ;
   iexp->start_module_iterator() ;
   for( ; iexp->is_valid_module() ;  iexp->go_next_module() )
   {
      PEL_ModuleExplorer* fexp = iexp->create_subexplorer( iexp ) ;
      std::string const& fname = fexp->string_data( "name" ) ;
      PDE_DiscreteField* f = df->set_of_discrete_fields()->item( fname ) ;
      apply_DOFs_imposed_value( f, df, fexp, false, verbose ) ;   
   }
   iexp->destroy() ; iexp=0 ;

   if( df_exp->has_module( "boundary_fields" ) )
   {
      iexp = df_exp->create_subexplorer( 0, "boundary_fields" ) ;
      iexp->start_module_iterator() ;
      for( ; iexp->is_valid_module() ; iexp->go_next_module() )
      {
         PEL_ModuleExplorer* fexp = iexp->create_subexplorer( iexp ) ;
         std::string const& fname = fexp->string_data( "name" ) ;
         PDE_DiscreteField* f = df->set_of_discrete_fields()->item( fname ) ;
         apply_DOFs_imposed_value( f, df, fexp, true, verbose ) ;
      }
      iexp->destroy() ; iexp=0 ;
   }
   
   if( verbose )
   {
      PEL::out() << std::endl ;
   }
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: set_discrete_fields_imposed_DOFs(
                                PDE_InterfaceAndFields const* interf,
                                PEL_ModuleExplorer const* interf_exp,
                                bool verbose ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: set_discrete_fields_imposed_DOFs" ) ;
   PEL_CHECK_PRE( interf!=0 ) ;
   PEL_CHECK_PRE( interf_exp!=0 ) ;
   PEL_CHECK_PRE( interf_exp->has_module( "fields" ) ) ;

   if( verbose )
   {
      PEL::out() << "*** Taking into account requests of "
                 << "\"DOFs_imposed_value\" modules" << std::endl ;
   }

   PEL_ModuleExplorer* exp =
                       interf_exp->create_subexplorer( 0, "fields" ) ;
   exp->start_module_iterator() ;
   for( ; exp->is_valid_module() ;  exp->go_next_module() )
   {
      PEL_ModuleExplorer* fexp = exp->create_subexplorer( exp ) ;
      std::string const& fname = fexp->string_data( "name" ) ;
      PDE_DiscreteField* f = interf->set_of_discrete_fields()->item( fname ) ;
      apply_DOFs_imposed_value( f, interf, fexp, verbose ) ;   
   }
   exp->destroy() ; exp=0 ;
   
   if( verbose )
   {
      PEL::out() << std::endl ;
   }
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: apply_DOFs_values(
                                      PDE_DiscreteField* f,
                                      PDE_DomainAndFields const* df,
                                      PEL_ModuleExplorer const* df_exp,
                                      PEL_ModuleExplorer const* f_exp,
                                      PEL_Vector const* readers,
                                      bool is_boundary_field,
                                      bool verbose ) const 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: apply_DOFs_values" ) ;
   PEL_CHECK( f!=0 ) ;
   PEL_CHECK( df!=0 ) ;
   PEL_CHECK( df_exp!=0 ) ;
   PEL_CHECK( f_exp!=0 ) ;
   PEL_CHECK( readers != 0 ) ;
   PEL_CHECK(
      FORALL(
         ( size_t i=0 ; i<readers->index_limit() ; ++i ),
         readers->at(i) != 0 &&
         dynamic_cast<PDE_ResultReader const*>( readers->at(i) ) != 0 ) ) ;

   PEL_ModuleExplorer* sexp = f_exp->create_subexplorer( 0, "DOFs_values" ) ;

   std::string mode = "set_by_value_at_node_location" ;
   if( sexp->has_entry( "mode" ) )
   {
      mode = sexp->string_data( "mode" ) ;
      if( mode!="set_by_value_at_node_location" &&
          mode!="set_by_L2_projection" )
      {
         std::string const val =
            "   - \"set_by_value_at_node_location\"\n"
            "   - \"set_by_L2_projection\"" ;
         PEL_Error::object()->raise_bad_data_value( sexp, "mode", val ) ;
      }
   }
   PDE_ResultReader const* reader = 0 ;
   for( size_t i=0 ; i<readers->index_limit() ; ++i )
   {
      PDE_ResultReader const* r =
                       static_cast<PDE_ResultReader const*>( readers->at(i) ) ;
      if( r->has_field( f->name() ) )
      {
         if( reader != 0 )
         {
            PEL_Error::object()->raise_plain(
               "*** Definition of PDE_ResultReader error :\n"
               "    discrete field \""+f->name()+"\" is defined twice" ) ;
         }
         reader = r ;
      }
   }
   if( mode=="set_by_L2_projection" )
   {
      initialize_by_L2_projection( f, df, df_exp, sexp, reader,
                                   is_boundary_field, verbose ) ;
   }
   else if( mode=="set_by_value_at_node_location" )
   {
      initialize_by_value_at_node_location( f, df, df_exp, sexp, reader,
                                            is_boundary_field, verbose ) ;
   }
   sexp->destroy() ; sexp = 0 ;
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: apply_DOFs_imposed_value( 
                                         PDE_DiscreteField* f,
                                         PDE_DomainAndFields const* df,
                                         PEL_ModuleExplorer const* f_exp,
                                         bool is_boundary_field,
                                         bool verbose ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: apply_DOFs_imposed_value" ) ;
   
   boolArray2D already_done( f->nb_nodes(), f->nb_components(), false ) ;

   if( f_exp->has_module( "DOFs_imposed_value" ) )
   {
      if( is_boundary_field )
      {
         PDE_DOFsInitializer_ERROR::n4( f->name() ) ;
      }
      PEL_ModuleExplorer* sexp =
                f_exp->create_subexplorer( 0, "DOFs_imposed_value" ) ;
      sexp->start_module_iterator() ;
      for( ; sexp->is_valid_module() ; sexp->go_next_module() )
      {
         PEL_ModuleExplorer* sse = sexp->create_subexplorer( 0 ) ;
         std::string const& location = sse->string_data( "location" ) ;
         if( location == "on_bounds" )
         {
            impose_on_bounds( f, df, sse, verbose, already_done ) ;
         }
         else if( location == "in_region" )
         {
            PDE_LocalFEcell* fe = df->create_LocalFEcell( 0 ) ;
            impose_in_region( f, df->set_of_vertices(), fe, sse, verbose, 
                              already_done ) ;
            fe->destroy() ; fe = 0 ;
         }
         else
         {
            PEL_Error::object()->raise_bad_data_value(
               sse, "location", "   - \"on_bounds\"\n   - \"in_region\"\n" ) ;
         }
         sse->destroy() ;
      }
      sexp->destroy() ; sexp = 0 ;
   }

   if( verbose )
   {
      size_t nb_imposed = 0 ;
      for( size_t n=0 ; n<f->nb_nodes() ; ++n )
         for( size_t ic=0 ; ic< f->nb_components() ; ++ic )
            if( f->DOF_has_imposed_value( n, ic ) ) ++nb_imposed ;
      PEL::out() << "    field \"" << f->name()
                 << "\" : " << nb_imposed << " imposed DOFs out of " 
                 << f->nb_nodes()*f->nb_components() << std::endl ;
   }
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: apply_DOFs_imposed_value( 
                                         PDE_DiscreteField* f,
                                         PDE_InterfaceAndFields const* interf,
                                         PEL_ModuleExplorer const* f_exp,
                                         bool verbose ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: apply_DOFs_imposed_value" ) ;

   boolArray2D already_done( f->nb_nodes(), f->nb_components(), false ) ;
   
   if( f_exp->has_module( "DOFs_imposed_value" ) )
   {
      PEL_ModuleExplorer* sexp =
                f_exp->create_subexplorer( 0, "DOFs_imposed_value" ) ;
      sexp->start_module_iterator() ;
      for( ; sexp->is_valid_module() ; sexp->go_next_module() )
      {
         PEL_ModuleExplorer* sse = sexp->create_subexplorer( 0 ) ;
         std::string const& location = sse->string_data( "location" ) ;
         if( location == "in_region" )
         {
            PDE_LocalFEmortarSide* fe = interf->create_LocalFEmortarSide( 0 ) ;
            impose_in_region( f, interf->set_of_vertices(), fe, sse, verbose,
                              already_done ) ;
            fe->destroy() ; fe = 0 ;
         }
         else
         {
            PEL_Error::object()->raise_bad_data_value(
               sse, "location", "   - \"in_region\"\n" ) ;
         }
         sse->destroy() ;
      }
      sexp->destroy() ; sexp = 0 ;
   }

   if( verbose )
   {
      size_t nb_imposed = 0 ;
      for( size_t n=0 ; n<f->nb_nodes() ; ++n )
         for( size_t ic=0 ; ic< f->nb_components() ; ++ic )
            if( f->DOF_has_imposed_value( n, ic ) ) ++nb_imposed ;
      PEL::out() << "    field \"" << f->name()
                 << "\" : " << nb_imposed << " imposed DOFs out of " 
                 << f->nb_nodes()*f->nb_components() << std::endl ;
   }
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: impose_on_bounds( PDE_DiscreteField* f,
                                        PDE_DomainAndFields const* df,
                                        PEL_ModuleExplorer const* exp,
                                        bool verbose,
                                        boolArray2D& already_done ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: impose_on_bounds" ) ;

   GE_Color const* color = 0 ;
   if( exp->has_entry( "color" ) )
   {
      color = GE_Color::object( exp->string_data( "color" ) ) ;
   }

   size_t nbc = PEL::bad_index() ;
   int comp = PEL::bad_int() ;
   read_components( exp, f, nbc, comp ) ;
   doubleVector val( nbc ) ;

   PDE_ValuesOnMeshing* cv =
      PDE_ValuesOnMeshing::create( 0,
                                   f->name(), nbc,
                                   exp,
                                   df->set_of_vertices(),
                                   0 ) ;
   
   PDE_LocalFEbound* fe = df->create_LocalFEbound( 0 ) ;
   fe->include_color( GE_Color::halo_color() ) ;
   fe->require_field_calculation( f, PDE_LocalFE::N ) ;

   size_t refi_level = 0 ;
   size_t refi_level_max = 0 ;
   do
   {
      for( fe->start() ; fe->is_valid() ; fe->go_next() )
      {
         if( !( ( fe->color() != GE_Color::halo_color() ) &&
                ( color==0 || color->is_matching( fe->color() ) ) ) ) continue ;
         
         cv->set_mesh( fe->polyhedron(), fe->color() ) ;
         for( size_t in = 0 ; in<fe->nb_local_nodes( f ) ; ++in )
         {
            size_t rl = fe->node_refinement_level( f, in ) ;
            if( rl > refi_level_max ) refi_level_max = rl ;
            if( rl == refi_level && fe->local_node_is_in_mesh( f, in ) )
            { 
               GE_Point const* pt_node = fe->local_node_location( f, in ) ;
               cv->compute_value( pt_node, val ) ;
               if( val.size() != nbc )
               {
                  PDE_DOFsInitializer_ERROR::n6( f->name(), nbc, val.size() ) ;
               }
               
               if( refi_level != 0 )
               {
                  fe->set_calculation_point( pt_node ) ;
                  for( size_t jn=0 ; jn<fe->nb_local_nodes( f ) ; ++jn )
                  {
                     size_t ll = fe->node_refinement_level( f, jn ) ;
                     if( ll < refi_level )
                     {
                        size_t jgn = fe->global_node( f, jn ) ;
                        double Npt = fe->N_at_pt( f, jn ) ;
                        if( PEL::abs( Npt ) > 1.e-5 ) //?????? 1.e-5
                        {
                           size_t ic = ( comp==-1 ) ? 0 : (size_t) comp ;
                           while( ic<f->nb_components() )
                           {
                              //???? mettre msg erreur circonstancié
                              PEL_ASSERT( f->DOF_has_imposed_value( jgn, ic )) ;
                              double& x = ( comp==-1 ) ? val( ic ) : val( 0 ) ;
                              x -= f->DOF_imposed_value( jgn, ic ) * Npt ;
                              ic = ( comp==-1 ) ? ic+1 : f->nb_components() ;
                           }
                        }
                     }
                  }
               }

               size_t node = fe->global_node( f, in ) ;
               modify_field_DOFs( f, node, pt_node, comp, val, already_done ) ;
            }
         }
      }
      ++refi_level ;
   }
   while( refi_level <= refi_level_max ) ;
   
   cv->destroy() ; cv = 0 ;
   fe->destroy() ; fe = 0 ;
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: impose_in_region( PDE_DiscreteField* f,
                                        GE_SetOfPoints const* vertices,
                                        PDE_LocalFE* fe,
                                        PEL_ModuleExplorer const* exp,
                                        bool verbose,
                                        boolArray2D& already_done ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: impose_in_region" ) ;
   
   size_t nbc = PEL::bad_index() ;
   int comp = PEL::bad_int() ;
   read_components( exp, f, nbc, comp ) ;
   doubleVector val( nbc ) ;

   PDE_ValuesOnMeshing* cv =
      PDE_ValuesOnMeshing::create( 0,
                                   f->name(), nbc,
                                   exp,
                                   vertices,
                                   0 ) ;

   if( ! exp->has_entry( "indicator" ) )
   {
      PDE_DOFsInitializer_ERROR::n7( f->name() ) ;
   }
   
   if( fe->is_excluded( GE_Color::halo_color() ) )
   {
      fe->include_color( GE_Color::halo_color() ) ;
   }
   fe->require_field_calculation( f, PDE_LocalFE::N ) ;

   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      bool ok_cv = false ;
      for( size_t in = 0 ; in<fe->nb_local_nodes( f ) ; in++ )
      {
         if( fe->node_refinement_level( f, in ) != 0 )
         {
            PDE_DOFsInitializer_ERROR::n8( f->name() ) ;
         }
         size_t node = fe->global_node( f, in ) ;
         GE_Point const* pt_node = fe->local_node_location( f, in ) ;
         COORDS->set( pt_node->coordinate_vector() ) ;
         if( exp->bool_data( "indicator", CTX ) )
         {
            if( !ok_cv )
            { 
               cv->set_mesh( fe->polyhedron(), fe->color() ) ;
               ok_cv = true ;
            }
            cv->compute_value( pt_node, val ) ;
            if( val.size() != nbc )
            {
               PDE_DOFsInitializer_ERROR::n6( f->name(), nbc, val.size() ) ;
            }

            modify_field_DOFs( f, node, pt_node, comp, val, already_done ) ;
         }
      }
   }
   cv->destroy() ; cv = 0 ;
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: modify_field_DOFs( PDE_DiscreteField* ff,
                                         size_t node,
                                         GE_Point const* pt_node,
                                         int comp,
                                         doubleVector const& val,
                                         boolArray2D& already_done )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: modify_field_DOFs" ) ;

   static double const dbl_eps = 1.E-6 ;
   static double const dbl_min = 1.E-6 ;
   
   size_t ic = ( comp==-1 ) ? 0 : (size_t) comp ;
   while( ic<ff->nb_components() )
   {
      double x = ( comp==-1 ) ? val( ic ) : val( 0 ) ;
      
      if( already_done( node, ic ) )
      {
         double xold = ff->DOF_imposed_value( node, ic ) ;
         if( !PEL::double_equality( x, xold, dbl_eps, dbl_min ) )
            PDE_DOFsInitializer_ERROR::n1( ff->name(), ic, pt_node, x, xold ) ;
      }
      else
      {
         already_done( node, ic ) = true ;
      }

      ff->modify_DOF( true, x, node, ic ) ;
      
      ic = ( comp==-1 ) ? ic+1 : ff->nb_components() ;
   }
}
   
//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: read_components( PEL_ModuleExplorer const* exp,
                                       PDE_DiscreteField const* ff,
                                       size_t& nbc,
                                       int& comp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: read_component" ) ;

   comp = -1 ;
   nbc = ff->nb_components() ;
   if( exp->has_entry( "component" ) ) 
   {
      comp = exp->int_data( "component" ) ;
      nbc = 1 ;
      if( comp<0 || comp>=(int) ff->nb_components() )
      {
         PDE_DOFsInitializer_ERROR::n5( exp, ff->nb_components() ) ;
      }
   }
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: initialize_by_L2_projection(
                                        PDE_DiscreteField* f,
                                        PDE_DomainAndFields const* df,
                                        PEL_ModuleExplorer const* df_exp,
                                        PEL_ModuleExplorer const* exp,
                                        PDE_ResultReader const* reader,
                                        bool is_boundary_field,
                                        bool verbose ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: initialize_by_L2_projection" ) ;
   PEL_CHECK( f!=0 ) ;
   PEL_CHECK( exp!=0 ) ;

   if( is_boundary_field )
   {
      PEL_Error::object()->raise_plain(
         "Initialization by L2 projection is only\n"
         "available with interior fields" ) ;
   }

   if( verbose )
   {
      PEL::out() << "    projecting field \""
                 << f->name() << "\"" << std::endl ;
   }

   PDE_ValuesOnMeshing* cv = 
               PDE_ValuesOnMeshing::create( 0, f->name(), f->nb_components(),
                                            exp,
                                            df->set_of_vertices(), reader ) ;
   
   PEL_ModuleExplorer* ee = df_exp->create_subexplorer( 0, "L2_projection" ) ;
   PDE_ProjectorForDOFsSetting* proj =
                 PDE_ValuesOnMeshingProjector::create( 0, f, 0, cv, df, ee ) ;
   
   proj->project_and_update_field() ;
   for( size_t i_level=1 ; i_level<f->storage_depth() ; ++i_level )
   {
      f->copy_DOFs_value( 0, i_level ) ;
   }

   cv->destroy() ; cv = 0 ;
   proj->destroy() ; proj = 0 ;
   ee->destroy() ; ee = 0 ;
}

//-------------------------------------------------------------------------
void
PDE_DOFsInitializer:: initialize_by_value_at_node_location(
                                        PDE_DiscreteField* f,
                                        PDE_DomainAndFields const* df,
                                        PEL_ModuleExplorer const* df_exp,
                                        PEL_ModuleExplorer const* exp,
                                        PDE_ResultReader const* reader,
                                        bool is_boundary_field,
                                        bool verbose ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFsInitializer:: initialize_by_value_at_node_location" ) ;

   if( verbose )
   {
      PEL::out() << "    field \""
                 << f->name()
                 << "\" : computing values at node locations"
                 << std::endl ;
   }

   double const dbl_eps = 1.E-6 ;
   double const dbl_min = 1.E-6 ;
   
   size_t nbc = f->nb_components() ;
   doubleVector val( nbc ) ;
   
   PDE_ValuesOnMeshing* cv =
      PDE_ValuesOnMeshing::create( 0,
                                   f->name(), nbc,
                                   exp,
                                   df->set_of_vertices(),
                                   reader ) ;

   doubleArray2D values( f->nb_components(), f->nb_nodes() ) ;
   values.set( PEL::bad_double() ) ;

   PDE_LocalFE* fe = 0 ;
   if( is_boundary_field )
   {
      fe = df->create_LocalFEbound( 0 ) ;
   }
   else
   {
      fe = df->create_LocalFEcell( 0 ) ;
   }
   fe->include_color( GE_Color::halo_color() ) ;

   fe->require_field_calculation( f, PDE_LocalFE::N ) ;
   
   size_t refi_level = 0 ;
   size_t refi_level_max = 0 ;
   do
   {
      for( fe->start() ; fe->is_valid() ; fe->go_next() )
      {
         cv->set_mesh( fe->polyhedron(), fe->color() ) ;
         for( size_t in = 0 ; in<fe->nb_local_nodes( f ) ; in++ )
         {
            size_t rl = fe->node_refinement_level( f, in ) ;
            if( rl > refi_level_max ) refi_level_max = rl ;
            if( rl == refi_level && fe->local_node_is_in_mesh( f, in ) )
            {
               GE_Point const* pt_node = fe->local_node_location( f, in ) ;
               cv->compute_value( pt_node, val ) ;
               if( val.size() != f->nb_components() )
               {
                  PDE_DOFsInitializer_ERROR::n0( f->name(), nbc, val.size() );
               }
               
               if( refi_level != 0 )
               {
                  fe->set_calculation_point( pt_node ) ;
                  for( size_t jn=0 ; jn<fe->nb_local_nodes( f ) ; ++jn )
                  {
                     size_t ll = fe->node_refinement_level( f, jn ) ;
                     if( ll < refi_level )
                     {
                        size_t jgn = fe->global_node( f, jn ) ; 
                        double Npt = fe->N_at_pt( f, jn ) ;
                        for( size_t ic=0 ; ic<nbc ; ++ic )
                        {
                           PEL_ASSERT( f->node_is_active( jgn ) ) ;
                           PEL_ASSERT( values( ic, jgn ) != PEL::bad_double() ) ;
                           val( ic ) -= values( ic, jgn ) * Npt ;
                        }
                     }
                  }
               }
               size_t node = fe->global_node( f, in ) ;
               for( size_t ic = 0 ; ic<val.size() ; ic++ )
               {
                  if( values( ic, node )==PEL::bad_double() )
                  {
                     values( ic, node ) = val( ic ) ;
                  }
                  else if( !PEL::double_equality( values( ic, node ), val( ic ),
                                         dbl_eps, dbl_min ) )
                  {
                     PDE_DOFsInitializer_ERROR::n1( f->name(), ic,
                                                    pt_node,
                                                    values( ic, node ),
                                                    val( ic ) ) ;
                  }
               }
            }
         }
      }
      ++refi_level ;
   } 
   while( refi_level <= refi_level_max ) ;

   for( size_t i=0 ; i<f->nb_nodes() ; ++i )
   {
      if( f->node_is_active( i ) )
      {
         for( size_t ic=0 ; ic<f->nb_components() ; ++ic )
         {
            if( values( ic, i )==PEL::bad_double() )
            {
               PDE_DOFsInitializer_ERROR::n3( f->name(), ic )  ;
            }
            if( !f->DOF_has_imposed_value( i, ic ) )
            {
               f->modify_DOF( false, values( ic, i ), i, ic ) ;
            }
         }
      }
   }

   cv->destroy() ; cv = 0 ;
   fe->destroy() ; fe = 0 ;
}

//internal--------------------------------------------------------------
void 
PDE_DOFsInitializer_ERROR:: n0( std::string const& name, size_t nbc,
                                size_t given_nbc )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "field \"" << name << "\" :" << std::endl
        << "   declared with " << nbc << " component(s)" << std::endl 
        << "   initialized with " << given_nbc << " component(s)" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_DOFsInitializer_ERROR:: n1( std::string const& name,
                                size_t ic,
                                GE_Point const* pt,
                                double val1, double val2 )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "field \"" << name << "\", component " << ic << " : " << std::endl ;
   mesg << "   conflict at node : " ;
   pt->print( mesg, 0 ) ; mesg << std::endl ;
   mesg << "   two possible values : " ;
   mesg << val1 << " and " << val2 << std::endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_DOFsInitializer_ERROR:: n3( std::string const& name,  size_t ic )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "field \"" << name << "\", component " << ic << " : " << std::endl ;
   mesg << "   it exists nodes without initial value defined" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_DOFsInitializer_ERROR:: n4( std::string const& name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "field \"" << name << "\" : " << std::endl ;
   mesg << "   no imposed DOF value for a boundary field" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_DOFsInitializer_ERROR:: n5( PEL_ModuleExplorer const* exp,
                                size_t nb_comps )
//internal--------------------------------------------------------------
{
   std::ostringstream ok ;
   ok << "{ 0" ;
   for( size_t i=1 ; i<nb_comps ; ++i )
   {
      ok << ", " << i ;
   }
   ok << " }" ;
   PEL_Error::object()->raise_bad_data_value( exp, "component", ok.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_DOFsInitializer_ERROR:: n6( std::string const& name, size_t nbc,
                                size_t given_nbc )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "DOFs_imposed_value for field \"" << name << "\" :" << std::endl
        << "   declared for " << nbc << " component(s)" << std::endl 
        << "   initialized with " << given_nbc << " component(s)" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_DOFsInitializer_ERROR:: n7( std::string const& name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "DOFs_imposed_value for field \"" << name << "\" :" << std::endl
        << "   an entry whose keyword is \"indicator\" and" << std::endl
        << "   and whose data is a boolean expression of $DV_X" << std::endl
        << "   is requested" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_DOFsInitializer_ERROR:: n8( std::string const& name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "DOFs_imposed_value for field \"" << name << "\" :" << std::endl
        << "   the \"in_region\" location is not" << std::endl
        << "   defined for multilevel approximations" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
