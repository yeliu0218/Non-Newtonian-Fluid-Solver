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

#include <PDE_DomainAndFields.hh>

#include <GE_Color.hh>
#include <GE_SetOfPoints.hh>

#include <PDE_AdapterCHARMS.hh>
#include <PDE_AdapterHN.hh>
#include <PDE_CrossProcessBFNumbering.hh>
#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainBuilder.hh>
#include <PDE_FieldComposition.hh>
#include <PDE_DOFsInitializer.hh>
#include <PDE_GridFE.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalFEinterface.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfBasisFunctions.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfFieldCompositions.hh>

#include <PEL.hh>
#include <PEL_Bool.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>
#include <PEL_Root.hh>
#include <PEL_String.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;

//-------------------------------------------------------------------------
PDE_DomainAndFields*
PDE_DomainAndFields:: create( PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   PDE_DomainAndFields* result = new PDE_DomainAndFields( a_owner, exp, 0 ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_distributed() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_DomainAndFields*
PDE_DomainAndFields:: create( PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp,
                              PEL_Communicator const* com )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( com != 0 ) ;

   PDE_DomainAndFields* result = new PDE_DomainAndFields( a_owner, exp, com ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->is_distributed() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_DomainAndFields:: PDE_DomainAndFields( PEL_Object* a_owner,
                                           PEL_ModuleExplorer const* exp,
                                           PEL_Communicator const* com )
//-------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , EXP( exp->create_clone( this ) )
   , DISC( exp->string_data( "type" ) )
   , BUILDER_FE( 0 )
   , RSAVER( 0 )
   , DECOHESION_EXP( 0 )
   , DA( 0 )
   , TA( 0 )
   , COM( com )
{
   std::string my_name = "unnamed" ;
   if( exp->has_entry( "name" ) ) my_name = exp->string_data( "name" ) ;

   BUILDER_FE = PDE_DomainBuilder::create( this, exp, my_name ) ;
   BUILDER_FE->set_facade( this ) ;

   if( COM != 0 )
   {
      PDE_SetOfDiscreteFields const* fields =
                                     BUILDER_FE->set_of_discrete_fields() ;
      for( fields->start() ; fields->is_valid() ; fields->go_next() )
      {
         PDE_DiscreteField* ff = fields->item() ;
         PDE_CrossProcessNodeNumbering* numbering =
                             create_CrossProcessNodeNumbering( ff ) ;
         ff->build_cross_process_globalization( numbering ) ;
      }

      PDE_SetOfBasisFunctions* bfs = BUILDER_FE->set_of_basis_functions() ;
      PDE_CrossProcessBFNumbering* bf_numbering =
                             create_CrossProcessBFNumbering( bfs ) ;
      bfs->build_cross_process_globalization( bf_numbering ) ;
   }

   PDE_DOFsInitializer::object()->initialize_discrete_fields( this, exp ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
PDE_DomainAndFields:: ~PDE_DomainAndFields( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
std::string const&
PDE_DomainAndFields:: name( void ) const
//-------------------------------------------------------------------------
{
   std::string const& result = BUILDER_FE->name() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainAndFields:: duplicate_field(
                                 std::string const& model_name,
                                 std::string const& name_of_new_field ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: duplicate_field(with copy)" ) ;
   PEL_CHECK_PRE( set_of_discrete_fields()->has( model_name ) ) ;
   PEL_CHECK_PRE( !set_of_discrete_fields()->has( name_of_new_field ) ) ;

   BUILDER_FE->duplicate_field( model_name, name_of_new_field ) ;
   if( COM != 0 )
   {
      PDE_DiscreteField* f =
         BUILDER_FE->set_of_discrete_fields()->item( name_of_new_field ) ;
      PDE_CrossProcessNodeNumbering* numbering =
                          create_CrossProcessNodeNumbering( f ) ;
      f->build_cross_process_globalization( numbering ) ;
   }

   PEL_CHECK_POST( set_of_discrete_fields()->has( name_of_new_field ) ) ;
   PEL_CHECK_POST(
      set_of_discrete_fields()->item( name_of_new_field )->owner()==
         set_of_discrete_fields() ) ;
}
//-------------------------------------------------------------------------
void
PDE_DomainAndFields:: duplicate_field(
                                 std::string const& model_name,
                                 std::string const& name_of_new_field,
                                 size_t nb_components_of_new_field,
                                 size_t storage_depth_of_new_field ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: duplicate_field" ) ;
   PEL_CHECK_PRE( set_of_discrete_fields()->has( model_name ) ) ;
   PEL_CHECK_PRE( !set_of_discrete_fields()->has( name_of_new_field ) ) ;
   PEL_CHECK_PRE( nb_components_of_new_field > 0 ) ;
   PEL_CHECK_PRE( storage_depth_of_new_field > 0 ) ;

   BUILDER_FE->duplicate_field( model_name, name_of_new_field,
                                nb_components_of_new_field,
                                storage_depth_of_new_field ) ;
   if( COM != 0 )
   {
      PDE_DiscreteField* f =
         BUILDER_FE->set_of_discrete_fields()->item( name_of_new_field ) ;
      PDE_CrossProcessNodeNumbering* numbering =
                          create_CrossProcessNodeNumbering( f ) ;
      f->build_cross_process_globalization( numbering ) ;
   }

   PEL_CHECK_POST( set_of_discrete_fields()->has( name_of_new_field ) ) ;
   PEL_CHECK_POST(
      set_of_discrete_fields()->item( name_of_new_field )->owner()==
         set_of_discrete_fields() ) ;
   PEL_CHECK_POST(
      set_of_discrete_fields()->item( name_of_new_field )->nb_components() ==
         nb_components_of_new_field ) ;
   PEL_CHECK_POST(
      set_of_discrete_fields()->item( name_of_new_field )->storage_depth() ==
         storage_depth_of_new_field ) ;
}

//----------------------------------------------------------------------
void
PDE_DomainAndFields:: append_fields( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: append_fields" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   BUILDER_FE->append_fields( exp ) ;
   if( COM != 0 )
   {
      PDE_SetOfDiscreteFields const* fields =
                                     BUILDER_FE->set_of_discrete_fields() ;
      for( fields->start() ; fields->is_valid() ; fields->go_next() )
      {
         PDE_DiscreteField* f = fields->item() ;
         if( ! f->is_distributed() )
         {
            PDE_CrossProcessNodeNumbering* numbering =
                                create_CrossProcessNodeNumbering( f ) ;
            f->build_cross_process_globalization( numbering ) ;
         }
      }
   }
   PDE_DOFsInitializer::object()->initialize_discrete_fields( this, exp ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainAndFields:: apply_requests_of_DOFs_values_modules(
                                                           bool verbose ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: apply_requests_of_DOFs_values_modules" ) ;

   PDE_DOFsInitializer::object()->set_discrete_fields_DOFs_values(
                                                      this, EXP, verbose ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainAndFields:: apply_requests_of_DOFs_imposed_value_modules(
                                                           bool verbose ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: apply_requests_of_DOFs_imposed_value_modules" ) ;

   PDE_DOFsInitializer::object()->set_discrete_fields_imposed_DOFs(
                                                      this, EXP, verbose ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_DomainAndFields:: nb_conformal_adjacent_domains( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: nb_conformal_adjacent_domains" ) ;

   size_t result = BUILDER_FE->nb_conformal_adjacencies() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_DomainAndFields const*
PDE_DomainAndFields:: conformal_adjacent_domain( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: conformal_adjacent_domain" ) ;
   PEL_CHECK_PRE( i < nb_conformal_adjacent_domains() ) ;

   PDE_DomainAndFields const* result =
                         BUILDER_FE->conformal_adjacency( i )->facade() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_CrossProcessNodeNumbering*
PDE_DomainAndFields:: create_CrossProcessNodeNumbering(
                                                 PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: create_CrossProcessNodeNumbering" ) ;
   PEL_CHECK_PRE( is_distributed() ) ;

   PDE_CrossProcessNodeNumbering* result =
            new PDE_CrossProcessNodeNumbering( a_owner, BUILDER_FE, COM ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->field() == 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_CrossProcessBFNumbering*
PDE_DomainAndFields:: create_CrossProcessBFNumbering(
                                               PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: create_CrossProcessBFNumbering" ) ;
   PEL_CHECK_PRE( is_distributed() ) ;

   PDE_CrossProcessBFNumbering* result =
            new PDE_CrossProcessBFNumbering( a_owner,  BUILDER_FE, COM ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_DomainAndFields:: is_distributed( void ) const
//-------------------------------------------------------------------------
{
      return( COM != 0 ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_DomainAndFields:: nb_space_dimensions( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: nb_space_dimensions" ) ;
   size_t result = BUILDER_FE->nb_space_dimensions() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_ResultSaver*
PDE_DomainAndFields:: result_saver( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: result_saver" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( RSAVER == 0 )
   {
      PEL_ModuleExplorer const* sexp =
                           EXP->create_subexplorer( 0, "PDE_ResultSaver" ) ;
      RSAVER = new PDE_ResultSaver( const_cast<PDE_DomainAndFields*>( this ),
                                    sexp ) ;
      sexp->destroy() ;
   }
   PDE_ResultSaver* result = RSAVER ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->attached_domain() == this ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_ResultSaver*
PDE_DomainAndFields:: novel_result_saver(
                                    PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: novel_result_saver" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_ResultSaver* result = new PDE_ResultSaver(
                             const_cast<PDE_DomainAndFields*>( this ), exp ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->attached_domain() == this ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_SetOfDiscreteFields*
PDE_DomainAndFields:: set_of_discrete_fields( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: set_of_discrete_fields" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_SetOfDiscreteFields* result = BUILDER_FE->set_of_discrete_fields() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST(
      FORALL( ( result->start() ; result->is_valid() ; result->go_next() ),
              result->item()->owner() == result ) ) ;
   PEL_CHECK_POST(
      FORALL( ( result->start() ; result->is_valid() ; result->go_next() ),
              EQUIVALENT( is_distributed(),
                          result->item()->is_distributed() ) ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_SetOfFieldCompositions const*
PDE_DomainAndFields:: set_of_field_compositions( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: set_of_field_compositions" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_SetOfFieldCompositions const* result =
                                     BUILDER_FE->set_of_field_compositions() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST(
      FORALL( ( result->start() ; result->is_valid() ; result->go_next() ),
              result->item()->owner() == result ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_DomainAndFields:: is_defined_on_boundary( std::string const& fname ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: is_defined_on_boundary" ) ;

   bool result = BUILDER_FE->is_defined_on_boundary( fname ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_SetOfBCs const*
PDE_DomainAndFields:: set_of_boundary_conditions( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: set_of_boundary_conditions" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_SetOfBCs const* result = BUILDER_FE->set_of_boundary_conditions() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
GE_SetOfPoints const*
PDE_DomainAndFields:: set_of_vertices( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: set_of_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_SetOfPoints* result = BUILDER_FE->set_of_vertices() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PEL_ModuleExplorer const*
PDE_DomainAndFields:: decohesion_explorer( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: decohesion_explorer" ) ;

   if( DECOHESION_EXP == 0 )
   {
      if( BUILDER_FE->has_explorer( "decohesion" ) )
             DECOHESION_EXP = BUILDER_FE->create_explorer(
                                 const_cast<PDE_DomainAndFields*>( this ),
                                 "decohesion" ) ;
   }
   PEL_ModuleExplorer const* result = DECOHESION_EXP ;

   PEL_CHECK_POST( IMPLIES( result!=0 , result->owner()==this ) ) ;
   PEL_CHECK_POST( IMPLIES( result!=0 , result->name()=="decohesion" ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_LocalFEcell*
PDE_DomainAndFields:: create_LocalFEcell( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: create_LocalFEcell" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_LocalFEcell* result =
      new PDE_LocalFEcell( a_owner, BUILDER_FE->finite_element_grid() ) ;
   result->exclude_color( GE_Color::halo_color() ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( !result->has_foot_finder() ) ;
   PEL_CHECK_POST( result->is_excluded( GE_Color::halo_color() ) ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_CursorFEside*
PDE_DomainAndFields:: create_CursorFEside( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: create_CursorFEside" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_GridFE const* gg = BUILDER_FE->finite_element_grid() ;

   PDE_LocalFEcell* cfe_1 = create_LocalFEcell( 0 ) ;
   cfe_1->include_color( GE_Color::halo_color() ) ;

   PDE_LocalFEcell* cfe_2 = create_LocalFEcell( 0 ) ;
   cfe_2->include_color( GE_Color::halo_color() ) ;

   PDE_CursorFEside* result = new PDE_CursorFEside( a_owner,
                                                    gg,
                                                    cfe_1,
                                                    cfe_2 ) ;
   result->exclude_color( GE_Color::halo_color() ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( result->is_excluded( GE_Color::halo_color() ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_LocalFEbound*
PDE_DomainAndFields:: create_LocalFEbound( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: create_LocalFEbound" ) ;
   PEL_CHECK_INV( invariant() ) ;

   PDE_LocalFEbound* result =
      new PDE_LocalFEbound( a_owner, BUILDER_FE->finite_element_grid() ) ;
   result->exclude_color( GE_Color::halo_color() ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( result->is_excluded( GE_Color::halo_color() ) ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_LocalFEinterface*
PDE_DomainAndFields:: create_LocalFEinterface( PEL_Object* a_owner,
                                            GE_Color const* inward_domain,
                                            GE_Color const* outward_domain ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: create_LocalFEinterface" ) ;
   PEL_CHECK_PRE( !inward_domain->is_matching( outward_domain ) ) ;

   PEL_CHECK_INV( invariant() ) ;

   PDE_LocalFEinterface* result =
                 new PDE_LocalFEinterface( a_owner,
                                           BUILDER_FE->finite_element_grid(),
                                           inward_domain, outward_domain ) ;
   result->exclude_color( GE_Color::halo_color() ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   PEL_CHECK_POST( result->is_excluded( GE_Color::halo_color() ) ) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainAndFields:: save_state( PEL_ObjectWriter* writer ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: save_state" ) ;
   PEL_CHECK_PRE( save_state_PRE( writer ) ) ;

   writer->start_new_object( "PDE_DomainAndFields" ) ;

   // Discretization status :
   writer->add_entry( "nb_space_dimensions",
                      PEL_Int::create( 0, nb_space_dimensions() ) ) ;
   writer->add_entry( "type", PEL_String::create( 0, DISC ) ) ;

   // Fields storing :
   set_of_discrete_fields()->save_state( writer ) ;

   // Vertex manager :
   set_of_vertices()->save_state( writer ) ;

   writer->finalize_object() ;

   PEL_CHECK_POST( save_state_POST( writer ) ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainAndFields:: restore_state( PEL_ObjectReader* reader )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: restore_state" ) ;
   PEL_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "PDE_DomainAndFields" ) ;

   // Does some checks
   PEL_ASSERT( reader->data_of_entry( "nb_space_dimensions" )->to_int()==
                                            (int) nb_space_dimensions() ) ;
   PEL_ASSERT( reader->data_of_entry( "type" )->to_string()==DISC ) ;

   // Retrieving fields :
   set_of_discrete_fields()->restore_state( reader ) ;

   // Retrieving vertex manager :
   BUILDER_FE->set_of_vertices()->restore_state( reader ) ;

   reader->end_object_retrieval() ;

   PEL_CHECK_POST( restore_state_POST( reader ) ) ;
}

//-------------------------------------------------------------------------
bool
PDE_DomainAndFields:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( nb_space_dimensions()==1 ||
               nb_space_dimensions()==2 ||
               nb_space_dimensions()==3 ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
PDE_AdapterCHARMS*
PDE_DomainAndFields:: adapter_CHARMS( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: adapter_CHARMS" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( DA==0 )
   {
      if( EXP->has_module( "PDE_AdapterCHARMS" ) )
      {
         PEL_ModuleExplorer const* sexp =
                    EXP->create_subexplorer( 0, "PDE_AdapterCHARMS" ) ;

         DA = new PDE_AdapterCHARMS( const_cast<PDE_DomainAndFields*>( this ),
                                     this,
                                     BUILDER_FE->finite_element_grid(),
                                     BUILDER_FE->set_of_basis_functions(),
                                     sexp ) ;
         sexp->destroy() ;
      }
   }
   PDE_AdapterCHARMS* result = DA ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( result!=0, result->owner()==this ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_AdapterHN*
PDE_DomainAndFields:: adapter_HN( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: adapter_HN" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( TA==0 )
   {
      if( EXP->has_module( "PDE_AdapterHN" ) )
      {
         PEL_ModuleExplorer const* sexp =
                    EXP->create_subexplorer( 0, "PDE_AdapterHN" ) ;

         TA = new PDE_AdapterHN( const_cast<PDE_DomainAndFields*>( this ),
                                 this,
                                 BUILDER_FE->finite_element_grid(),
                                 BUILDER_FE->set_of_basis_functions(),
                                 sexp ) ;
         sexp->destroy() ;
      }
   }
   PDE_AdapterHN* result = TA ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( IMPLIES( result!=0, result->owner()==this ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainAndFields:: print_grid( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   BUILDER_FE->finite_element_grid()->print( os, indent_width ) ;
}
