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

#include <PDE_DomainBuilder.hh>

#include <GE_Color.hh>
#include <GE_Meshing.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Transform.hh>

#include <PEL.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>

#include <PDE_BasisFunctionBound.hh>
#include <PDE_BasisFunctionCell.hh>
#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_FieldComposition.hh>
#include <PDE_GridFE.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_ResultReader.hh>
#include <PDE_SetOfBasisFunctions.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfFieldCompositions.hh>

#include <doubleArray2D.hh>
#include <stringVector.hh>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::string ; using std::ostringstream ; using std::setw ;
using std::map ;

struct PDE_DomainBuilder_ERROR
{
   static void n1( std::string const& fname, std::string const& ename,
		   std::string const& mname ) ;
   static void n2( std::string const& name ) ;
   static void n5( std::string const& fname, std::string const& mname ) ;
   static void n6( PEL_ModuleExplorer const* exp ) ;
   static void n7( PEL_ModuleExplorer const* exp ) ;
   static void n8( PEL_ModuleExplorer const* exp ) ;
   static void n9( std::string const& fname ) ;
} ;

//----------------------------------------------------------------------
std::map< std::string, PDE_DomainBuilder* >&
PDE_DomainBuilder::INSTANCES( void )
//----------------------------------------------------------------------
{
   static std::map< std::string, PDE_DomainBuilder* > result ;
   return result ;
}

//----------------------------------------------------------------------
PDE_DomainBuilder*
PDE_DomainBuilder:: object( std::string const& a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: object" ) ;

   PDE_DomainBuilder* result = 0 ;

   map<string,PDE_DomainBuilder*>::const_iterator it =
                                                  INSTANCES().find( a_name ) ;
   if( it == INSTANCES().end() )
   {
      std::ostringstream ss ;
      ss << "domain \"" << a_name << "\" does not exists" ;
      PEL_Error::object()->raise_plain( ss.str() ) ;
   }
   else
   {
      result = (*it).second ;
   }
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DomainBuilder*
PDE_DomainBuilder:: create( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp,
                            std::string const& a_name )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: create" ) ;

   if( INSTANCES().count( a_name ) != 0 )
   {
      std::ostringstream ss ;
      ss << "domain \"" << a_name << "\" already exists" ;
      PEL_Error::object()->raise_plain( ss.str() ) ;
   }

   PDE_DomainBuilder* result = new PDE_DomainBuilder( a_owner, exp, a_name ) ;

   INSTANCES()[ a_name ] = result ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DomainBuilder:: PDE_DomainBuilder( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp,
                                       std::string const& a_name )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , EXP( exp->create_clone( this ) )
   , NAME( a_name )
   , DOM( 0 )
   , VERB( exp->int_data( "verbose_level" ) )
   , DIM( exp->int_data( "nb_space_dimensions" ) )
   , GRID( 0 )
   , FIELDS( PDE_SetOfDiscreteFields::create( this ) )
   , FIELD_COMPOS( 0 )
   , BCS( 0 )
   , BFS( PDE_SetOfBasisFunctions::create( this ) )
   , VMAN( 0 )
   , FIELDS_ON_BOUNDARY( 0 )
   , W_PT( 0 )
{
   PEL_LABEL( "PDE_DomainBuilder:: PDE_DomainBuilder" ) ;

   exp->test_data_in( "nb_space_dimensions", "1,2,3" ) ;
   
   W_PT = GE_Point::create( this, DIM ) ;
   VMAN = GE_SetOfPoints::create( this, DIM ) ;
   
   build_all() ;
}

//----------------------------------------------------------------------
PDE_DomainBuilder:: ~PDE_DomainBuilder( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: ~PDE_DomainBuilder" ) ;

   map<string,PDE_DomainBuilder*>::iterator it = INSTANCES().find( NAME ) ;
   PEL_CHECK( it != INSTANCES().end() ) ;
   INSTANCES().erase( it ) ;
}

//----------------------------------------------------------------------
std::string const&
PDE_DomainBuilder:: name( void ) const
//----------------------------------------------------------------------
{
   return( NAME ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainBuilder:: duplicate_field(
                                 std::string const& model_name,
                                 std::string const& name_of_new_field ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: duplicate_field(with copy)" ) ;
   PEL_CHECK_PRE( set_of_discrete_fields()->has( model_name ) ) ;
   PEL_CHECK_PRE( !set_of_discrete_fields()->has( name_of_new_field ) ) ;

   PDE_SetOfDiscreteFields* sdf = set_of_discrete_fields() ;
   
   size_t i = sdf->index_of( model_name ) ;
   
   PDE_DiscreteField* model_f = sdf->item( i ) ;
   bool on_bounds = sdf->boundary_located( i ) ;
   PEL_Vector const* model_elms = sdf->reference_elements( i ) ;
   
   PDE_DiscreteField* new_f =
                      model_f->create_duplication( 0, name_of_new_field ) ;
   sdf->append( new_f,
                on_bounds,
                model_elms->create_clone( 0 ) ) ;

   GRID->duplicate_discretization( model_f, new_f, on_bounds ) ;
   
   PEL_CHECK_POST( set_of_discrete_fields()->has( name_of_new_field ) ) ;
   PEL_CHECK_POST(
      set_of_discrete_fields()->item( name_of_new_field )->owner() ==
         set_of_discrete_fields() ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainBuilder:: duplicate_field(
                                 std::string const& model_name,
                                 std::string const& name_of_new_field,
                                 size_t nb_components_of_new_field,
                                 size_t storage_depth_of_new_field ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: duplicate_field" ) ;
   PEL_CHECK_PRE( set_of_discrete_fields()->has( model_name ) ) ;
   PEL_CHECK_PRE( !set_of_discrete_fields()->has( name_of_new_field ) ) ;
   PEL_CHECK_PRE( nb_components_of_new_field > 0 ) ;
   PEL_CHECK_PRE( storage_depth_of_new_field > 0 ) ;

   PDE_SetOfDiscreteFields* sdf = set_of_discrete_fields() ;

   size_t i = sdf->index_of( model_name ) ;
   
   PDE_DiscreteField const* model_f = sdf->item( i ) ;
   bool on_bounds = sdf->boundary_located( i ) ;
   PEL_Vector const* model_elms = sdf->reference_elements( i ) ;
   
   PDE_DiscreteField* new_f = model_f->create_duplication( 0, 
                                              name_of_new_field,
                                              nb_components_of_new_field,
                                              storage_depth_of_new_field ) ;
   
   
   new_f->set_nb_nodes( model_f->nb_nodes() ) ;
   sdf->append( new_f,
                on_bounds,
                model_elms->create_clone( 0 ) ) ;

   GRID->duplicate_discretization( model_f, new_f, on_bounds ) ;

   PEL_CHECK_POST( set_of_discrete_fields()->has( name_of_new_field ) ) ;
   PEL_CHECK_POST(
      set_of_discrete_fields()->item( name_of_new_field )->owner() ==
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
PDE_DomainBuilder:: append_fields( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: append_fields" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   build_set_of_discrete_fields( exp ) ;
   
   build_discretizations( exp ) ;
}

//----------------------------------------------------------------------
bool
PDE_DomainBuilder:: has_explorer( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   return( EXP->has_module( path_and_name ) ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PDE_DomainBuilder:: create_explorer( PEL_Object* a_owner,
                                     std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   return( EXP->create_subexplorer( a_owner, path_and_name ) ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainBuilder:: set_facade( PDE_DomainAndFields* a_dom )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: set_facade" ) ;
   PEL_CHECK_PRE( facade() == 0 ) ;

   DOM = a_dom ;

   PEL_CHECK_POST( facade() == a_dom ) ;
}

//-------------------------------------------------------------------------
PDE_DomainAndFields*
PDE_DomainBuilder:: facade( void ) const
//-------------------------------------------------------------------------
{
   return( DOM ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_DomainBuilder:: nb_conformal_adjacencies( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: nb_conformal_adjacent_domains" ) ;

   size_t result = ADJ.size() ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_DomainBuilder*
PDE_DomainBuilder:: conformal_adjacency( size_t i ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: conformal_adjacent_domain" ) ;
   PEL_CHECK_PRE( i < nb_conformal_adjacencies() ) ;

   PDE_DomainBuilder* result = ADJ[i] ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_DomainBuilder:: insert_conformal_adjacency( PDE_DomainBuilder* other )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields:: insert_conformal_adjacency" ) ;

   PEL_ASSERT( std::find( ADJ.begin(), ADJ.end(), other ) == ADJ.end() ) ;
   ADJ.push_back( other ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DomainBuilder:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   return( DIM ) ;
}

//----------------------------------------------------------------------
PDE_SetOfDiscreteFields*
PDE_DomainBuilder:: set_of_discrete_fields( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: set_of_discrete_fields" ) ;

   PDE_SetOfDiscreteFields* result = FIELDS ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   PEL_CHECK_POST(
      FORALL( ( result->start() ; result->is_valid() ; result->go_next() ),
              result->item()->owner() == result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SetOfFieldCompositions const*
PDE_DomainBuilder:: set_of_field_compositions( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: set_of_field_compositions" ) ;

   PDE_SetOfFieldCompositions const* result = FIELD_COMPOS ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   PEL_CHECK_POST(
      FORALL( ( result->start() ; result->is_valid() ; result->go_next() ),
              result->item()->owner() == result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_DomainBuilder:: is_defined_on_boundary( std::string const& fname ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: is_defined_on_boundary" ) ;

   bool result = ( FIELDS_ON_BOUNDARY.has( fname ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SetOfBCs const*
PDE_DomainBuilder:: set_of_boundary_conditions( void ) const

//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: set_of_boundary_conditions" ) ;

   PDE_SetOfBCs const* result = BCS ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_SetOfPoints*
PDE_DomainBuilder:: set_of_vertices( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: set_of_vertices" ) ;

   GE_SetOfPoints* result = VMAN ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SetOfBasisFunctions*
PDE_DomainBuilder:: set_of_basis_functions( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: set_of_basis_functions" ) ;

   PDE_SetOfBasisFunctions* result = BFS ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_GridFE*
PDE_DomainBuilder:: finite_element_grid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: finite_element_grid" ) ;
   if( GRID == 0 )
   {
      PEL_Error::object()->raise_plain( "Finite Element Grid not available" ) ;
   }
   PDE_GridFE* result = GRID ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_DomainBuilder:: build_all( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: build_all" ) ;

   if( VERB!=0 && name()!="unnamed" )
      PEL::out() << "*** Domain : \"" << name() << "\"" << endl << endl ;

   GE_Meshing* tria = create_meshing( 0 ) ;
   GRID = PDE_GridFE::create( this, VMAN, tria->nb_cells(), tria->nb_faces() ) ;
      
   build_set_of_discrete_fields( EXP ) ;
   
   build_grid( tria ) ;
   
   tria->destroy() ; tria = 0 ;

   build_discretizations( EXP ) ;

   BCS = PDE_SetOfBCs::create( this, EXP, FIELDS ) ;

   build_field_compositions() ;
}

//----------------------------------------------------------------------
GE_Meshing*
PDE_DomainBuilder:: create_meshing( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: create_meshing" ) ;
   
   if( EXP->has_entry( "check_meshing_consistency" ) )
   {
      bool check_consistency = EXP->bool_data( "check_meshing_consistency" ) ;
      if( check_consistency )
         GE_Mpolyhedron::set_check_consistency() ;
      else
         GE_Mpolyhedron::unset_check_consistency() ;
   }
   
   PEL_ModuleExplorer const* se = EXP->create_subexplorer( 0, "GE_Meshing" ) ;

   if( VERB!=0 ) PEL::out()<< "*** Building grid" << endl ;
   GE_Meshing* result = GE_Meshing::create( 0, se, DIM ) ;
   if( VERB!=0 ) result->print(  PEL::out(), 7 ) ;

   if( result->nb_cells() == 0 )
   {
      PDE_DomainBuilder_ERROR::n8( se ) ;
   }
   else if( result->nb_faces() == 0 || result->nb_vertices() == 0 )
   {
      PEL_Error::object()->raise_internal( "*** PDE_DomainBuilder error:\n"
                                           "    inconsistent meshing" ) ;
   }
   se->destroy() ; se=0 ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_DomainBuilder:: build_set_of_discrete_fields( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: build_set_of_discrete_fields" ) ;
   
   if( !exp->has_module( "interior_fields" ) &&
       !exp->has_module( "boundary_fields" ) )
   {
      PDE_DomainBuilder_ERROR::n2( exp->name() ) ;
   }
   
   // fields defined on the interior of the domain
   // --------------------------------------------
   if( exp->has_module( "interior_fields" ) )
   {
      PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "interior_fields" ) ;
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer const* se = ee->create_subexplorer( 0 ) ;
         PDE_DiscreteField* ff = 
             PDE_DiscreteField::create( 0,
                                        se->string_data( "name" ),
                                        se->int_data( "nb_components" ),
                                        se->int_data( "storage_depth" ) ) ;
         
         PEL_Vector* elms = create_elements( 0, se ) ;
         FIELDS->append( ff, false, elms ) ;
         GRID->extend_mesh_discretizations( ff, false, elms ) ;
         
         se->destroy() ; se = 0 ;
      }
      ee->destroy() ; ee=0 ;
   }
   
   // fields defined on the boundary of the domain
   // --------------------------------------------
   if( exp->has_module( "boundary_fields" ) )
   {
      PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "boundary_fields" ) ;
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer const* se = ee->create_subexplorer( 0 ) ; 
         PDE_DiscreteField* ff = 
             PDE_DiscreteField::create( 0,
                                        se->string_data( "name" ),
                                        se->int_data( "nb_components" ),
                                        se->int_data( "storage_depth" ) ) ;
         PEL_Vector* elms = create_elements( 0, se ) ;
         FIELDS->append( ff, true, elms ) ;
         GRID->extend_mesh_discretizations( ff, true, elms ) ;
         se->destroy() ; se = 0 ;
      }
      ee->destroy() ; ee=0 ;
   }
}

//----------------------------------------------------------------------
void
PDE_DomainBuilder:: build_grid( GE_Meshing* tria )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: build_grid" ) ;
   
   // vertices
   // --------
   if( VERB>1 ) PEL::out() << "       building vertices..." << endl ;
   tria->start_vertex_iterator() ;
   for( ; tria->valid_vertex() ; tria->go_next_vertex() )
   {
      W_PT->set_coordinates( tria->vertex_coordinates() ) ;
      VMAN->append( W_PT, tria->vertex_color() ) ;
   }

   // faces
   // -----
   if( VERB>1 ) PEL::out() << "       building faces..." << endl ;
   size_t nb_meshes = 0 ;
   tria->start_face_iterator() ;
   for( ; tria->valid_face() ; tria->go_next_face() )
   {
      GRID->build_face( nb_meshes++, 
                        tria->face_polyhedron_name(),
                        tria->face_vertices(), 
                        tria->face_color() ) ;
   }

   // cells
   // -----
   if( VERB>1 ) PEL::out() << "       building cells..." << endl ;
   nb_meshes = 0 ;
   tria->start_cell_iterator() ;
   for( ; tria->valid_cell() ; tria->go_next_cell() )
   {
      GRID->build_cell( nb_meshes++, 
                        tria->cell_polyhedron_name(),
                        tria->cell_vertices(), 
                        tria->cell_faces(), 
                        tria->cell_color() ) ;
   }

   // done here because the colors defined by the meshing must be built before
   // ---------
   build_special_colors() ;
   GRID->set_periodicity_requests( EXP ) ;

   // sides and bounds
   // ----------------
   if( VERB>1 ) PEL::out() << "       building bounds..." << endl ;
   GRID->separate_internal_and_external_faces() ;      

   NVERT = VMAN->nb_points() ;

   PEL_ASSERT( GRID->nb_cells() == tria->nb_cells() ) ;

   if( VERB!=0 )
   {
      PEL::out()<< "    Number of vertices : " << NVERT << endl ;
      PEL::out()<< "    Number of cells  : " << GRID->nb_cells()  << endl ;
      PEL::out()<< "    Number of sides  : " << GRID->nb_sides()  << endl;
      PEL::out()<< "    Number of bounds : " << GRID->nb_bounds() << endl ;
      PEL::out()<< endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_DomainBuilder:: build_special_colors( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: build_special_colors" ) ;

   // macro-colors
   // ------------
   if( EXP->has_module( "macro_colors" ) )
   {
      PEL_ModuleExplorer* se = EXP->create_subexplorer( 0, "macro_colors" ) ;
      se->start_entry_iterator() ;
      for( ; se->is_valid_entry() ; se->go_next_entry() )
      {
         string const& nn = se->keyword() ;
         stringVector const& name_list = se->data(this)->to_string_vector(0) ;
         GE_Color::extend( nn, name_list ) ;
      }
      se->destroy() ;
   }
}

//----------------------------------------------------------------------
void
PDE_DomainBuilder:: build_discretizations( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: build_discretizations" ) ;

   GRID->set_decohesion_requests( exp ) ;
   
   if( VERB!=0 )
   {
      PEL::out()<< "*** Building fields " << endl ;
      PEL::out()<< "    |        Field        | Number of nodes |" << endl ;
   }

   // fields 
   //    not yet fully built
   //    defined on the interior of the domain 
   // ----------------------------------------
   for( size_t i=0 ; i<FIELDS->nb_fields() ; ++i )
   {
      if( FIELDS->item( i )->nb_nodes() != 0 ) continue ;
      
      if( FIELDS->boundary_located( i ) ) continue ;
      
      GRID->define_mesh_as_cell_of_cluster() ;
      
      build_one_interior_field( FIELDS->item( i ),
                                FIELDS->reference_elements( i ) ) ;
   }

   // fields 
   //    not yet fully built
   //    defined on the boundary of the domain 
   // ----------------------------------------
   for( size_t i=0 ; i<FIELDS->nb_fields() ; ++i )
   {
      if( FIELDS->item( i )->nb_nodes() != 0 ) continue ;
      
      if( !FIELDS->boundary_located( i ) ) continue ;
      
      GRID->define_mesh_as_bound_of_cluster() ;
      
      build_one_boundary_field( FIELDS->item( i ),
                                FIELDS->reference_elements( i ) ) ;
   }
   
   if( VERB!=0 )
   {
      PEL::out()<< endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_DomainBuilder:: build_field_compositions( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: build_field_compositions" ) ;

   FIELD_COMPOS = PDE_SetOfFieldCompositions::create( this ) ;
   if( EXP->has_module( "field_compositions" ) )
   {
      PEL_ModuleExplorer* e =
                 EXP->create_subexplorer( 0, "field_compositions" ) ;
      FIELD_COMPOS->build_compositions( e, DIM, FIELDS ) ;
      e->destroy() ; e=0 ;
   }
}

//----------------------------------------------------------------------
void
PDE_DomainBuilder:: build_one_interior_field( PDE_DiscreteField* ff,
                                              PEL_Vector const* elms )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: build_one_interior_field" ) ;

   if( VERB != 0 ) PEL::out() << "    |" << setw( 21 ) << ff->name() << "|" ;

   // ***
   size_t elts_index ;
   BFS->extend_ref_elts_grp( elms, elts_index ) ;

   // determine nodes
   // ---------------
   
   size_t max_elm_nodes = 0 ;
   for( size_t i=0 ; i<elms->index_limit() ; ++i )
   {
      PDE_ReferenceElement const* elm = 
          static_cast< PDE_ReferenceElement const* >( elms->at( i ) ) ;
      if( elm->nb_nodes() > max_elm_nodes ) max_elm_nodes = elm->nb_nodes() ;
   }
   
   size_t nb_nodes = 0 ;
   size_t_vector nb_meshnodes( GRID->nb_cells() ) ;
   size_t_array2D node_of_mesh( max_elm_nodes, GRID->nb_cells() ) ;
   
   build_field_nodes( ff, elms, nb_nodes, nb_meshnodes, node_of_mesh ) ;

   // ***  
   ff->set_nb_nodes( nb_nodes ) ;
   
   // set local discretizations
   // -------------------------
   
   std::vector< PDE_BasisFunctionCell* >
                         built_bfs( nb_nodes, (PDE_BasisFunctionCell*)0 ) ;
   for( size_t im=0 ; im<node_of_mesh.index_bound(1) ; im++ )
   {
      PDE_CellFE* m = static_cast< PDE_CellFE* >( GRID->cells()->at( im ) ) ;
      if( nb_meshnodes( im ) != 0 )
      {
         size_t_vector locNode( nb_meshnodes( im ) ) ;
         for( size_t n=0 ; n<locNode.size() ; n++ )
         {
            locNode( n ) = node_of_mesh( n, im ) ;
         }
         PDE_ReferenceElement const* elm = ref_element( ff->name(),
                                                        elms,
                                                        m->polyhedron() ) ;
         size_t ee = m->index_of_element( elm ) ;
         PEL_ASSERT( ee != PEL::bad_index() ) ;
         PEL_ASSERT( elm->nb_nodes() == nb_meshnodes( im ) ) ;
         for( size_t ln=0 ; ln<nb_meshnodes( im ) ; ++ln )
         {
            size_t gn = node_of_mesh( ln, im ) ;
            PDE_BasisFunctionCell* bf = m->basis_function( ee, ln ) ;
            if( bf == 0 )
            {
               bf = built_bfs[ gn ] ;
               if( bf == 0 )
               {
                  bf = PDE_BasisFunctionCell::create( 0 ) ;
                  bf->set_active() ;
                  built_bfs[ gn ] = bf ;
                  bf->append_one_field( ff, gn ) ;
                  BFS->add( bf, elts_index ) ;
               }
               bf->extend_pieces( m, ee, ln ) ;
               m->set_basis_function( ee, ln, bf ) ;
            }
            else if( built_bfs[ gn ] == 0 )
            {
               bf->append_one_field( ff, gn ) ;
               built_bfs[ gn ] = bf ;
               PEL_ASSERT( bf->node_of_DOF( ff ) == gn ) ;
            }
         }
      }
      else if( m->color() != GE_Color::halo_color() )
      {
         string mesg = "Incomplete initialization of field: " + ff->name() ;
         PEL_Error::object()->raise_internal( mesg ) ;
      }
   }

   if( VERB != 0 ) PEL::out() << setw( 17 ) << ff->nb_nodes() << "|" << endl ;
}

//----------------------------------------------------------------------
void
PDE_DomainBuilder:: build_one_boundary_field( PDE_DiscreteField* ff,
                                              PEL_Vector const* elms )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: build_one_boundary_field" ) ;

   if( VERB != 0 ) PEL::out() << "    |" << setw( 21 ) << ff->name() << "|" ;

   //??? this info is now located in PDE_SetOfDiscreteFields
   //??? I (B. Piar) tried to remove it, but without success
   FIELDS_ON_BOUNDARY.extend( ff->name() ) ;

   // *** determine nodes
   size_t max_elm_nodes = 0 ;
   for( size_t i=0 ; i<elms->index_limit() ; ++i )
   {
      PDE_ReferenceElement const* elm = 
          static_cast< PDE_ReferenceElement const* >( elms->at( i ) ) ;
      if( elm->nb_nodes() > max_elm_nodes ) max_elm_nodes = elm->nb_nodes() ;
   }
   
   size_t nb_nodes = 0 ;
   size_t_vector nb_meshnodes( GRID->nb_bounds() ) ;
   size_t_array2D node_of_mesh( max_elm_nodes, GRID->nb_bounds() ) ;
   
   build_field_nodes( ff, elms, nb_nodes, nb_meshnodes, node_of_mesh ) ;

   // ***
   ff->set_nb_nodes( nb_nodes ) ;
   
   // *** set local discretizations
   
   std::vector< PDE_BasisFunctionBound* >
                         built_bfs( nb_nodes, (PDE_BasisFunctionBound*)0 ) ;
   for( size_t im=0 ; im<node_of_mesh.index_bound(1) ; im++ )
   {
      PDE_BoundFE* m = static_cast< PDE_BoundFE* >( GRID->bounds()->at(im) );
      if( nb_meshnodes( im ) != 0 )
      {
         size_t_vector locNode( nb_meshnodes( im ) ) ;
         for( size_t n=0 ; n<locNode.size() ; n++ )
         {
            locNode( n ) = node_of_mesh( n, im ) ;
         }

         PDE_ReferenceElement const* elm = ref_element( ff->name(),
                                                        elms,
                                                        m->polyhedron() ) ;
         size_t ee = m->index_of_element( elm ) ;
         PEL_ASSERT( elm->nb_nodes() == nb_meshnodes( im ) ) ;
         for( size_t ln=0 ; ln<nb_meshnodes( im ) ; ++ln )
         {
            size_t gn = node_of_mesh( ln, im ) ;
            PDE_BasisFunctionBound* bf = m->basis_function( ee, ln ) ;
            if( bf == 0 )
            {
               bf = built_bfs[ gn ] ;
               if( bf == 0 )
               {
                  bf = PDE_BasisFunctionBound::create( GRID ) ;
                  bf->set_active() ;
                  built_bfs[ gn ] = bf ;
                  bf->append_one_field( ff, gn ) ;
               }
               bf->extend_pieces( m, ee, ln ) ;
               m->set_basis_function( ee, ln, bf ) ;
            }
            else if( built_bfs[ gn ] == 0 )
            {
               bf->append_one_field( ff, gn ) ;
               built_bfs[ gn ] = bf ;
               PEL_ASSERT( bf->node_of_DOF( ff ) == gn ) ;
            }
         }
      }
      else if( m->color() != GE_Color::halo_color() )
      {
         string mesg = "Incomplete initialization of field: " + ff->name() ;
         PEL_Error::object()->raise_internal( mesg ) ;
      }
   }

   if( VERB != 0 ) PEL::out() << setw( 17 ) << ff->nb_nodes() << "|" << endl ;
}

//----------------------------------------------------------------------
PEL_Vector*
PDE_DomainBuilder:: create_elements( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: create_elements" ) ;
   
   bool elem  = ( exp->has_entry( "element_name" ) ) ;
   bool elems = ( exp->has_entry( "element_names" ) ) ;
   if(  elem &&  elems ) PDE_DomainBuilder_ERROR::n6( exp ) ;
   if( !elem && !elems ) PDE_DomainBuilder_ERROR::n7( exp ) ;

   //??? tests should be done: eg eqach reference polyhedron should
   //??? appear once and only once
   
   PEL_Vector* result = 0 ;
   if( elem )
   {
      result = PEL_Vector::create( a_owner, 1 ) ;
      string const& nn = exp->string_data( "element_name" ) ;
      PDE_ReferenceElement const* elm = PDE_ReferenceElement::object( nn ) ;
      result->set_at( 0, const_cast< PDE_ReferenceElement* >( elm ) ) ;
   }
   if( elems )
   {
      stringVector const& nn = exp->stringVector_data( "element_names" ) ;
      result = PEL_Vector::create( a_owner, nn.size() ) ;
      for( size_t i=0 ; i<nn.size() ; ++i )
      {
         PDE_ReferenceElement const* elm =
                                     PDE_ReferenceElement::object( nn( i ) ) ;
         result->set_at( i, const_cast< PDE_ReferenceElement* >( elm ) ) ;
      }
   }
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->index_limit() >= 1 ) ;
   PEL_CHECK_POST( result->count() == result->index_limit() ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<result->index_limit() ; ++i ),
         dynamic_cast< PDE_ReferenceElement const* >( result->at( i ) ) 
         != 0  ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_ReferenceElement const*
PDE_DomainBuilder:: ref_element( std::string const& field_name,
                                 PEL_Vector const* elms,
                                 GE_Mpolyhedron const* poly ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: ref_element" ) ;
   PEL_CHECK( poly != 0 ) ;
   
   PDE_ReferenceElement const* result = 0 ;
   if( elms->index_limit() == 1 )
   {
      result = static_cast< PDE_ReferenceElement* >( elms->at( 0 ) ) ;
      if( result->reference_polyhedron() != poly->reference_polyhedron() )
      {
         PDE_DomainBuilder_ERROR::n1( field_name, 
                                      result->name(), poly->name() ) ;
      }
   }
   else
   {
      for( size_t i=0 ; result == 0 && i<elms->index_limit() ; ++i )
      {
         PDE_ReferenceElement const* elm =
             static_cast< PDE_ReferenceElement* >( elms->at( i ) ) ;
         if( elm->reference_polyhedron() == poly->reference_polyhedron() )
         {
            result = elm ;
         }
      }
      if( result == 0 )
      {
         PDE_DomainBuilder_ERROR::n5( field_name, poly->name() ) ;
      }
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->reference_polyhedron() ==
                                        poly->reference_polyhedron() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_DomainBuilder:: build_field_nodes( PDE_DiscreteField const* ff,
                                       PEL_Vector const* elms,
                                       size_t& nb_nodes,
                                       size_t_vector& nb_meshnodes,
                                       size_t_array2D& node_of_mesh ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainBuilder:: build_field_nodes" ) ;
   PEL_CHECK( node_of_mesh.index_bound( 1 ) == 
                 GRID->nb_meshes_for_all_clusters() ) ;
   PEL_CHECK( nb_nodes == 0 ) ;
   PEL_CHECK( GRID->mesh_is_defined() ) ;
   
   // nb_nodes               : total number of nodes
   // nb_meshnodes( im )     : number of local nodes in the im-th mesh
   // node_of_mesh( in, im ) : global node number for the in-th local node 
   //                          of the im-th mesh 
   
   size_t node_start = 0 ;
   for( size_t icl=0 ; icl<GRID->nb_clusters() ; ++icl )
   {
      // STEP 1
      // ------
      //    construction of the set of all node points in the
      //    current cluster (the node points are recorded by PDE_GridFE)
      //
      //    idx_pt_of_mesh( in, im )
      //    ************************
      //       for the im-th mesh (which belongs to the current cluster)
      //       index (in PDE_GridFE) of the in-th node point

      size_t_array2D idx_pt_of_mesh( node_of_mesh ) ;
      
      GRID->start_nodepts_recording() ;

      GRID->start_mesh_iterator( icl ) ;
      for( ; GRID->valid_mesh() ; GRID->go_next_mesh() )
      {
         PDE_MeshFE const* m = GRID->current_mesh() ;
         size_t im  = m->id_number() ;
      
         GE_Color const* mcol = m->color() ;
         GE_Mpolyhedron const* poly = m->polyhedron() ;
         if( poly->dimension() == poly->nb_space_dimensions() ||
             mcol != GE_Color::halo_color() )
         {
            PDE_ReferenceElement const* elm = ref_element( ff->name(),
                                                           elms,
                                                           poly ) ;
            size_t nb_loc = elm->nb_nodes() ;
            nb_meshnodes( im ) = nb_loc ;

            for( size_t in = 0 ; in<nb_loc ; in++ )
            {               
               size_t idx = PEL::bad_index() ;
               bool is_new = false ;
               GRID->record_nodept( poly, elm->node_location( in ), 
                                    idx, is_new ) ;
               PEL_ASSERT( is_new ) ;
               idx_pt_of_mesh( in, im ) = idx ;
            }
         }
      }
      
      GRID->terminate_nodepts_recording() ;

      // STEP 2
      // ------
      //    assign an index to each node point (two different node points
      //    may represent the same basis function, eg due to periodicity
      //    requests)
      //
      //    node_start
      //    **********
      //       cross-cluster index of the first node nodes of the current
      //       cluster
            
      for( size_t im=0 ; im<node_of_mesh.index_bound( 1 ) ; ++im )
      {
         if( ! GRID->mesh_belongs_to_cluster( im, icl ) ) continue ; // <---
         
         size_t_vector built( 0 ) ;
         for( size_t in=0 ; in<nb_meshnodes( im ) ; ++in )
         {
            size_t nn = GRID->inferred_node_of_recorded_nodept( 
                                               idx_pt_of_mesh( in, im ) ) ;
            PEL_ASSERT( nn != PEL::bad_index() ) ;
            nn += node_start ;
            if( built.has( nn ) )
            {
               PDE_DomainBuilder_ERROR::n9( ff->name() ) ;
            }
            else
            {
               built.append( nn ) ;
               node_of_mesh( in, im ) = nn ;
            }
         }
      }
      
      node_start += GRID->nb_inferred_nodes() ;
   }
   
   nb_nodes = node_start ;
}

//internal--------------------------------------------------------------
void
PDE_DomainBuilder_ERROR:: n1( std::string const& fname,
                              std::string const& ename,
                              std::string const& mname )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "field \"" << fname << "\" :" << endl
        << "   a mesh \"" << mname << "\" and " << endl
	<< "   a reference element \"" << ename << "\"" << endl
        << "   are incompatible" << endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_DomainBuilder_ERROR:: n2( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "module \"" << name << "\" :" << endl
        << "   there should be at least a submodule called" << endl
        << "   \"interior_fields\" or" << endl
        << "   \"boundary_fields\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_DomainBuilder_ERROR:: n5( std::string const& fname,
                              std::string const& mname )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "field \"" << fname << "\" :" << endl
        << "   no reference element found for" << endl
        << "   the following cell: \"" << mname << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_DomainBuilder_ERROR:: n6( PEL_ModuleExplorer const* exp )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "the entries of keyword: " << endl ;
   mesg << "   \"element_name\" and \"element_names\"" << endl ;
   mesg << "are mutually exclusive" << endl ;
   mesg << "(one and only one of them should be present)." ;
   PEL_Error::object()->raise_module_error( exp, mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_DomainBuilder_ERROR:: n7( PEL_ModuleExplorer const* exp )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "an entry of keyword:" << endl ;
   mesg << "   \"element_name\" or \"element_names\"" << endl ;
   mesg << "is expected." ;
   PEL_Error::object()->raise_module_error( exp, mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_DomainBuilder_ERROR:: n8( PEL_ModuleExplorer const* exp )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "empty meshing encountered." ;
   PEL_Error::object()->raise_module_error( exp, mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_DomainBuilder_ERROR:: n9( std::string const& fname )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "field \"" << fname << "\" :" << endl
        << "   there exists a cell where two nodes are found" << endl
        << "   identical (probably due to periodic boundary conditions)" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
