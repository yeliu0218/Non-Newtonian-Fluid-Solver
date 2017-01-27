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

#include <PDE_GridFE.hh>

#include <PEL.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_IndexSet.hh>
#include <PEL_Iterator.hh>
#include <PEL_ListIdentity.hh>
#include <PEL_Vector.hh>
#include <boolVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_ReferencePolyhedronRefiner.hh>
#include <GE_SetOfPoints.hh>
#include <GE_SegmentSegment_INT.hh>
#include <GE_SimplePolygon2D.hh>
#include <GE_Transform.hh>

#include <PDE_BasisFunctionCell.hh>
#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DiscOnMeshFE.hh>
#include <PDE_FaceFE.hh>
#include <PDE_FacesOfCellFE.hh>
#include <PDE_ReferenceElement.hh>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::cout ; using std::endl ;
using std::setw ;
using std::ostringstream ;

struct PDE_GridFE_ERROR
{
  static void n0( GE_Color const* color ) ;
  static void n1( PDE_FaceFE const* face, PEL_Iterator* it_tr ) ;
  static void n2( PDE_FaceFE const* face1, PDE_FaceFE const* face2,
                  GE_Transform const* tr ) ;
  static void n3( PDE_FaceFE const* face1, GE_Transform const* tr ) ;
  static void n4( PDE_FaceFE const* face ) ;
  static void n5( PDE_FaceFE const* face, PEL_Iterator* it_tr,
                  size_t rank, size_t nb_ranks ) ;
  static void n6( void ) ;
  static void n7( std::string const& color_name ) ;
  static void n8( std::string const& color_name,
                  std::string const& col1, std::string const& col2 ) ;
} ;

//--------------------------------------------------------------------
PDE_GridFE*
PDE_GridFE:: create( PEL_Object* a_owner,
                     GE_SetOfPoints* a_vertex_manager,
                     size_t a_nb_cells,
                     size_t a_nb_faces )
//--------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: create" ) ;
   PEL_CHECK_PRE( a_vertex_manager!=0 ) ;
   PEL_CHECK_PRE( a_vertex_manager->dimension()==1 ||
                  a_vertex_manager->dimension()==2 ||
                  a_vertex_manager->dimension()==3 ) ;

   PDE_GridFE* result =
      new PDE_GridFE( a_owner, a_vertex_manager, a_nb_cells, a_nb_faces ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->set_of_vertices()==a_vertex_manager ) ;
   PEL_CHECK_POST( result->nb_space_dimensions()==a_vertex_manager->dimension() ) ;
   PEL_CHECK_POST( result->set_of_vertices()->has_as_observer( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_GridFE:: PDE_GridFE( PEL_Object* a_owner,
                         GE_SetOfPoints* a_vertex_manager,
                         size_t a_nb_cells,
                         size_t a_nb_faces )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , VERTS( a_vertex_manager )
   , C_DISCS( PEL_Vector::create( this, GE_ReferencePolyhedron::nb_objects() ) )
   , B_DISCS( PEL_Vector::create( this, GE_ReferencePolyhedron::nb_objects() ) )
   , CELLS( PEL_Vector::create( this, a_nb_cells ) )
   , FACES( PEL_Vector::create( 0, a_nb_faces ) )
   , SIDES( PEL_Vector::create( this, a_nb_faces ) )
   , BOUNDS( PEL_Vector::create( this, 0 ) )
   , MESH_IT( 0 )
   , OK_MESH_IT( 0 )
   , I_CURRENT_CLUSTER( PEL::bad_index() )
   , NB_MESHES( PEL::bad_index() )
   , NEW_NODEPTS( 0 )
   , PRE_NODEPTS( 0 )
   , ON_PRE( false )
   , NB_PRE_NODEPTS( 0 )
   , NB_NEW_NODEPTS( 0 ) 
   , PT_2_NODE( 0 )
   , NB_NODES( 0 ) 
   , GEO_STATE( 0 )
   , TRANSFOS( 0 )
   , TRANSFOS_IT( 0 )
   , W_PT( GE_Point::create( this, VERTS->dimension() ) )
   , FULL_DECOHESION( false )
   , NB_CLUSTERS( 1 )
   , CLUSTER_OF_MESH( 0 )
   , IDX_TEST( PEL_IndexSet::create( this ) ) 
   , VERB( 0 )
{
   PEL_LABEL( "PDE_GridFE::PDE_GridFE" ) ;
   PEL_CHECK_INV( invariant() ) ;

   VERTS->attach_observer( this ) ;
}

//----------------------------------------------------------------------
PDE_GridFE:: ~PDE_GridFE( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE::~PDE_GridFE" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( FACES != 0 )
   {
      FACES->destroy() ;
      FACES = 0 ;
   }

   VERTS->detach_observer( this ) ;   
}

//-----------------------------------------------------------------------------
void
PDE_GridFE:: update( void )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: update" ) ;
   PEL_SAVEOLD( size_t, geometric_state, geometric_state() ) ;

   // Update cells:
   for( size_t i=0 ; i<CELLS->index_limit() ; ++i )
   {
      PDE_MeshFE const* m = static_cast<PDE_MeshFE const*>( CELLS->at(i) ) ;
      if( m!=0 ) m->polyhedron()->update() ;
   }
   
   // Update sides:
   for( size_t i=0 ; i<SIDES->index_limit() ; ++i )
   {
      PDE_MeshFE const* m = static_cast<PDE_MeshFE const*>( SIDES->at(i) ) ;
      if( m!=0 ) m->polyhedron()->update() ;
   }
   
   // Update bounds:
   for( size_t i=0 ; i<BOUNDS->index_limit() ; ++i )
   {
      PDE_MeshFE const* m = static_cast<PDE_MeshFE const*>( BOUNDS->at(i) ) ;
      if( m!=0 ) m->polyhedron()->update() ;
   }

   change_geometric_state() ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( geometric_state()!=OLD(geometric_state) ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: extend_mesh_discretizations( PDE_DiscreteField const* ff,
                                          bool on_bounds,
                                          PEL_Vector const* elms )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: extend_cell_discretizations" ) ;
   
   bool elms_must_exist = false ;
   if( on_bounds )
   {
      if( BOUNDS->count() != 0 ) elms_must_exist = true ;
      extend_mesh_discretizations( elms_must_exist, B_DISCS, ff, elms ) ;
   }
   else
   {
      if( CELLS->count() != 0 ) elms_must_exist = true ;
      extend_mesh_discretizations( elms_must_exist, C_DISCS, ff, elms ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_GridFE:: duplicate_discretization( PDE_DiscreteField const* model_f,
                                       PDE_DiscreteField* new_f,
                                       bool on_bounds )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: duplicate_discretization" ) ;
   //??? preconditions
   
   if( on_bounds )
   {
      duplicate_discretization( model_f, new_f, B_DISCS, BOUNDS ) ;
   }
   else
   {
      duplicate_discretization( model_f, new_f, C_DISCS, CELLS ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_GridFE:: duplicate_discretization( PDE_DiscreteField const* model_f,
                                       PDE_DiscreteField* new_f,
                                       PEL_Vector const* discs,
                                       PEL_Vector const* meshes )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: duplicate_discretization" ) ;
   
   PEL_VectorIterator* it = discs->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_DiscOnMeshFE* disc = 
              static_cast< PDE_DiscOnMeshFE* >( it->item() ) ;
      if( disc->has_discretization( model_f ) )
      {
         disc->duplicate_discretization( model_f, new_f ) ;
      }
   }
   it->destroy() ; it = 0 ;
   
   it = meshes->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_MeshFE* mesh = static_cast< PDE_MeshFE* >( it->item() ) ;
      mesh->duplicate_discretization( model_f, new_f ) ;
   }
   it->destroy() ; it = 0 ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: extend_mesh_discretizations( bool elm_must_exist,
                                          PEL_Vector* discs,
                                          PDE_DiscreteField const* ff,
                                          PEL_Vector const* elms )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: extend_mesh_discretizations" ) ;
   
   for( size_t e=0 ; e<elms->index_limit() ; ++e )
   {
      PDE_ReferenceElement const* elm = 
          static_cast< PDE_ReferenceElement const* >( elms->at( e ) ) ;

      GE_ReferencePolyhedron const* rpoly = elm->reference_polyhedron() ;
      size_t i = rpoly->id_number() ;
   
      PEL_Object* oo = discs->at( i ) ;
      if( oo == 0 )
      {
         PEL_ASSERT( !elm_must_exist ) ;
         discs->set_at( i, PDE_DiscOnMeshFE::create( discs, ff, elm ) ) ;
      }
      else
      {
         PDE_DiscOnMeshFE* disc = static_cast<PDE_DiscOnMeshFE*>( oo ) ;
         if( elm_must_exist )
         {
            PEL_ASSERT( disc->index_of_element( elm ) != PEL::bad_index() ) ;
         }
         disc->add_discretization( ff, elm ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_GridFE:: build_face( size_t number,
                         std::string const& polyhedron_name,
                         size_t_vector const& face_vertices,
                         GE_Color const* color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: build_face" ) ;
   
   GE_Mpolyhedron* poly = GE_Mpolyhedron:: create( polyhedron_name,
                                                   VERTS,
                                                   face_vertices ) ;
   size_t rlevel = 0 ;
   PDE_FaceFE* s = PDE_FaceFE::create( VERTS, number, poly, color, rlevel ) ;
   FACES->set_at( number, s ) ; 
}

//----------------------------------------------------------------------
void
PDE_GridFE:: build_cell( size_t number,
                         std::string const& polyhedron_name,
                         size_t_vector const& cell_vertices,
                         size_t_vector const& cell_faces,
                         GE_Color const* color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: build_cell" ) ;

   GE_Mpolyhedron* poly = GE_Mpolyhedron:: create( polyhedron_name,
                                                   VERTS,
                                                   cell_vertices ) ;
   PEL_Object* oo = C_DISCS->at( poly->reference_polyhedron()->id_number() ) ;
   PDE_DiscOnMeshFE* disc = static_cast< PDE_DiscOnMeshFE* >( oo ) ;
   size_t rlevel = 0 ;
   PDE_CellFE* cell =  PDE_CellFE::create( VERTS, number, poly, color, 
                                           rlevel, disc ) ;
   
   //??? may be a precondition on cell_faces
   PEL_ASSERT( poly->nb_faces() == cell_faces.size() ) ;
   for( size_t i=0 ; i<cell_faces.size() ; ++i )
   {
      size_t iface = cell_faces( i ) ;
      PDE_FaceFE* face = static_cast< PDE_FaceFE* >( FACES->at( iface ) ) ;
      cell->append_face( face ) ;
      face->insert_adjacent_cell( cell ) ;
   }
   
   //??? must be equivalent to add_cell( cell )
   CELLS->set_at( number, cell ) ; 
}

//----------------------------------------------------------------------
void
PDE_GridFE:: set_periodicity_requests( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: set_periodicity_requests" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   
   if( exp->has_module( "domain_periodicity" ) )
   {
      TRANSFOS = PEL_List::create( this ) ;
      PEL_ModuleExplorer* ee =
                          exp->create_subexplorer( 0, "domain_periodicity" ) ;
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer const* eee = ee->create_subexplorer( 0 ) ;
         GE_Transform* a_tr = GE_Transform::create( this, 
                                                    nb_space_dimensions(), 
                                                    eee ) ;
         TRANSFOS->append( a_tr ) ;
         eee->destroy() ; eee = 0 ;
      }
      ee->destroy() ; ee = 0 ;

      if( TRANSFOS->index_limit()==0 ) PDE_GridFE_ERROR::n6() ;

      TRANSFOS_IT = TRANSFOS->create_iterator( this ) ;

      if( VERB!=0 )
      {
         PEL::out()<< "*** Periodicity" << endl ;
         TRANSFOS_IT->start() ;
         for( ; TRANSFOS_IT->is_valid() ; TRANSFOS_IT->go_next() )
         {
            TRANSFOS_IT->item()->print( PEL::out(), 6 ) ;
         }
         PEL::out() << endl ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_GridFE:: separate_internal_and_external_faces( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: separate_internal_and _external_faces" ) ;
   
   if( TRANSFOS_IT != 0 )
   {
      separate_sides_and_bounds( TRANSFOS_IT ) ;
   }
   else
   {
      separate_sides_and_bounds() ;      
   }
}

//----------------------------------------------------------------------
void
PDE_GridFE:: set_decohesion_requests( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: set_decohesion_requests" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   
   if( exp->has_module( "decohesion" ) )
   {
      PEL_ModuleExplorer const* se =
                         exp->create_subexplorer( 0, "decohesion" ) ;
      std::string const& type = se->string_data( "type" ) ;
      if( type == "everywhere" )
      {
         FULL_DECOHESION = true ;
      }
      else if( type == "inter_cell_aggregates" )
      {
         stringVector const& cn = se->stringVector_data( "aggregate_colors" ) ;
         NB_CLUSTERS = cn.size() ;
         for( size_t i=0 ; i<NB_CLUSTERS ; ++i )
         {
            CLUSTER_COLOR.push_back( GE_Color::object( cn( i ) ) ) ;
         }
      }
      else
      {
         PEL_Error::object()->raise_bad_data_value( se, "type",
                         "   \"everywhere\"\n   \"inter_cell_aggregates\"" ) ;
      }
      se->destroy() ;
   }
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: nb_clusters( void ) const
//----------------------------------------------------------------------
{
   return( NB_CLUSTERS ) ;   
}

//----------------------------------------------------------------------
void
PDE_GridFE:: define_mesh_as_cell_of_cluster( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: define_mesh_as_cell_of_cluster" ) ;
   
   build_cluster_to_mesh_connectivity( CELLS ) ;
   
   PEL_CHECK_POST( mesh_is_defined() ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: define_mesh_as_bound_of_cluster( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: define_mesh_as_bound_of_cluster" ) ;
   
   build_cluster_to_mesh_connectivity( BOUNDS ) ;
   
   PEL_CHECK_POST( mesh_is_defined() ) ;
}

//----------------------------------------------------------------------
bool
PDE_GridFE:: mesh_is_defined( void ) const
//----------------------------------------------------------------------
{
   return( MESH_IT != 0 ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: nb_meshes_for_all_clusters( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: nb_meshes_for_all_clusters" ) ;
   PEL_CHECK_PRE( mesh_is_defined() ) ;
   
   return( NB_MESHES ) ;   
}

//----------------------------------------------------------------------
bool
PDE_GridFE:: mesh_belongs_to_cluster( size_t id_number_of_mesh, 
                                      size_t i_cluster ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: cluster_of_mesh" ) ;
   PEL_CHECK_PRE( mesh_is_defined() ) ;
   PEL_CHECK_PRE( id_number_of_mesh < nb_meshes_for_all_clusters() ) ;
   
   bool result = ( CLUSTER_OF_MESH( id_number_of_mesh ) == i_cluster ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: start_mesh_iterator( size_t i_cluster )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: start_mesh_iterator" ) ;
   PEL_CHECK_PRE( mesh_is_defined() ) ;
   PEL_CHECK_PRE( i_cluster < nb_clusters() ) ;
   
   I_CURRENT_CLUSTER = i_cluster ;
   
   MESH_IT->start() ;
   while( MESH_IT->is_valid() && rejected_mesh() )
   {
      MESH_IT->go_next() ;
   }
   OK_MESH_IT = MESH_IT->is_valid() ;
}

//----------------------------------------------------------------------
bool
PDE_GridFE:: valid_mesh( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: valid_mesh" ) ;
   
   return( OK_MESH_IT ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: go_next_mesh( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: go_next_mesh" ) ;
   PEL_CHECK_PRE( valid_mesh() ) ;
   
   MESH_IT->go_next() ;
   while( MESH_IT->is_valid() && rejected_mesh() )
   {
      MESH_IT->go_next() ;
   }
   OK_MESH_IT = MESH_IT->is_valid() ;
}

//----------------------------------------------------------------------
PDE_MeshFE*
PDE_GridFE:: current_mesh( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: current_mesh" ) ;
   PEL_CHECK_PRE( valid_mesh() ) ;
   
   return( the_mesh( MESH_IT ) ) ;
}

//----------------------------------------------------------------------
void 
PDE_GridFE:: start_nodepts_recording( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: start_nodepts_recording" ) ;
   PEL_CHECK_PRE( !ongoing_nodepts_recording() ) ;

   NEW_NODEPTS = GE_SetOfPoints::create( 0, VERTS->dimension() ) ;
   NB_NEW_NODEPTS = 0 ;
   
   PEL_CHECK_POST( ongoing_nodepts_recording() ) ;
}

//----------------------------------------------------------------------
bool
PDE_GridFE:: ongoing_nodepts_recording( void ) const
//----------------------------------------------------------------------
{
   bool result = ( NEW_NODEPTS != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void 
PDE_GridFE:: start_preexisting_nodepts_recording( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: start_nodepts_recording" ) ;
   PEL_CHECK_PRE( ongoing_nodepts_recording() ) ;
   PEL_CHECK_PRE( !ongoing_preexisting_nodepts_recording() ) ;

   ON_PRE = true ;
   PRE_NODEPTS = GE_SetOfPoints::create( 0, VERTS->dimension() ) ;
   NB_PRE_NODEPTS = 0 ;
   
   PEL_CHECK_POST( ongoing_nodepts_recording() ) ;
   PEL_CHECK_POST( ongoing_preexisting_nodepts_recording() ) ;
}

//----------------------------------------------------------------------
bool
PDE_GridFE:: ongoing_preexisting_nodepts_recording( void ) const
//----------------------------------------------------------------------
{
   bool result = ( ON_PRE ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: record_preexisting_nodept( PDE_BasisFunctionCell const* bf,
                                        size_t& bf_index )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: record_excluded_nodept" ) ;
   PEL_CHECK_PRE( ongoing_nodepts_recording() ) ;
   PEL_CHECK_PRE( ongoing_preexisting_nodepts_recording() ) ;
   PEL_CHECK_PRE( bf != 0 ) ;
   
   bf->geometrical_node( W_PT ) ;
   PRE_NODEPTS->extend( W_PT ) ;
   bf_index = PRE_NODEPTS->index( W_PT ) ;
   
   if( VERB > 1 )
   {
      PEL::out() << "   old nodept: " << bf_index ;
      W_PT->print( PEL::out(), 3 ) ; PEL::out() << endl ;
   }
   
   PEL_CHECK_POST( bf_index != PEL::bad_index() ) ;
}

//----------------------------------------------------------------------
void 
PDE_GridFE:: terminate_preexisting_nodepts_recording( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: terminate_preexisting_nodepts_recording" ) ;
   PEL_CHECK_PRE( ongoing_nodepts_recording() ) ;
   PEL_CHECK_PRE( ongoing_preexisting_nodepts_recording() ) ;
      
   NB_PRE_NODEPTS = PRE_NODEPTS->nb_points() ;
   ON_PRE = false ;
   
   PEL_CHECK_POST( ongoing_nodepts_recording() ) ;
   PEL_CHECK_POST( !ongoing_preexisting_nodepts_recording() ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: record_nodept( GE_Mpolyhedron const* poly,
                            GE_Point const* pt_ref,
                            size_t& nodept_index, bool& is_new )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: record_nodept" ) ;
   PEL_CHECK_PRE( ongoing_nodepts_recording() ) ;
   PEL_CHECK_PRE( poly != 0 ) ;
   PEL_CHECK_PRE( pt_ref != 0 ) ;
   PEL_CHECK_PRE( pt_ref->nb_coordinates() == poly->dimension() ) ;
   
   nodept_index = PEL::bad_index() ;
   is_new = true ;
   if( !FULL_DECOHESION )
   {
      poly->apply_mapping( pt_ref, W_PT ) ;
      if( ( PRE_NODEPTS != 0 ) && ( PRE_NODEPTS->has( W_PT ) ) )
      {
         is_new = false ;
         nodept_index = PRE_NODEPTS->index( W_PT ) ;
      }
      else if( NEW_NODEPTS->has( W_PT ) )
      {
         nodept_index = NEW_NODEPTS->index( W_PT ) + NB_PRE_NODEPTS ;
      }
   }
   if( is_new && ( nodept_index == PEL::bad_index() ) )
   {
      nodept_index = NB_NEW_NODEPTS + NB_PRE_NODEPTS ;
      ++NB_NEW_NODEPTS ;
      if( !FULL_DECOHESION )
      {
         NEW_NODEPTS->append( W_PT ) ;
      }
   }
   
   PEL_CHECK_POST( nodept_index != PEL::bad_index() ) ;
}

//----------------------------------------------------------------------
void 
PDE_GridFE:: terminate_nodepts_recording( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: terminate_nodepts_recording" ) ;
   PEL_CHECK_PRE( ongoing_nodepts_recording() ) ;
   PEL_CHECK_PRE( !ongoing_preexisting_nodepts_recording() ) ;
      
   NB_NODES = NB_NEW_NODEPTS ;
   if( NB_NEW_NODEPTS != 0 )
   {
      PT_2_NODE.re_initialize( 0 ) ;
      if( !FULL_DECOHESION && TRANSFOS_IT != 0 )
      {
         size_t nb_pts = NEW_NODEPTS->nb_points() ;
         
         PT_2_NODE.re_initialize( nb_pts, PEL::bad_index() ) ;
         
         size_t_vector ass( nb_pts ) ;
         for( size_t ip=0 ; ip<nb_pts ; ++ip )
         {
            ass( ip ) = associated_point( NEW_NODEPTS, 
                                          TRANSFOS_IT, ip, W_PT ) ;
         }
         size_t current_node = 0 ;
         for( size_t ip=0 ; ip<nb_pts ; ++ip )
         {
            if( PT_2_NODE( ip ) == PEL::bad_index() )
            {
               fill_pt_2_node( ip, ass,
                               PT_2_NODE, current_node ) ;
            }
         }
         NB_NODES = current_node ;
      }
   }
   NEW_NODEPTS->destroy() ; 
   NEW_NODEPTS = 0 ;
   
   if( PRE_NODEPTS != 0 ) 
   {
      PEL_ASSERT( NB_PRE_NODEPTS == PRE_NODEPTS->nb_points() ) ;
      PRE_NODEPTS->destroy() ;
      PRE_NODEPTS = 0 ; 
   }
   
   PEL_CHECK_POST( !ongoing_nodepts_recording() ) ;
   PEL_CHECK_POST( !ongoing_preexisting_nodepts_recording() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: nb_recorded_nodepts( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: nb_recorded_nodepts" ) ;
   PEL_CHECK_PRE( !ongoing_nodepts_recording() ) ;
   PEL_CHECK_PRE( !ongoing_preexisting_nodepts_recording() ) ;
   
   return( NB_NEW_NODEPTS + NB_PRE_NODEPTS ) ;  
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: nb_inferred_nodes( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: nb_inferred_nodes" ) ;
   PEL_CHECK_PRE( !ongoing_nodepts_recording() ) ;
   PEL_CHECK_PRE( !ongoing_preexisting_nodepts_recording() ) ;
   
   return( NB_NODES ) ;  
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: inferred_node_of_recorded_nodept( size_t nodept_index ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: nodept_2_node" ) ;
   PEL_CHECK_PRE( !ongoing_nodepts_recording() ) ;
   PEL_CHECK_PRE( !ongoing_preexisting_nodepts_recording() ) ;
   PEL_CHECK_PRE( nodept_index < nb_recorded_nodepts() ) ;
   
   size_t result = PEL::bad_index() ;
   if( nodept_index >= NB_PRE_NODEPTS )
   {
      result = nodept_index - NB_PRE_NODEPTS ;
      if( !FULL_DECOHESION && TRANSFOS_IT != 0 )
      {
         result = PT_2_NODE( result ) ;
      }
   }
   
   PEL_CHECK_POST( ( result == PEL::bad_index() ) || 
                   ( result < nb_inferred_nodes() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: bf_index_of_recorded_nodept( size_t nodept_index ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: bf_index_of_recorded_nodept" ) ;
   PEL_CHECK_PRE( !ongoing_nodepts_recording() ) ;
   PEL_CHECK_PRE( !ongoing_preexisting_nodepts_recording() ) ;
   PEL_CHECK_PRE( nodept_index < nb_recorded_nodepts() ) ;
   
   size_t result = PEL::bad_index() ;
   if( nodept_index >= NB_PRE_NODEPTS )
   {
      result = nodept_index - NB_PRE_NODEPTS ;
      if( !FULL_DECOHESION && TRANSFOS_IT != 0 )
      {
         result = PT_2_NODE( result ) ;
      }
      result += NB_PRE_NODEPTS ;         
   }
   else
   {
      result = nodept_index ;
   }

   PEL_CHECK_POST( result < nb_recorded_nodepts() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   return( VERTS->dimension() ) ;
}

//----------------------------------------------------------------------
GE_SetOfPoints*
PDE_GridFE:: set_of_vertices( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: set_of_vertices" ) ;

   GE_SetOfPoints* result = VERTS ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->dimension()==nb_space_dimensions() ) ;
   PEL_CHECK_POST( result->has_as_observer( this ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------
size_t
PDE_GridFE:: nb_cells( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: nb_cells" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( CELLS->count() ) ;
}

//----------------------------------------------------------------------
PEL_Vector const*
PDE_GridFE:: cells( void ) const
//----------------------------------------------------------------------
{
   PEL_Vector const* result = CELLS ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: refine_cell( PDE_CellFE* ccell,
                          GE_ReferencePolyhedronRefiner const* crf )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: refine_cell" ) ;
   PEL_CHECK_PRE( ccell != 0 ) ;
   PEL_CHECK_PRE( cells()->at( ccell->id_number() ) == ccell ) ;
   PEL_CHECK_PRE( ccell->is_active() ) ;
   PEL_CHECK_PRE( ccell->nb_childs() == 0 ) ; // never been refined
   
   if( VERB > 1 ) PEL::out() << "splitting of cell "
                             << ccell->id_number() << endl ;

   // GE_ReferencePolyhedronRefiner defines vertices, subfaces, and
   //    a local numbering for them
   //
   // `idx_in_VERTS(i)' : index in `VERTS' of the vertex 
   //                     numbered`i' in GE_ReferencePolyhedronRefiner
   //
   // Among all the subfaces of ccell, some may already exist (from
   //    neighboring cells) 
   // `subfaces' : the subfaces of ccell in order of construction
   //    the first `old_subfaces_limit' items are those already existing
   //    the other are created by the refinement of ccell
   // `idx_in_subfaces(i)' : index in `subfaces' of the subface 
   //    numbered `i' in GE_ReferencePolyhedronRefiner
   //
   // `tree' : a tree of PEL_IndexSet used to compare subfaces
   //
   PEL_Vector* subfaces = PEL_Vector::create( 0, 0 ) ;
   size_t_vector idx_in_VERTS( crf->nb_vertices() ) ;
   size_t_vector idx_in_subfaces( crf->nb_subfaces() ) ;

   // Step 1: building of the new vertices and of idx_in_VERTS
   build_new_vertices( crf, ccell, idx_in_VERTS ) ;

   // Step 2: recording of the already existing subfaces and then
   //         building of the new subfaces
   // `tree' : work data structure (tree of PEL_IndexSet) used to 
   //          compare subfaces
   PEL_BalancedBinaryTree* tree = PEL_BalancedBinaryTree::create( 0 ) ;
   
   record_already_existing_faces( crf, ccell, subfaces, tree ) ;
   size_t old_subfaces_limit = subfaces->index_limit() ;
   
   build_new_faces( crf, ccell, idx_in_VERTS, tree, subfaces, idx_in_subfaces ) ;
   PEL_ASSERT( subfaces->index_limit() == crf->nb_subfaces() ) ;
   tree->destroy() ; tree = 0 ;

   // Step 3: building of the new subcells
   build_new_cells( crf, ccell, idx_in_VERTS, 
                    subfaces, old_subfaces_limit, idx_in_subfaces ) ;
   
   subfaces->destroy() ; subfaces = 0 ;

   if( VERB > 1 ) PEL::out() << "deactivating cell "<< ccell->id_number()
                                                    << endl << endl ;
   
   ccell->set_inactive() ;
}

//--------------------------------------------------------------------
size_t
PDE_GridFE:: nb_sides( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: nb_sides" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( SIDES->count() ) ;
}

//----------------------------------------------------------------------
PEL_Vector const*
PDE_GridFE:: sides( void ) const
//----------------------------------------------------------------------
{
   PEL_Vector const* result = SIDES ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: nb_bounds( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: nb_bounds" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( BOUNDS->count() ) ;
}

//----------------------------------------------------------------------
PEL_Vector const*
PDE_GridFE:: bounds( void ) const
//----------------------------------------------------------------------
{
   PEL_Vector const* result = BOUNDS ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: set_verbose_level( size_t a_verbose_level )
//----------------------------------------------------------------------
{
   VERB = a_verbose_level ;
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: geometric_state( void ) const
//----------------------------------------------------------------------
{
   return( GEO_STATE ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: change_geometric_state( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: change_geometric_state" ) ;
   PEL_SAVEOLD( size_t, geometric_state, geometric_state() ) ;
   
   ++GEO_STATE ;
   
   PEL_CHECK_POST( geometric_state() != OLD(geometric_state) ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: print" ) ;

   std::string space( indent_width, ' ' ) ;
   os << space << "VERTICES" << endl ;
   VERTS->print( os, indent_width+3 ) ;
   os << space << "SIDES" << endl ;
   for( size_t i=0 ; i<nb_sides() ; ++i )
   {
      PDE_FaceFE* side = static_cast<PDE_FaceFE*>( SIDES->at( i ) ) ;
      os << space << "   ------" << endl ;
      side->print( os, indent_width+3 ) ;
   }
   os << space << "BOUNDS" << endl ;
   for( size_t i=0 ; i<nb_bounds() ; ++i )
   {
      PDE_BoundFE* bound = static_cast<PDE_BoundFE*>( BOUNDS->at( i ) ) ;
      os << space << "   ------" << endl ;
      bound->print( os, indent_width+3 ) ;
   }
   os << space << "CELLS" << endl ;
   for( size_t i=0 ; i<nb_cells() ; ++i )
   {
      PDE_CellFE* cell = static_cast<PDE_CellFE*>( CELLS->at( i ) ) ;
      os << space << "   ------" << endl ;
      cell->print( os, indent_width+3 ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_GridFE:: add_cell( PDE_CellFE* a_cell )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: append_cell" ) ;
   PEL_CHECK_PRE( a_cell != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t ii = a_cell->id_number() ;
   if( ii >= CELLS->index_limit() )
   {
      CELLS->resize( ii+1 ) ;
   }
   PEL_ASSERT( CELLS->at( ii ) == 0 ) ;
   CELLS->set_at( ii, a_cell ) ;

   PEL_CHECK_POST( 
      static_cast< PDE_CellFE* >( cells()->at( a_cell->id_number() ) ) == 
      a_cell ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: add_bound( PDE_BoundFE* bound ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: append_bound" ) ;
   PEL_CHECK_PRE( bound != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   BOUNDS->append( bound ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: add_side( PDE_FaceFE* a_side )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: append_side" ) ;
   PEL_CHECK_PRE( a_side != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t ii = a_side->side_id() ;
   if( ii == SIDES->index_limit() )
   {
      SIDES->resize( ii+1 ) ;
   }
   SIDES->set_at( ii, a_side ) ;

   PEL_CHECK_POST( 
      static_cast< PDE_FaceFE* >( sides()->at( a_side->side_id() ) ) == 
      a_side ) ;
}

//??? may be useless
//??? transmit directly the face instead of its number
//----------------------------------------------------------------------
void
PDE_GridFE:: build_bound( size_t number,
                          size_t side_number,
                          GE_Color const* color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: build_bound" ) ;

   PDE_FaceFE* side = static_cast<PDE_FaceFE*>( FACES->at( side_number ) ) ;
   GE_Mpolyhedron const* side_poly = side->polyhedron() ;

   size_t_vector vertices( side_poly->nb_vertices() ) ;
   for( size_t iv=0 ; iv<side_poly->nb_vertices() ; ++iv )
   {
      GE_Point const* current_vertex = side_poly->vertex( iv ) ;
      vertices( iv ) = VERTS->index( current_vertex ) ;
      PEL_ASSERT( VERTS->point( vertices( iv ) ) == current_vertex ) ;
   }
   GE_Mpolyhedron* poly = GE_Mpolyhedron:: create( side_poly->name(),
                                                   VERTS,
                                                   vertices ) ;
   
   PEL_Object* oo = B_DISCS->at( poly->reference_polyhedron()->id_number() ) ;
//   PEL_ASSERT( oo != 0 ) ; //???????????????????????????????,,
   PDE_DiscOnMeshFE* disc = static_cast< PDE_DiscOnMeshFE* >( oo ) ;
   size_t rlevel = 0 ;
   PDE_BoundFE* bdm = PDE_BoundFE::create( VERTS, number, poly, color, rlevel,
                                           disc ) ;

   add_bound( bdm ) ;

   PEL_ASSERT( side->nb_adjacent_cells() == 1 ) ;
   PDE_CellFE* mesh = side->adjacent_cell( 0 ) ;
      
   bdm->insert_adjacent_cell( mesh ) ;
   side->insert_adjacent_bound( bdm ) ;
}

//-----------------------------------------------------------------------
void
PDE_GridFE:: build_new_vertices( GE_ReferencePolyhedronRefiner const* crf,
                                 PDE_CellFE const* ccell,
                                 size_t_vector& idx_in_VERTS )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: build_new_vertices" ) ;

   bool newvert = false ;
   if( VERB > 1 ) PEL::out() << "   vertices involved :" << endl ;

   GE_Mpolyhedron const* cpoly = ccell->polyhedron() ;

   for( size_t i=0 ; i<crf->nb_vertices() ; ++i )
   {
      cpoly->apply_mapping( crf->vertex( i ), W_PT ) ;
      if( VERB > 1 ) newvert = !VERTS->has( W_PT ) ;
      VERTS->extend( W_PT ) ;
      idx_in_VERTS( i ) = VERTS->index( W_PT ) ;
      if( VERB > 1 )
      {
         PEL::out() << setw( 9 ) << idx_in_VERTS( i ) ;
         W_PT->print( PEL::out(), 2 ) ;
         if( newvert ) PEL::out() << "   NEW" ;
         PEL::out() << endl ;
      }
   }
}

//-----------------------------------------------------------------------
void
PDE_GridFE:: record_already_existing_faces(
                                GE_ReferencePolyhedronRefiner const* crf,
                                PDE_CellFE const* ccell,
                                PEL_Vector* subfaces,
                                PEL_BalancedBinaryTree* tree )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: record_already_existing_faces" ) ;
   
   PDE_FacesOfCellFE* itcs = ccell->create_faces_iterator( 0 ) ;
   for( itcs->start() ; itcs->is_valid() ; itcs->go_next() )
   {
      PDE_FaceFE* face = itcs->item() ;
      if( face->refinement_level() == ccell->refinement_level()+1 )
      {
         if( VERB > 1 )
                  PEL::out() << "   already existing subside of vertices :" ;
         add_index_set( tree, face, subfaces->index_limit() ) ;
         subfaces->append( face ) ;
      }
      else if( face->refinement_level() == ccell->refinement_level() &&
               face->nb_childs() == 0 )
      {
         face->set_nb_childs( crf->nb_subfaces_per_face() ) ;
         if( face->is_periodic() )
         {
            PDE_FaceFE* pface = face->periodic_neighbour() ;
            PEL_ASSERT( pface->nb_childs() == 0 ) ;
            pface->set_nb_childs( crf->nb_subfaces_per_face() ) ;
         }
      }
   }

   itcs->destroy() ;
}

//-----------------------------------------------------------------------
void
PDE_GridFE:: build_new_faces( GE_ReferencePolyhedronRefiner const* crf,
                              PDE_CellFE const* ccell,
                              size_t_vector const& idx_in_VERTS,
                              PEL_BalancedBinaryTree const* tree,
                              PEL_Vector* subfaces,
                              size_t_vector& idx_in_subfaces )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: build_new_faces" ) ;

   if( VERB > 1 ) PEL::out() << "   subsides" << endl ;

   size_t i_side = nb_sides() ; 
   size_t i_face = i_side + nb_bounds() ;
   
   PDE_FaceFE* aside =
               static_cast<PDE_FaceFE*>( ccell->faces()->at( 0 ) ) ;
   for( size_t is=0 ; is<crf->nb_subfaces() ; ++is )
   {
      if( VERB > 1 ) PEL::out() << "      vertices: " ;
      size_t nbvs = crf->subface_reference_polyhedron()->nb_vertices() ;
      size_t_vector idx_verts( nbvs ) ;
      for( size_t iv=0 ; iv<nbvs ; ++iv )
      {
         idx_verts( iv ) = idx_in_VERTS( crf->subface_vertex( is, iv ) ) ;
         if( VERB > 1 ) PEL::out() << setw( 4 ) << idx_verts( iv ) << " " ;
      }
      IDX_TEST->re_initialize( idx_verts, 0 ) ;

      if( tree->has( IDX_TEST ) )
      {
         // the refined face already existed before the refinement of ccell
         PEL_IndexSet const* old_idx_set = 
                      static_cast< PEL_IndexSet* >( tree->item( IDX_TEST ) ) ;
         if( VERB > 1 ) PEL::out() << "                   " ;
         idx_in_subfaces( is ) = old_idx_set->id() ;
      }
      else
      {
         if( VERB > 1 ) PEL::out() << "   NEW " ;
         GE_Mpolyhedron* rpoly =
            GE_Mpolyhedron::create( aside->polyhedron()->name(), VERTS, 
                                    idx_verts ) ;
         PDE_FaceFE* rface = PDE_FaceFE::create( VERTS,
                                                 i_face++,
                                                 rpoly,
                                                 ccell->color(),
                                                 ccell->refinement_level()+1 ) ;
         
         idx_in_subfaces( is ) = subfaces->index_limit() ;
         subfaces->append( rface ) ;

         size_t idx_parent = crf->subface_parent( is ) ;
         if( idx_parent != PEL::bad_index() )
         {
            PDE_FaceFE* ss =
                static_cast<PDE_FaceFE*>( ccell->faces()->at( idx_parent ) ) ;
            if( VERB > 1 ) PEL::out() << " parent:" << setw( 4 )
                                            << ss->id_number() ;
            rface->set_color( ss->color() ) ;
            rface->set_parent( ss ) ;
            ss->append_child( rface ) ;
            
            if( ss->has_adjacent_bound() )
            {
               PDE_BoundFE* cbound = ss->adjacent_bound() ;

               GE_Mpolyhedron* pp = 
                  GE_Mpolyhedron::create( rpoly->name(), VERTS, idx_verts ) ;
               
               //????? ca se répète ??? passer directement le vecteur ???
               PEL_Object* oo = B_DISCS->at( 
                                pp->reference_polyhedron()->id_number() ) ;
               PDE_DiscOnMeshFE* disc = 
                          static_cast< PDE_DiscOnMeshFE* >( oo ) ;
               PDE_BoundFE* rbound =
                  PDE_BoundFE::create( VERTS,
                                       nb_bounds(),
                                       pp,
                                       cbound->color(),
                                       cbound->refinement_level()+1,
                                       disc ) ;
               add_bound( rbound ) ;
               rface->insert_adjacent_bound( rbound ) ;
               rbound->set_parent( cbound ) ;
            }
            else
            {
               rface->set_side_id( i_side++ ) ;
               add_side( rface ) ;               
               if( ss->is_periodic() )
               {
                  PDE_FaceFE* pss = ss->periodic_neighbour() ; //??? pss ?
                  GE_Transform const* tr = ss->periodic_transform() ;
                  // création de la sousface périodique
                  for( size_t iv=0 ; iv<rpoly->nb_vertices() ; ++iv )
                  {
                     W_PT->set( rpoly->vertex( iv ) ) ;
                     tr->apply( W_PT ) ;
                     VERTS->extend( W_PT ) ;
                     idx_verts( iv ) = VERTS->index( W_PT ) ;
                  }
                  GE_Mpolyhedron* pp = 
                        GE_Mpolyhedron::create( rpoly->name(), VERTS, idx_verts ) ;
                  PDE_FaceFE* prface = 
                           PDE_FaceFE::create( VERTS,
                                               i_face++,
                                               pp,
                                               pss->color(),
                                               rface->refinement_level() ) ;
                  prface->set_parent( pss ) ;
                  pss->append_child( prface ) ;
                  
                  prface->insert_adjacent_cell( pss->adjacent_cell( 0 ) ) ;
                     
                  rface->set_periodicity( prface, tr ) ;
                  prface->set_periodicity( rface, tr->inverse() ) ;
               }
            }
         }
         else
         {
            if( VERB > 1 ) PEL::out() << setw( 12 ) << "no parent" ;
            rface->set_side_id( i_side++ ) ;
            add_side( rface ) ;
         }
      }
      if( VERB > 1 )
      {
         PEL::out() << "   id:" << setw( 4 ) ;
         PDE_FaceFE const* face = static_cast< PDE_FaceFE* >( 
                                  subfaces->at( idx_in_subfaces( is ) ) ) ;
         PEL::out() << face->id_number() << endl ;
      }
   }
}

//-----------------------------------------------------------------------
void
PDE_GridFE:: build_new_cells( GE_ReferencePolyhedronRefiner const* crf,
                              PDE_CellFE* ccell,
                              size_t_vector const& idx_in_VERTS,
                              PEL_Vector const* subfaces,
                              size_t old_subfaces_limit,
                              size_t_vector const& idx_in_subfaces )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: build_new_cells" ) ;

   if( VERB > 1 ) PEL::out() << "   subcells" << endl ;

   GE_Mpolyhedron const* cpoly = ccell->polyhedron() ;

   ccell->set_nb_childs( crf->nb_subcells() ) ;
   for( size_t ic=0 ; ic<crf->nb_subcells() ; ++ic )
   {
      if( VERB > 1 ) PEL::out() << "      vertices: " ;
      size_t nbvs = crf->subcell_reference_polyhedron()->nb_vertices() ;
      size_t_vector idx( nbvs ) ;
      for( size_t iv=0 ; iv<nbvs ; ++iv )
      {
         idx( iv ) = idx_in_VERTS( crf->subcell_vertex( ic, iv ) ) ;
         if( VERB > 1 ) PEL::out() << setw( 4 ) << idx( iv ) << " " ;
      }
      GE_Mpolyhedron* rpoly =
         GE_Mpolyhedron::create( cpoly->name(), VERTS, idx ) ;

      PDE_CellFE* rcell = PDE_CellFE::create_with_discretizations_pattern(
                                              VERTS,
                                              nb_cells(),
                                              rpoly,
                                              ccell->color(),
                                              ccell->refinement_level()+1,
                                              ccell ) ;

      if( VERB > 1 ) PEL::out() << "     id:" << rcell->id_number() << endl ;

      ccell->set_child( ic, rcell ) ;
      rcell->set_parent( ccell ) ;
      add_cell( rcell ) ;

      //-- lien avec les sides raffinées
      size_t nbss = crf->subcell_reference_polyhedron()->nb_faces() ;
      for( size_t is=0 ; is<nbss ; ++is )
      {
         size_t iii = idx_in_subfaces( crf->subcell_face( ic, is ) ) ;
         PDE_FaceFE* face = static_cast< PDE_FaceFE* >( subfaces->at( iii ) ) ;
         if( VERB > 1 )
            PEL::out() << "         subface id_number: " << face->id_number() ;
         rcell->append_face( face ) ;
         if( iii < old_subfaces_limit )
         {
            // refined face already existing before the refinement of ccell
            //    => belongs to the boundary of ccell
            //    => does not belong to the boundary of the domain
            if( VERB > 1 ) PEL::out() << "   side_id: " << endl ;
            PEL_ASSERT( !face->has_adjacent_bound() ) ;
            PDE_FaceFE* prs = face->parent() ;
            if( prs!=0 && prs->is_active() )
            {
               PEL_ASSERT( !face->is_active() ) ;
               PDE_CellFE* ocell = prs->adjacent_cell_other_than( ccell ) ;
               if( VERB > 1 )
                  PEL::out() << "            adj cell: " << ocell->id_number()
                             << "replaced by cell: " << rcell->id_number()
                             << endl ;
               face->set_adjacent_cells( rcell, ocell ) ;
            }
            else
            {
               if( VERB > 1 )
                  PEL::out() << "            set adj cell: "
                             << ccell->id_number() << "  cell: "
                             << rcell->id_number() << endl ;
               replace_adjacent_cell( face, ccell, rcell ) ;
            }
         }
         else
         {
            // refined face newly created by the refinement of ccell
            if( VERB > 1 ) PEL::out() << "   NEW" << endl ;
            if( face->has_adjacent_cell( ccell ) )
            {
               face->replace_adjacent_cell( ccell, rcell ) ;
            }
            else
            {
               face->insert_adjacent_cell( rcell ) ;
            }
            if( face->nb_adjacent_cells() == 1 )
            {
               PDE_FaceFE* prs = face->parent() ;
               if( prs != 0 && !prs->has_adjacent_bound() )
               {
                  // the refined face is on the boundary of ccell
                  if( !prs->is_periodic() )
                  {
                     PDE_CellFE* ocell = prs->adjacent_cell_other_than( ccell ) ;
                     face->insert_adjacent_cell( ocell ) ;
                     if( ocell->parent() != 0 )
                     {
                        PEL_ASSERT( ocell->parent() != ccell ) ;
                     }
                  }
                  else
                  {
                     //???? necessairement, face est periodique
                     //???? mais avec qui ???
                  }
               }
            }
            
            if( face->has_adjacent_bound() )
            {
               PDE_BoundFE* bound = face->adjacent_bound() ;
               bound->insert_adjacent_cell( face->adjacent_cell( 0 ) ) ;
            }
         }
      }
   }
}

//-----------------------------------------------------------------------
void
PDE_GridFE:: replace_adjacent_cell( PDE_FaceFE* face,
                                    PDE_CellFE* old_cell,
                                    PDE_CellFE* new_cell )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: replace_adjacent_cell" ) ;

   if( face->is_active() )
   {
      face->replace_adjacent_cell( old_cell, new_cell ) ;
   }
   else
   {
      PEL_ASSERT( face->nb_childs() != 0 ) ;
      for( size_t i=0 ; i<face->nb_childs() ; ++i )
      {
         replace_adjacent_cell( face->child( i ), old_cell, new_cell ) ;
      }
   }
}

//-----------------------------------------------------------------------
void
PDE_GridFE:: add_index_set( PEL_BalancedBinaryTree* tree,
                            PDE_MeshFE const* a_mesh,
                            size_t an_id )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: add_index_set" ) ;
   
   GE_Mpolyhedron const* poly = a_mesh->polyhedron() ;
   size_t_vector idx( poly->nb_vertices() ) ;
   for( size_t iv=0 ; iv<poly->nb_vertices() ; ++iv )
   {
      GE_Point const* pt = poly->vertex( iv ) ;
      PEL_ASSERT( VERTS->has( pt ) ) ;
      idx( iv ) = VERTS->index( pt ) ;
      if( VERB > 1 ) PEL::out() << setw( 4 ) << idx( iv ) ;
   }
   if( VERB > 1 ) PEL::out() << endl ;
   PEL_IndexSet* idx_set = PEL_IndexSet::create( tree, idx, an_id ) ;
   PEL_ASSERT( !tree->has( idx_set ) ) ;
   tree->extend( idx_set ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: separate_sides_and_bounds( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: separate_sides_and_bounds" ) ;
   
   PEL_ASSERT( TRANSFOS_IT == 0 ) ;
   
   size_t is = 0 ;
   size_t ib = 0 ;
   PEL_VectorIterator* its = PEL_VectorIterator::create( 0, FACES ) ;
   for( its->start(); its->is_valid(); its->go_next() )
   {
      PDE_FaceFE* side = static_cast<PDE_FaceFE*>( its->item() ) ;
      if( side->nb_adjacent_cells() == 1 )
      {
         build_bound( ib++, side->id_number(), side->color() ) ;
      }
      else
      {
         side->set_side_id( is ) ;
         add_side( side ) ;
         is++ ;
      }
   }
   its->destroy() ; its = 0 ;

   // to end the building process
   // may be formalize in "start_building" and "end_building" methods
   FACES->destroy() ; FACES = 0 ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: separate_sides_and_bounds( PEL_Iterator* it_tr )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: separate_sides_and_bounds" ) ;
   
   size_t rank = PEL_Exec::communicator()->rank() ;
   size_t nb_ranks = PEL_Exec::communicator()->nb_ranks() ;
   
   PEL_ListIdentity* source_cols = PEL_ListIdentity::create( 0 ) ;
   PEL_ListIdentity* target_cols = PEL_ListIdentity::create( 0 ) ;
   for( it_tr->start() ; it_tr->is_valid() ; it_tr->go_next() )
   {
      GE_Transform const* tr = static_cast< GE_Transform* >( it_tr->item() ) ;
      GE_Color const* c1 = tr->source_color() ;
      if( source_cols->has( c1 ) || target_cols->has( c1 ) )
      {
         PDE_GridFE_ERROR::n0( c1 ) ;
      }
      source_cols->append( const_cast< GE_Color* >( c1 ) ) ;
      
      GE_Color const* c2 = tr->target_color() ;
      PEL_ASSERT( !source_cols->has( c2 ) && !target_cols->has( c2 ) ) ;
      target_cols->append( const_cast< GE_Color* >( c2 ) ) ;
   }
   
   GE_Point* pt_tr =  GE_Point::create( 0, nb_space_dimensions() ) ;
   PEL_BalancedBinaryTree* tree = PEL_BalancedBinaryTree::create( 0 ) ;
   
   size_t nb_faces = FACES->count() ;
   
   boolVector is_bound( nb_faces ) ;
   is_bound.set( false ) ;
   size_t is = 0 ;
   for( size_t i=0 ; i < nb_faces ; ++i )
   {
      PDE_FaceFE* face = static_cast<PDE_FaceFE*>( FACES->at( i ) ) ;
      PEL_ASSERT( face != 0 ) ;
      if( face->nb_adjacent_cells() == 1 )
      {
         add_index_set( tree, face, i ) ;
         is_bound( i ) = true ;
      }
      else
      {
         GE_Color const* col = face->color() ;
         if( source_cols->has( col ) || target_cols->has( col ) )
         {
            PDE_GridFE_ERROR::n1( face, it_tr ) ;
         }
         face->set_side_id( is ) ;
         add_side( face ) ;
         is++ ;         
      }
   }
   
   size_t_vector pair_1( 0 ) ;
   size_t_vector pair_2( 0 ) ;
   PEL_Vector* tr_12 = PEL_Vector::create( 0, 0 ) ;
   PEL_BalancedBinaryTreeIterator* it_bd = tree->create_iterator( 0 ) ;
   for( it_bd->start() ; it_bd->is_valid() ; it_bd->go_next() )
   {
      PEL_IndexSet const* ss = static_cast< PEL_IndexSet* >( it_bd->item() ) ;
      size_t_vector const& verts = ss->elements() ;
      for( it_tr->start() ; it_tr->is_valid() ; it_tr->go_next() )
      {
         GE_Transform const* tr = static_cast< GE_Transform* >( it_tr->item() ) ;
         size_t_vector idx( verts.size() ) ;
         bool associated = true ;
         for( size_t iv=0 ; iv<verts.size() ; ++iv )
         {
            pt_tr->set( VERTS->point( verts( iv ) ) ) ;
            tr->apply( pt_tr ) ;
            if( VERTS->has( pt_tr ) )
            {
               idx( iv ) = VERTS->index( pt_tr ) ;
            }
            else
            {
               associated = false ;
               break ;
            }
         }
         if( associated )
         {
            size_t i1 = ss->id() ;
            PDE_FaceFE* face1 = static_cast<PDE_FaceFE*>( FACES->at( i1 ) ) ;
            
            IDX_TEST->re_initialize( idx, 0 ) ;
            if( !tree->has( IDX_TEST ) )
            {
               PDE_GridFE_ERROR::n3( face1, tr ) ;
            }
            PEL_IndexSet const* ss_tr = 
                static_cast< PEL_IndexSet* >( tree->item( IDX_TEST ) ) ;
            
            size_t i2 = ss_tr->id() ;
            PDE_FaceFE* face2 = static_cast<PDE_FaceFE*>( FACES->at( i2 ) ) ;
            if( ( face1->color() != tr->source_color() ) ||
                ( face2->color() != tr->target_color() ) )
            {
               PDE_GridFE_ERROR::n2( face1, face2, tr ) ;
            }
            
            PEL_ASSERT( is_bound( i1 ) && is_bound( i2 ) ) ;
            is_bound( i1 ) = false ;
            is_bound( i2 ) = false ;
            
            if( pair_1.has( i1 ) || pair_2.has( i1 ) )
            {
               PDE_GridFE_ERROR::n4( face1 ) ;               
            }
            if( pair_1.has( i2 ) || pair_1.has( i2 ) )
            {
               PDE_GridFE_ERROR::n4( face2 ) ;               
            }
            pair_1.append( i1 ) ;
            pair_2.append( i2 ) ;
            tr_12->append( const_cast< GE_Transform* >( tr ) ) ;
         }
      }
   }
   it_bd->destroy() ; it_bd = 0 ;
   
   for( size_t i=0 ; i<pair_1.size() ; ++i )
   {
      PDE_FaceFE* side1 = static_cast<PDE_FaceFE*>( FACES->at( pair_1( i ) ) ) ;
      PDE_FaceFE* side2 = static_cast<PDE_FaceFE*>( FACES->at( pair_2( i ) ) ) ;
      GE_Transform const* tr = static_cast<GE_Transform*>( tr_12->at( i ) ) ;
      side1->set_side_id( is ) ;
      side1->set_periodicity( side2, tr ) ;
      add_side( side1 ) ;
      side2->set_periodicity( side1, tr->inverse() ) ;
      if( side1->adjacent_cell(0)->color() == GE_Color::halo_color() )
      {
         side1->set_color( GE_Color::halo_color() ) ;
         side2->set_color( GE_Color::halo_color() ) ;
      }
      is++ ;         
   }
      
   size_t ib = 0 ;
   for( size_t i=0 ; i<nb_faces ; ++i )
   {
      PDE_FaceFE* face = static_cast<PDE_FaceFE*>( FACES->at( i ) ) ;
      GE_Color const* fcolor = face->color() ;
      if( is_bound( i ) )
      {
         if( source_cols->has( fcolor ) || target_cols->has( fcolor ) )
         {
            PDE_GridFE_ERROR::n5( face, it_tr, rank, nb_ranks ) ;
         }
         build_bound( ib++, face->id_number(), face->color() ) ;
      }
      if( source_cols->has( fcolor ) || target_cols->has( fcolor ) )
      {
         if( !face->is_periodic() )
         {
            PDE_GridFE_ERROR::n5( face, it_tr, rank, nb_ranks ) ;
         }
      }
   }
   
   tr_12->destroy() ;
   pt_tr->destroy() ;
   tree->destroy() ;
   source_cols->destroy() ;
   target_cols->destroy() ;
   
   // to end the building process
   // may be formalize in "start_building" and "end_building" methods
   FACES->destroy() ; FACES = 0 ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: build_cluster_to_mesh_connectivity( PEL_Vector const* meshes )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: build_cluster_to_mesh_connectivity" ) ;
   
   CLUSTER_OF_MESH.re_initialize( meshes->count() ) ;
   CLUSTER_OF_MESH.set( PEL::bad_index() ) ;

   NB_MESHES = meshes->count() ;
   if( MESH_IT != 0 ) destroy_possession( MESH_IT ) ;
   MESH_IT = meshes->create_iterator( this ) ;
   for( MESH_IT->start() ; MESH_IT->is_valid() ; MESH_IT->go_next() )
   {
      PDE_MeshFE const* m = static_cast<PDE_MeshFE*>( MESH_IT->item() ) ;
      GE_Color const* mcol = m->color() ;

      size_t icl = 0 ;
      if( NB_CLUSTERS != 1 )
      {
         bool found = false ;
         for( size_t i=0 ; i<NB_CLUSTERS ; ++i )
         {
            if( CLUSTER_COLOR[i]->is_matching( mcol ) )
            {
               if( found )
                  PDE_GridFE_ERROR::n8( mcol->name(),
                                         CLUSTER_COLOR[icl]->name(),
                                         CLUSTER_COLOR[i]->name() ) ;
               found = true ;
               icl = i ;
            }
         }
         if( !found ) PDE_GridFE_ERROR::n7( mcol->name() ) ;
      }
      CLUSTER_OF_MESH( m->id_number() ) = icl ;
   }
}

//----------------------------------------------------------------------
size_t
PDE_GridFE:: associated_point( GE_SetOfPoints const* points,
                               PEL_ListIterator* tr_it,
                               size_t ip,
                               GE_Point* work_pt )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: associated_point" ) ;

   size_t result = PEL::bad_index() ;

   GE_Point const* pt = points->point( ip ) ;
   for( tr_it->start() ; tr_it->is_valid() ; tr_it->go_next() )
   {
      GE_Transform const* tr = static_cast< GE_Transform* >( tr_it->item() ) ;
      work_pt->set( pt ) ;
      tr->apply( work_pt ) ;
      if( points->has( work_pt ) )
      {
         result = points->index( work_pt ) ;
         break ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_GridFE:: fill_pt_2_node(
          size_t ip, size_t_vector const& ass,
          size_t_vector& pt_2_node, size_t& current_node )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: fill_pt_2_node" ) ;

   if( ass( ip ) == PEL::bad_index() )
   {
      if( pt_2_node( ip ) == PEL::bad_index() )
      {
         pt_2_node( ip ) = current_node ;
         ++current_node ;
      }
   }
   else
   {
      if( pt_2_node( ass( ip ) ) == PEL::bad_index() )
      {
         fill_pt_2_node( ass( ip ), ass,
                         pt_2_node, current_node ) ;
      }
      pt_2_node( ip ) = pt_2_node( ass( ip ) ) ;
   }
}

//----------------------------------------------------------------------
bool
PDE_GridFE:: rejected_mesh( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: rejected_mesh" ) ;
   
   bool result = ( CLUSTER_OF_MESH( the_mesh( MESH_IT )->id_number() ) != 
                   I_CURRENT_CLUSTER ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_MeshFE*
PDE_GridFE:: the_mesh( PEL_VectorIterator const* it )
//---------------------------------------------------------------------------
{
   return( static_cast<PDE_MeshFE*>( it->item() ) ) ;
}

//----------------------------------------------------------------------
PDE_CellFE*
PDE_GridFE:: get_mesh( GE_Point const* pt,
                       GE_SegmentSegment_INT const* seg_seg_intersector,
                       PDE_CellFE* guess_mesh ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridFE:: get_mesh" ) ;
   PEL_CHECK_PRE( nb_space_dimensions()==2 ) ;
   PEL_CHECK_PRE( pt!=0 && pt->nb_coordinates()==2 ) ;
   PEL_CHECK_PRE( seg_seg_intersector!=0 ) ;
   PEL_CHECK_PRE( IMPLIES( guess_mesh!=0,
                           guess_mesh->polyhedron()->dimension()==2 &&
                           guess_mesh->polyhedron()->nb_space_dimensions()==2 ) ) ;
   
   if( guess_mesh==0 )
   {
      guess_mesh = static_cast<PDE_CellFE*>( CELLS->at( 0 ) ) ;
   }

   PDE_CellFE* result = guess_mesh ;
   
   if( !result->polyhedron()->contains( pt ) )
   {

      // Standard search
      result = search_mesh( pt, seg_seg_intersector, guess_mesh, nb_cells() ) ;
      if( result==0 ) 
      {
         // Search failed : initialization with the closest bounds
         double d = PEL::max_double() ;
         PDE_CellFE* guess_mesh1 = 0 ;
         PDE_CellFE* guess_mesh2 = 0 ;
         PEL_VectorIterator* bd_it = PEL_VectorIterator::create( 0, BOUNDS ) ;
         for( bd_it->start() ; bd_it->is_valid() ; bd_it->go_next() )
         {
            PDE_BoundFE const* bound = 
                               static_cast<PDE_BoundFE*>( bd_it->item() ) ;
            // BUG : FAULT=00129
            // double dd = pt->distance( bound->polyhedron()->center() ) ;
            double dd =
               PEL::min( pt->distance( bound->polyhedron()->vertex(0) ),
                         pt->distance( bound->polyhedron()->vertex(1) ) ) ;
            if( PEL::toler( PEL::abs( dd-d ) ) )
            {
               guess_mesh2 = bound->adjacent_cell() ;
            }
            else if( dd<d )
            {
               d = dd ;
               guess_mesh1 = bound->adjacent_cell() ;
               guess_mesh2 = 0 ;
            }
         }
         bd_it->destroy() ; bd_it = 0 ;

         // Search :
         result = search_mesh( pt, seg_seg_intersector, guess_mesh1, nb_cells() ) ;
         if( result==0 && guess_mesh2!=0 )
         {
            result = search_mesh( pt, seg_seg_intersector, guess_mesh2, nb_cells() ) ;
         }
      }
   }

   PEL_CHECK_POST(
      IMPLIES( guess_mesh!=0 && guess_mesh->polyhedron()->contains( pt ),
               result==guess_mesh ) ) ;
   PEL_CHECK_POST(
      IMPLIES( result!=0, result->polyhedron()->contains( pt ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_CellFE*
PDE_GridFE:: search_mesh(
                       GE_Point const* pt,
                       GE_SegmentSegment_INT const* seg_seg_intersector,
                       PDE_CellFE* guess_mesh,
                       size_t recursive_max ) const
//----------------------------------------------------------------------
{
   PEL_CHECK( nb_space_dimensions()==2 ) ;
   PEL_CHECK( pt!=0 && pt->nb_coordinates()==2 ) ;
   PEL_CHECK( seg_seg_intersector!=0 ) ;
   PEL_CHECK( guess_mesh!=0 &&
              guess_mesh->polyhedron()->dimension()==2 &&
              guess_mesh->polyhedron()->nb_space_dimensions()==2 ) ;

   PDE_CellFE * result = 0 ;
   
   if( recursive_max!=0 )
   {
      result = guess_mesh ;
      if( !result->polyhedron()->contains( pt ) )
      {
         // Looking for the side which crosses the segment delimited by
         // pt and mesh barycenter
         PDE_FaceFE const* side = 0 ;
         PEL_VectorIterator* itSide = 
                  PEL_VectorIterator::create( 0, result->faces() ) ;
         while( itSide->is_valid() && side==0 )
         {
            PDE_FaceFE const* ss = static_cast<PDE_FaceFE*>( itSide->item() ) ;
            if( seg_seg_intersector->has_intersection(
               ss->polyhedron()->vertex(0),
               ss->polyhedron()->vertex(1),
               result->polyhedron()->center(),
               pt ) )
            {
               side = ss ;
            }
            else
            {
               itSide->go_next() ;
            }
         }
         itSide->destroy() ; itSide = 0 ;
         PEL_ASSERT( side!=0 ) ;
 
         if( side->has_adjacent_bound() )
         {
            result = 0 ;
         }
         else
         {
            // Find the mesh linked to this side (different from the initial mesh)
            PDE_CellFE* other_mesh = side->adjacent_cell_other_than( result ) ;
            result = search_mesh( pt, seg_seg_intersector, other_mesh,
                                  recursive_max-1 ) ;
         }
      }
   }

   PEL_CHECK( IMPLIES( guess_mesh->polyhedron()->contains( pt ),
                       result==guess_mesh ) ) ;
   PEL_CHECK( IMPLIES( result!=0, result->polyhedron()->contains( pt ) ) ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
PDE_GridFE_ERROR:: n0( GE_Color const* color )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "   the color \"" << color->name() << "\" is" << endl
        << "   involved in more than one periodicity request"  ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_GridFE_ERROR:: n1( PDE_FaceFE const* face, PEL_Iterator* it_tr )
//internal--------------------------------------------------------------
{
   GE_Transform const* tr = 0 ;
   GE_Color const* fcol = face->color() ;
   bool found = false ;
   for( it_tr->start() ; it_tr->is_valid() ; it_tr->go_next() )
   {
      tr = static_cast< GE_Transform* >( it_tr->item() ) ;
      if( ( tr->source_color() == fcol ) || ( tr->target_color() == fcol ) ) 
      {
         found = true ;
         break ;
      }
   }
   PEL_ASSERT( found ) ;
   
   ostringstream mesg ;
   mesg << "*** Domain periodicity error." << endl ;
   mesg << "    The following face" << endl ;
   face->polyhedron()->print( mesg, 7 ) ;
   mesg << "    has the color \"" << fcol->name() << "\", so it is" << endl ;
   mesg << "    involved in the periodicity request by" << endl ;
   mesg << "    the transformation" << endl ;
   tr->print( mesg, 7 ) ;
   mesg << "    but this face does not belong to the domain boundary." ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_GridFE_ERROR:: n2( PDE_FaceFE const* face1, PDE_FaceFE const* face2,
                              GE_Transform const* tr )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** Domain periodicity error." << endl ;
   mesg << "    It has been found that the transformation" << endl ;
   tr->print( mesg, 7 ) ;
   mesg << "    associates the following face of color \"" 
        << face1->color()->name() << "\"" << endl ;
   face1->polyhedron()->print( mesg, 7 ) ;
   mesg << "    to the following face of color \"" 
        << face2->color()->name() << "\"" << endl ;
   face2->polyhedron()->print( mesg, 7 ) ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_GridFE_ERROR:: n3( PDE_FaceFE const* face1,
                              GE_Transform const* tr )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** Domain periodicity error." << endl ;
   mesg << "    It has been found that the transformation" << endl ;
   tr->print( mesg, 7 ) ;
   mesg << "    associates the following face" << endl ;
   face1->polyhedron()->print( mesg, 7 ) ;
   mesg << "    to a face that is not on the domain boundary." ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_GridFE_ERROR:: n4( PDE_FaceFE const* face )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** Domain periodicity error." << endl ;
   mesg << "    The following face" << endl ;
   face->polyhedron()->print( mesg, 7 ) ;
   mesg << "    is involved in more than one periodic association." ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_GridFE_ERROR:: n5( PDE_FaceFE const* face, 
                              PEL_Iterator* it_tr,
                              size_t rank,
                              size_t nb_ranks )
//internal--------------------------------------------------------------
{
   GE_Transform const* tr = 0 ;
   GE_Color const* fcol = face->color() ;
   bool found = false ;
   for( it_tr->start() ; it_tr->is_valid() ; it_tr->go_next() )
   {
      tr = static_cast< GE_Transform* >( it_tr->item() ) ;
      if( ( tr->source_color() == fcol ) || ( tr->target_color() == fcol ) ) 
      {
         found = true ;
         break ;
      }
   }
   PEL_ASSERT( found ) ;
   
   ostringstream mesg ;
   mesg << "*** Domain periodicity error." << endl ;
   mesg << "    The following face" ;
   if( nb_ranks > 1 )
   {
      mesg << " of the process " << rank << endl ;
   }
   else
   {
      mesg << endl ;
   }
   face->polyhedron()->print( mesg, 7 ) ;
   mesg << "    has the color \"" << fcol->name() << "\"" << endl ;
   mesg << "    so it should be associated by the transformation" << endl ;
   tr->print( mesg, 7 ) ;
   if( nb_ranks > 1 )
   {
      mesg << "    to another face of the process " << rank << endl ;
      mesg << "    but such is not the case." ;      
   }
   else
   {
      mesg << "    to another face, but such is not the case." ;
   }
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_GridFE_ERROR:: n6( void )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** Domain periodicity error." << endl ;
   mesg << "    No valid transformation has been found in" << endl ;
   mesg << "    MODULE periodic_domain." ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_GridFE_ERROR:: n7( std::string const& color_name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "module \"decohesion \" :" << endl
        << "   cells of color \"" << color_name
        << "\" are not handled" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
PDE_GridFE_ERROR:: n8( std::string const& color_name,
                              std::string const& col1,
                              std::string const& col2 )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "module \"decohesion \" :" << endl
        << "   cells of color \"" << color_name << "\" are found to be" << endl
        << "   in the aggregate of color \"" << col1 << "\" and" << endl
        << "   in the aggregate of color \"" << col2 << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

