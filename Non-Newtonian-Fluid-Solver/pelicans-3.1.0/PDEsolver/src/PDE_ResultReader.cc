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

#include <PDE_ResultReader.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_ReferenceElement.hh>

#include <PEL.hh>
#include <PEL_Context.hh>
#include <PEL_DataOnMeshingReader.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>
#include <doubleArray2D.hh>
#include <intArray2D.hh>

#include <iostream>
#include <sstream>

struct PDE_ResultReader_ERROR
{
   static void n0( std::string const& mesg ) ;
} ;

//----------------------------------------------------------------------
PDE_ResultReader*
PDE_ResultReader:: create( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: create" ) ;
   PEL_CHECK_PRE( exp!=0 ) ;

   PEL_ModuleExplorer const* e =
      exp->create_subexplorer( 0, "PEL_DataOnMeshingReader" ) ;
   PEL_DataOnMeshingReader* data_reader =
      PEL_DataOnMeshingReader::make( 0,
                                     e->string_data( "concrete_name" ),
                                     e ) ;
   e->destroy() ; e = 0 ;

   // Check data reader:
   if( data_reader->meshing()==0 )
   {
      PDE_ResultReader_ERROR::n0( "No stored meshing found" ) ;
   }
   if( data_reader->fields()==0 )
   {
      PDE_ResultReader_ERROR::n0( "No stored fields found" ) ;
   }
   
   PDE_ResultReader* result = new PDE_ResultReader( a_owner,
                                                    exp, data_reader ) ;

   data_reader->set_owner( result ) ;

   PEL_CHECK_POST( result!=0 && result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_ResultReader:: PDE_ResultReader(
                           PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp,
                           PEL_DataOnMeshingReader const* data_reader )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , DEBUG( exp->has_entry( "debug" ) ? exp->bool_data( "debug" ) : false )
   , READER( data_reader )
   , NB_SP_DIMS( data_reader->meshing()->
                       doubleArray2D_data( "vertices" ).index_bound(0) )
   , NB_VERTS( data_reader->meshing()->
                       doubleArray2D_data( "vertices" ).index_bound(1) )
   , VERTS( data_reader->meshing()->doubleArray2D_data( "vertices" ) )
   , NB_CELLS( data_reader->meshing()->
                       intArray2D_data( "cell2vertex" ).index_bound(1) )
   , CELL2VERT( data_reader->meshing()->intArray2D_data( "cell2vertex" ) )
   , FIELDS( 0 )
   , FNAMES( 0 )
   , COMP( 0 )
   , VERT2CELL( 0, 0 )
   , NB_CELL_PER_VERT(0)
   , EPS_VERT( 0. )
   , CLOSEST_VERT(0)
   , VISITED_VERT( 0 )
   , VERT( 0 )
   , POLY( 0 )
   , POLY_VERTS( PEL_Vector::create( this, 0 ) )
   , REF_ELM( 0 )
   , PT( 0 )
   , PT_REF( 0 )
{
   PEL_LABEL( "PDE_ResultReader:: PDE_ResultReader" ) ;

   // Grid connectivity :   
   size_t nb_verts_per_cell = CELL2VERT.index_bound( 0 ) ;
   std::string poly_type = "" ;
   std::string elem_type = "" ;
   if( NB_SP_DIMS==0 && nb_verts_per_cell==1 )
   {
      poly_type = "GE_Mpoint" ;
      elem_type = "PDE_0D_Q0_1node" ;
   }
   else if( NB_SP_DIMS==1 && nb_verts_per_cell==2 )
   {
      poly_type = "GE_Segment" ;
      elem_type = "PDE_1D_P1_2nodes" ;
   }
   else if( NB_SP_DIMS==2 && nb_verts_per_cell==3 )
   {
      poly_type = "GE_Triangle" ;
      elem_type = "PDE_2D_P1_3nodes" ;
   }
   else if( NB_SP_DIMS==2 && nb_verts_per_cell==4 )
   {
      poly_type = "GE_Quadrilateral" ;
      elem_type = "PDE_2D_Q1_4nodes" ;
   }
   else if( NB_SP_DIMS==3 && nb_verts_per_cell==8 )
   {
      poly_type = "GE_Cuboid" ;
      elem_type = "PDE_3D_Q1_8nodes" ;
   }
   else if( NB_SP_DIMS==3 && nb_verts_per_cell==4 )
   {
      poly_type = "GE_Tetrahedron" ;
      elem_type = "PDE_3D_P1_4nodes" ;
   }
   else
   {
      PDE_ResultReader_ERROR::n0( "inconsistent cells" ) ;
   }
   if( exp->has_entry( "polyhedron_type" ) )
   {
      poly_type = exp->string_data( "polyhedron_type" ) ;
   }
   VERT.re_initialize( NB_SP_DIMS ) ;
   VISITED_VERT.re_initialize( NB_VERTS ) ;
   NB_CELL_PER_VERT.re_initialize( NB_VERTS ) ;
   NB_CELL_PER_VERT.set(0) ;
   for( size_t i=0 ; i<NB_CELLS ; ++i )
   {
      for( size_t j=0 ; j<CELL2VERT.index_bound(0) ; ++j )
      {
         NB_CELL_PER_VERT( CELL2VERT(j,i) )++ ;
      }
   }
   size_t nb_m = 0 ;
   for( size_t i=0 ; i<NB_VERTS ; ++i )
   {
      nb_m = PEL::max( nb_m, NB_CELL_PER_VERT(i) ) ;
   }
   VERT2CELL.re_initialize( NB_VERTS, nb_m ) ;
   NB_CELL_PER_VERT.set(0) ;
   for( size_t i=0 ; i<NB_CELLS ; ++i )
   {
      for( size_t j=0 ; j<CELL2VERT.index_bound(0) ; ++j )
      {
         size_t const i_vert = CELL2VERT(j,i) ;
         VERT2CELL( i_vert, NB_CELL_PER_VERT( i_vert ) ) = i ;
         NB_CELL_PER_VERT( i_vert )++ ;
      }
   }
   // Poly :
   for( size_t j=0 ; j<CELL2VERT.index_bound(0) ; ++j )
   {
      size_t i_vert = CELL2VERT(j,0) ;
      VERTS.extract_section( 1, i_vert, VERT ) ;
      GE_Point* p_pt = GE_Point::create( POLY_VERTS, VERT ) ;
      POLY_VERTS->append( p_pt ) ;
   }
   POLY = GE_Mpolyhedron::create( this, poly_type, POLY_VERTS ) ;
   REF_ELM = PDE_ReferenceElement::object( elem_type ) ;
   PT = GE_Point::create( this, NB_SP_DIMS ) ;
   PT_REF = GE_Point::create( this, NB_SP_DIMS ) ;
   if( POLY->nb_vertices()!=CELL2VERT.index_bound(0) &&
       exp->has_entry( "polyhedron_type" ) )
   {
      PEL_Error::object()->raise_data_error(
         exp,
         "polyhedron_type",
         "this definition polyhedron is not compatible\""
         "with the stored meshing" ) ;
   }
   if( POLY->nb_vertices()!=CELL2VERT.index_bound(0) ||
       REF_ELM->nb_nodes()!=CELL2VERT.index_bound(0) )
   {
      PEL_Error::object()->raise_internal( "inconsistent polyhedron" ) ;
   }

   // Meshing superimposition :

   if( exp->has_module( "meshing_superimposition" ) )
   {

      PEL_ModuleExplorer* mesh_exp 
                = exp->create_subexplorer( 0, "meshing_superimposition" ) ;   
      EPS_VERT =  mesh_exp->double_data( "tolerance_for_vertices" ) ;
      mesh_exp->destroy() ; mesh_exp = 0 ;
   }
   else 
   {
      EPS_VERT = 1.E-7 ;
   }

   // Fields :
   PEL_ModuleExplorer* fields_exp = exp->create_subexplorer( 0, "fields" ) ;
   for( fields_exp->start_module_iterator() ;
        fields_exp->is_valid_module() ;
        fields_exp->go_next_module() )
   {
      PEL_ModuleExplorer* f_exp = fields_exp->create_subexplorer( 0 ) ;
      std::string const& field_name = f_exp->string_data( "field" ) ;
      std::string const& storage_name = f_exp->string_data( "entry_name" ) ;
      if( READER->fields()->has_module( storage_name ) )
      {
         FNAMES.append( storage_name ) ;
         FIELDS.append( field_name ) ;
         if( f_exp->has_entry( "component" ) )
         {
            COMP.append( (size_t) f_exp->int_data( "component" ) ) ;
         }
         else
         {
            COMP.append( PEL::bad_index() ) ;
         }
      }
      f_exp->destroy() ; f_exp = 0 ;
   }
   fields_exp->destroy() ; fields_exp = 0 ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_ResultReader:: ~PDE_ResultReader( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: ~PDE_ResultReader" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultReader:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   return( NB_SP_DIMS ) ;
}

//----------------------------------------------------------------------
bool
PDE_ResultReader:: is_in_grid( GE_Point const* pt ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: is_in_grid" ) ;
   PEL_CHECK_PRE( pt!=0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates()==nb_space_dimensions() ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = false ;
   size_t i_vert = closest_vertex( pt, result ) ;

   // Faster algorithm
   if( !result )
   {
      size_t i_cell = cell_index( pt, i_vert, POLY, POLY_VERTS ) ;
      result = ( i_cell!=PEL::bad_index() ) ;
   }
   // In failure case, global search in grid
   if( !result )
   {
      i_vert = closest_vertex_global_search( pt, result ) ;
      if( !result )
      {
         size_t i_cell = cell_index( pt, i_vert, POLY, POLY_VERTS ) ;
         result = ( i_cell!=PEL::bad_index() ) ;
      }
   }

   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_ResultReader:: has_integration_domain( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: has_integration_domain" ) ;
   bool result = ( READER->integration_domain()!=0 ) ;
   PEL_CHECK_POST( IMPLIES( nb_space_dimensions()!=2, !result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
PDE_ResultReader:: inner_boundary( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: inner_boundary" ) ;
   PEL_CHECK_PRE( has_integration_domain() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   doubleArray2D const& result =
      READER->integration_domain()->doubleArray2D_data( "inner_boundary" ) ;
   
   PEL_CHECK_POST( result.index_bound(0)==nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_ResultReader:: has_field( std::string const& field_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: has_field" ) ;
   PEL_CHECK_PRE( !field_name.empty() ) ;
   PEL_CHECK_INV( invariant() ) ;
   bool result =  FIELDS.has( field_name ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
doubleVector const&
PDE_ResultReader::  field_value( std::string const& field_name,
                                 GE_Point const* pt ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader::  field_value" ) ;
   PEL_CHECK_PRE( !field_name.empty() ) ;
   PEL_CHECK_PRE( pt!=0 ) ;
   PEL_CHECK_PRE( pt->nb_coordinates()==nb_space_dimensions() ) ;
   PEL_CHECK_PRE( has_field( field_name ) ) ;
   PEL_CHECK_PRE( is_in_grid( pt ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   static doubleVector result(0) ;

   // Cell containing the current point :
   bool is_a_vertex = false ;
   size_t i_cell = PEL::bad_index() ;

   // Faster algorithm
   size_t i_vert = closest_vertex( pt, is_a_vertex ) ;
   if( !is_a_vertex )
   {
      i_cell = cell_index( pt, i_vert, POLY, POLY_VERTS ) ;
      // If faster algorithm fails, try a global search
      if( i_cell == PEL::bad_index() )
      {
	 i_vert = closest_vertex_global_search( pt, is_a_vertex ) ;
	 if( !is_a_vertex )
	 {
	    i_cell = cell_index( pt, i_vert, POLY, POLY_VERTS ) ;
	 }
	 else
	 {
	    i_cell = PEL::bad_index() ;
	 }
      }
   }

   // Value of the field :
   size_t i_field = FIELDS.index_of( field_name ) ;
   std::string stored_name = FNAMES( i_field ) ;
   size_t i_comp = COMP( i_field ) ;
   std::string const& location =
          READER->fields()->string_data( stored_name+"/location" ) ;
   doubleArray2D const& value_table =
      READER->fields()->doubleArray2D_data( stored_name+"/value" ) ;
   size_t nb_comps = value_table.index_bound(0) ;
   if( i_comp==PEL::bad_index() )
   {
      result.re_initialize( nb_comps ) ;
   }
   else
   {
      if( i_comp>=nb_comps )
      {
         PDE_ResultReader_ERROR::n0( 
            "Invalid number of components for the field \""+field_name+"\"" ) ;
      }
      result.re_initialize( 1 ) ;
   }
   result.set( -PEL::max_double() ) ;
   if( location=="at_vertices" )
   {
      if( value_table.index_bound(1)!=NB_VERTS )
      {
         PDE_ResultReader_ERROR::n0(
            "Invalid number of values for the field \""+field_name+"\"" ) ;
      }
      if( is_a_vertex )
      {
         if( i_comp==PEL::bad_index() )
         {
            for( size_t i=0 ; i<nb_comps ; ++i )
            {
               result(i) = value_table( i, i_vert ) ;
            }
         }
         else
         {
            result(0) = value_table( i_comp, i_vert ) ;
         }
      }
      else
      {
         set_value_at_vertices( pt, i_cell, POLY, i_comp, value_table,
                                result ) ;
      }
   }
   else if( location=="at_cell_centers" )
   {
      if( value_table.index_bound(1)!=NB_CELLS )
      {
         PDE_ResultReader_ERROR::n0(
            "Invalid number of values for the field \""+field_name+"\"" ) ;
      }
      if( is_a_vertex )
      {
         set_value_at_cell_centers( pt, i_vert, i_comp,
                                    POLY, POLY_VERTS,
                                    value_table, result ) ;
      }
      else
      {
         if( i_comp==PEL::bad_index() )
         {
            for( size_t i=0 ; i<nb_comps ; ++i )
            {
               result(i) = value_table( i, i_cell ) ;
            }
         }
         else
         {
            result(0) = value_table( i_comp, i_cell ) ;
         }
      }
   }

   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultReader:: closest_vertex_global_search( GE_Point const* pt,
                                                 bool& is_a_vertex ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader::closest_vertex_global_search" ) ;
   PEL_CHECK( pt!=0 && pt->nb_coordinates()==NB_SP_DIMS ) ;

   is_a_vertex = false ;
   size_t result = PEL::bad_index() ;
   double d_closest = PEL::max_double() ;

   for( size_t i=0 ; i<NB_VERTS && !is_a_vertex; ++i )
   {
       VERTS.extract_section( 1, i, VERT ) ;
       double d_min = distance( pt, VERT ) ;
       if( d_min<d_closest )
       {
          result = i ;
          d_closest = d_min ;
          if( d_min<EPS_VERT ) is_a_vertex = true ;
       }
   }

   if( DEBUG )
   {
      VERTS.extract_section( 1, result, VERT ) ;
      PEL::out() << " Closest vertex of: " ;
      pt->print( PEL::out(),0 ) ;
      PEL::out() << " is " << VERT << ", is_a_vertex: " << is_a_vertex << std::endl ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultReader:: closest_vertex( GE_Point const* pt,
                                   bool& is_a_vertex ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader::closest_vertex" ) ;
   PEL_CHECK( pt!=0 && pt->nb_coordinates()==NB_SP_DIMS ) ;

   VISITED_VERT.set( false ) ;

   CLOSEST_VERT = 0 ;
   VERTS.extract_section( 1, CLOSEST_VERT, VERT ) ;


   double d_min = distance( pt, VERT ) ;
   
   bool found = ( d_min<EPS_VERT ) ;
   while( !found )
   {
      size_t next = index_of_neighbour( pt, CLOSEST_VERT, d_min ) ;
      if( next!=PEL::bad_index() )
      {
         CLOSEST_VERT = next ;
      }
      else
      {
         found = true ;
      }
   }
   is_a_vertex = ( d_min<EPS_VERT ) ;
   
   return( CLOSEST_VERT ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultReader:: cell_index( GE_Point const* pt,
                               size_t i_closest_vertex,
                               GE_Mpolyhedron* cell_poly,
                               PEL_Vector* cell_poly_verts ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: cell_index" ) ;
   PEL_CHECK( pt!=0 && pt->nb_coordinates()==NB_SP_DIMS ) ;
   PEL_CHECK( i_closest_vertex<NB_VERTS ) ;
   PEL_CHECK( cell_poly!=0 && cell_poly->dimension()==NB_SP_DIMS ) ;
   PEL_CHECK( cell_poly_verts!=0 &&
              cell_poly_verts->count()==cell_poly_verts->index_limit() ) ;
   PEL_CHECK( cell_poly->nb_vertices()==cell_poly_verts->count() ) ;
   PEL_CHECK( FORALL( ( size_t i=0 ; i<cell_poly->nb_vertices() ; ++i ),
                      cell_poly_verts->has( cell_poly->vertex(i) ) ) ) ;

   if( DEBUG )
   {
      PEL::out() << " NB_CELL_PER_VERT( i_closest_vertex ) : "
                 << NB_CELL_PER_VERT( i_closest_vertex ) << std::endl ;
   }

   // Find if pt is in the meshes connected to `i_closest_vertex' :
   size_t result = PEL::bad_index() ;
   for( size_t i=0 ;
        result==PEL::bad_index() && i<NB_CELL_PER_VERT( i_closest_vertex ) ;
        ++i )
   {
      size_t i_cell = VERT2CELL( i_closest_vertex, i ) ;
      for( size_t j=0 ; j<cell_poly_verts->count() ; ++j )
      {
         size_t i_vert = CELL2VERT( j, i_cell ) ;
         GE_Point* v = static_cast<GE_Point*>( cell_poly_verts->at(j) ) ;
         VERTS.extract_section( 1, i_vert, VERT ) ;
         v->set_coordinates( VERT ) ;
      }
      cell_poly ->update() ;

      if( DEBUG ) cell_poly->print( PEL::out(), 2 ) ;

      if( cell_poly->contains( pt ) )
      {
         result = i_cell ;
      }
   }
   
   PEL_CHECK( IMPLIES( result!=PEL::bad_index(), result<NB_CELLS ) ) ;
   PEL_CHECK( IMPLIES( result!=PEL::bad_index(), cell_poly->contains( pt ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_ResultReader::  set_value_at_vertices(
                         GE_Point const* pt,
                         size_t i_cell, GE_Mpolyhedron const* cell_poly,
                         size_t i_comp,
                         doubleArray2D const& value_table,
                         doubleVector& val ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: set_value_at_vertices" ) ;
   PEL_CHECK( pt!=0 && pt->nb_coordinates()==NB_SP_DIMS ) ;
   PEL_CHECK( i_cell<NB_CELLS ) ;
   PEL_CHECK( IMPLIES( i_comp!=PEL::bad_index(),
                       val.size()==1 ) ) ;
   PEL_CHECK( IMPLIES( i_comp!=PEL::bad_index(),
                       i_comp<value_table.index_bound(0) ) ) ;
   PEL_CHECK( IMPLIES( i_comp==PEL::bad_index(),
                       value_table.index_bound(0)==val.size() ) ) ;
   PEL_CHECK( value_table.index_bound(1)==NB_VERTS ) ;
   PEL_CHECK( cell_poly!=0 && cell_poly->contains( pt ) ) ;
   PEL_CHECK( REF_ELM->nb_nodes()==cell_poly->nb_vertices() ) ;

   val.set( 0. ) ;
   cell_poly->apply_inverse_mapping( pt, PT_REF ) ;
   for( size_t i=0 ; i<REF_ELM->nb_nodes() ; ++i )
   {
      // Index of the vertex corresponding to the node in the cell :
      cell_poly->apply_mapping( REF_ELM->node_location(i), PT ) ;
      doubleVector const& node_coords = PT->coordinate_vector() ;
      size_t i_vert = PEL::bad_index() ;
      double d_min = PEL::max_double() ;
      for( size_t ii=0 ; ii<cell_poly->nb_vertices() ; ++ii )
      {
         double d = distance( cell_poly->vertex(ii), node_coords ) ;
         if( d<d_min )
         {
            d_min = d ;
            i_vert = CELL2VERT( ii, i_cell ) ;
         }
      }

      // Finite element value :
      double const phi = REF_ELM->N_local( i, PT_REF ) ;
      if( i_comp==PEL::bad_index() )
      {
         for( size_t j=0 ; j<val.size() ; ++j )
         {
            val(j) += value_table( j, i_vert )*phi ;
         }
      }
      else
      {
         val(0) += value_table( i_comp, i_vert )*phi ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_ResultReader::  set_value_at_cell_centers(
                               GE_Point const* pt,
                               size_t i_vert,
                               size_t i_comp,
                               GE_Mpolyhedron* dummy_cell_poly,
                               PEL_Vector* dummy_cell_poly_verts,
                               doubleArray2D const& value_table,
                               doubleVector& val ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: set_value_at_cell_centers" ) ;
   PEL_CHECK( pt!=0 && pt->nb_coordinates()==NB_SP_DIMS ) ;
   PEL_CHECK( i_vert<NB_VERTS ) ;
   PEL_CHECK( IMPLIES( i_comp!=PEL::bad_index(),
                       val.size()==1 ) ) ;
   PEL_CHECK( IMPLIES( i_comp!=PEL::bad_index(),
                       i_comp<value_table.index_bound(0) ) ) ;
   PEL_CHECK( IMPLIES( i_comp==PEL::bad_index(),
                       value_table.index_bound(0)==val.size() ) ) ;
   PEL_CHECK( value_table.index_bound(1)==NB_CELLS ) ;
   PEL_CHECK(
      dummy_cell_poly!=0 && dummy_cell_poly->dimension()==NB_SP_DIMS ) ;
   PEL_CHECK(
      dummy_cell_poly_verts!=0 &&
      dummy_cell_poly_verts->count()==dummy_cell_poly_verts->index_limit() ) ;
   PEL_CHECK( dummy_cell_poly->nb_vertices()==dummy_cell_poly_verts->count() ) ;
   PEL_CHECK(
      FORALL(
         ( size_t i=0 ; i<dummy_cell_poly->nb_vertices() ; ++i ),
         dummy_cell_poly_verts->has( dummy_cell_poly->vertex(i) ) ) ) ;
   
   val.set( 0. ) ;

   double w_tot = 0. ;
   for( size_t i=0 ; i<NB_CELL_PER_VERT( i_vert ) ; ++i )
   {
      size_t i_cell = VERT2CELL( i_vert, i ) ;
      for( size_t j=0 ; j<dummy_cell_poly_verts->count() ; ++j )
      {
         size_t i_v = CELL2VERT( j, i_cell ) ;
         GE_Point* v = static_cast<GE_Point*>( dummy_cell_poly_verts->at(j) ) ;
         VERTS.extract_section( 1, i_v, VERT ) ;
         v->set_coordinates( VERT ) ;
      }
      dummy_cell_poly ->update() ;
      double const w = dummy_cell_poly->measure() ;
      w_tot += w ;
      if( i_comp==PEL::bad_index() )
      {
         for( size_t j=0 ; j<val.size() ; ++j )
         {
            val(j) += w*value_table( j, i_cell ) ;
         }
      }
      else
      {
         val(0) += w*value_table( i_comp, i_cell ) ;
      }
   }
   for( size_t j=0 ; j<val.size() ; ++j )
   {
      val(j) /= w_tot ;
   }
}

//----------------------------------------------------------------------
double
PDE_ResultReader:: distance( GE_Point const* pt1,
                             doubleVector const& pt2 ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: distance" ) ;
   PEL_CHECK( pt1!=0 && pt1->nb_coordinates()==NB_SP_DIMS ) ;
   PEL_CHECK( pt2.size()==NB_SP_DIMS ) ;
   double result = 0. ;
   for( size_t i=0 ; i<NB_SP_DIMS ; ++i )
   {
      result += PEL::abs( pt1->coordinate(i)-pt2(i) ) ;
//      result += ( pt1->coordinate(i)-pt2(i) )*( pt1->coordinate(i)-pt2(i) ) ;
   }
//   return( PEL::sqrt(result) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ResultReader:: index_of_neighbour(
                       GE_Point const* pt,
                       size_t const current_vert,
                       double& d_min ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ResultReader:: index_of_neighbour" ) ;
   PEL_CHECK( pt!=0 && pt->nb_coordinates()==NB_SP_DIMS ) ;
   PEL_CHECK( current_vert<NB_VERTS ) ;
   size_t result = PEL::bad_index() ;
   size_t n = CELL2VERT.index_bound(0) ;
   for( size_t i=0 ; i<NB_CELL_PER_VERT(current_vert) ; ++i )
   {
      size_t const i_mesh = VERT2CELL( current_vert, i ) ;
      for( size_t j=0 ; j<n ; ++j )
      {
         size_t const i_vert = CELL2VERT(j,i_mesh) ;
         if( !VISITED_VERT( i_vert ) )
         {
            VISITED_VERT( i_vert ) = true ;
            VERTS.extract_section( 1, i_vert, VERT ) ;
            double const d = distance( pt, VERT ) ;
            if( d<d_min )
            {
               d_min = d ;
               result = i_vert ;
            }
         }
      }
   }
   PEL_CHECK_POST( IMPLIES( result!=PEL::bad_index(), result<NB_VERTS ) ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ResultReader_ERROR:: n0( std::string const& mesg )
//internal--------------------------------------------------------------
{
   std::ostringstream err ;
   err << "*** " << "PDE_ResultReader:: invalid input file" << std::endl ;
   err << "*** " << mesg << std::endl ;
   PEL_Error::object()->raise_plain( err.str() ) ;
}
