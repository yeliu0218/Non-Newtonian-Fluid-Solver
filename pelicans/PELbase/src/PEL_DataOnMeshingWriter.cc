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

#include <PEL_DataOnMeshingWriter.hh>

#include <PEL_DoubleArray2D.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_IntVector.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_String.hh>
#include <PEL_StringVector.hh>
#include <PEL_assertions.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <intArray2D.hh>
#include <intVector.hh>
#include <stringVector.hh>

#include <limits>
#include <iostream>
#include <sstream>

using std::string ;
using std::endl ;
using std::ostringstream ;

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingWriter:: is_parallel_writer( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
PEL_DataOnMeshingWriter*
PEL_DataOnMeshingWriter:: make( PEL_Object* a_owner,
                                std::string const& name,
                                PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataOnMeshingWriter:: make" ) ;
   PEL_CHECK_PRE( !name.empty() ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   PEL_DataOnMeshingWriter const* proto =
      static_cast<PEL_DataOnMeshingWriter const*>(
                                    plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
      
   PEL_DataOnMeshingWriter* result = proto->create_replica( a_owner, exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_DataOnMeshingWriter:: PEL_DataOnMeshingWriter( std::string const& name )
//----------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
{
   PEL_LABEL( "PEL_DataOnMeshingWriter:: PEL_DataOnMeshingWriter" ) ;

   plugins_map()->register_item( name, this ) ;
   
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
PEL_DataOnMeshingWriter:: PEL_DataOnMeshingWriter( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
{
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//----------------------------------------------------------------------
int
PEL_DataOnMeshingWriter:: index_for_trash( void )
//----------------------------------------------------------------------
{
   return( -1 ) ;
}

//----------------------------------------------------------------------
double
PEL_DataOnMeshingWriter:: undefined_value( void )
//----------------------------------------------------------------------
{
   return( PEL::bad_double() ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_DataOnMeshingWriter:: create_meshing_module( PEL_Object* a_owner,
          std::string const& meshing_name,
          int nb_sp_dims, doubleArray2D const& vertices,
          intVector const& cell_nb_vertices, intArray2D const& cell2vertex,
          intVector const& cell_nb_faces,    intArray2D const& cell2face,
          intVector const& face_nb_vertices, intArray2D const& face2vertex,
          intArray2D const& face2cell,
          stringVector const& color_table,
          std::string const& halo_color_name,
          intArray2D const& macro_colors,
          intVector const& vertex_color,
          intVector const& cell_color,
          intVector const& face_color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataOnMeshingWriter:: create_meshing_module" ) ;
   PEL_CHECK_PRE( !meshing_name.empty() ) ;
   PEL_CHECK_PRE( vertices.index_bound( 0 ) == (size_t)nb_sp_dims ) ;
   PEL_CHECK_PRE( vertices.index_bound( 1 ) != 0 ) ;
   PEL_CHECK_PRE( cell2vertex.index_bound( 0 ) != 0 ) ;
   PEL_CHECK_PRE( cell2vertex.index_bound( 1 ) != 0 ) ;
   PEL_CHECK_PRE( cell2vertex.index_bound( 1 ) == cell_nb_vertices.size() ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_nb_vertices.size() ; ++ic ),
         ( cell_nb_vertices( ic ) > 0 ) &&
         ( (size_t)cell_nb_vertices( ic ) <= cell2vertex.index_bound( 0 ) ))) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_nb_vertices.size() ; ++ic ),
         FORALL( ( size_t iv=0 ; iv<(size_t)cell_nb_vertices( ic ) ; ++iv ),
            cell2vertex( iv, ic ) >= 0  ) ) ) ;
   PEL_CHECK_PRE( cell2face.index_bound( 0 ) != 0 ) ;
   PEL_CHECK_PRE( cell2face.index_bound( 1 ) == cell2vertex.index_bound( 1 ) );
   PEL_CHECK_PRE( cell2face.index_bound( 1 ) == cell_nb_faces.size() ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_nb_faces.size() ; ++ic ),
         ( cell_nb_faces( ic ) > 0 ) &&
         ( (size_t)cell_nb_faces( ic ) <= cell2face.index_bound( 0 ) ) ) )  ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_nb_faces.size() ; ++ic ),
         FORALL( ( size_t f=0 ; f<(size_t)cell_nb_faces( ic ) ; ++f ),
            cell2face( f, ic ) >= 0  ) ) ) ;
   PEL_CHECK_PRE( face2vertex.index_bound( 0 ) != 0 ) ;
   PEL_CHECK_PRE( face2vertex.index_bound( 1 ) != 0 ) ;
   PEL_CHECK_PRE( face2vertex.index_bound( 1 ) == face_nb_vertices.size() ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_nb_vertices.size() ; ++f ),
         ( face_nb_vertices( f ) > 0 ) && 
         ( (size_t)face_nb_vertices( f ) <= face2vertex.index_bound( 0 ) ) )) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_nb_vertices.size() ; ++f ),
         FORALL( ( size_t iv=0 ; iv<(size_t)face_nb_vertices( f ) ; ++iv ),
            face2vertex( iv, f ) >= 0  ) ) ) ;
   PEL_CHECK_PRE( face2cell.index_bound( 0 ) == 2 ) ;
   PEL_CHECK_PRE( face2cell.index_bound( 1 ) == face_nb_vertices.size() ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_nb_vertices.size() ; ++f ),
         face2cell( 0, f ) >= 0  ) ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_nb_vertices.size() ; ++f ),
         ( face2cell( 1, f ) == index_for_trash() ) ||
         ( face2cell( 1, f ) >= 0 ) ) ) ;
   PEL_CHECK_PRE( color_table.size() != 0 ) ;
   PEL_CHECK_PRE(
      EXISTS( ( size_t i=0 ; i<color_table.size() ; ++i ),
              halo_color_name == color_table( i ) ) ) ;
   PEL_CHECK_PRE( vertex_color.size() == vertices.index_bound( 1 ) ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t iv=0 ; iv<vertex_color.size() ; ++iv ),
         ( vertex_color( iv ) >= 0 ) &&
         ( (size_t)vertex_color( iv ) < color_table.size() ) ) ) ;
   PEL_CHECK_PRE( cell_color.size() == cell_nb_vertices.size( ) ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t ic=0 ; ic<cell_color.size() ; ++ic ),
         ( cell_color( ic ) >= 0 ) &&
         ( (size_t)cell_color( ic ) < color_table.size() ) ) ) ;
   PEL_CHECK_PRE( face_color.size() == face_nb_vertices.size( ) ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t f=0 ; f<face_color.size() ; ++f ),
         ( face_color( f ) >= 0 ) &&
         ( (size_t)face_color( f ) < color_table.size() ) ) ) ;

   PEL_Module* result = PEL_Module::create( a_owner, "meshing" ) ;

   result->add_entry( "name", PEL_String::create( result, meshing_name ) ) ;
   result->add_entry( "nb_sp_dims", PEL_Int::create( result, nb_sp_dims ) ) ;
   result->add_entry( "vertices", 
                      PEL_DoubleArray2D::create( result, vertices ) ) ;
   result->add_entry( "type", PEL_String::create( result, "meshing" ) ) ;
   result->add_entry( "cell_nb_vertices", 
                      PEL_IntVector::create( result, cell_nb_vertices ) ) ;
   result->add_entry( "cell2vertex", 
                      PEL_IntArray2D::create( result, cell2vertex ) ) ;

   result->add_entry( "cell_nb_faces", 
                      PEL_IntVector::create( result, cell_nb_faces ) ) ;
   result->add_entry( "cell2face",
                      PEL_IntArray2D::create( result, cell2face ) ) ;
   result->add_entry( "face_nb_vertices", 
                      PEL_IntVector::create( result, face_nb_vertices ) ) ;
   result->add_entry( "face2vertex", 
                      PEL_IntArray2D::create( result, face2vertex ) ) ;
   result->add_entry( "face2cell",
                      PEL_IntArray2D::create( result, face2cell ) ) ;
   result->add_entry( "color_table", 
                      PEL_StringVector::create( result, color_table ) ) ;
   result->add_entry( "color_table_connectivity", 
                      PEL_IntArray2D::create( result, macro_colors ) ) ;   
   result->add_entry( "halo_color_name", 
                      PEL_String::create( result, halo_color_name ) ) ;
   result->add_entry( "vertex_color",
                      PEL_IntVector::create( result, vertex_color ) ) ;
   result->add_entry( "cell_color",
                      PEL_IntVector::create( result, cell_color ) ) ;
   result->add_entry( "face_color",
                      PEL_IntVector::create( result, face_color ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == "meshing" ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_DataOnMeshingWriter:: create_field_module(
                              PEL_Object* a_owner,
                              std::string const& field_name,
                              std::string const& meshing, 
                              std::string const& location,
                              doubleArray2D const& value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataOnMeshingWriter:: create_field_module" ) ;

   PEL_Module* result = PEL_Module::create( a_owner, field_name ) ;
   result->add_entry( "name", PEL_String::create( result, field_name ) ) ;
   result->add_entry( "value", PEL_DoubleArray2D::create( result, value ) ) ;
   result->add_entry( "type", PEL_String::create( result, "field" ) ) ;
   result->add_entry( "location", PEL_String::create( result, location ) ) ;
   result->add_entry( "meshing", PEL_String::create( result, meshing ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == field_name ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DataOnMeshingWriter:: ~PEL_DataOnMeshingWriter( void  )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingWriter:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
void
PEL_DataOnMeshingWriter:: raise_field_location_error(
                   std::string const& field_name,
                   std::string const& location,
                   stringVector const& allowed_locations ) const
//----------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** " << type_name() << " error:" << endl << endl ;
   mesg << "    unable to save field of name \""
        << field_name << "\"" << endl ;
   mesg << "    at location \"" << location << "\"." << endl << endl ;
   mesg << "    allowed locations are:" ;
   for( size_t i=0 ; i<allowed_locations.size() ; ++i )
   {
      mesg << endl
           << "       - \"" << allowed_locations(i) << "\"" ;
   }
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingWriter:: write_cycle_PRE(
                                        PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( exp != 0 ) ;
   PEL_ASSERT( exp->name().substr(0,6) == "cycle_" ) ;
   PEL_ASSERT( exp->has_entry( "cycle_number" ) ) ;
   PEL_ASSERT( IMPLIES( exp->has_module( "meshing" ),
                        exp->has_entry( "meshing/nb_sp_dims" ) &&
                        exp->has_entry( "meshing/vertices" ) &&
                        exp->has_entry( "meshing/cell2vertex" ) ) ) ;
   PEL_ASSERT( IMPLIES( exp->has_module( "integration_domain" ),
                        exp->has_entry( "integration_domain/inner_boundary" ) &&
                        exp->has_entry( "integration_domain/polygon" ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingWriter:: create_replica_PRE( 
                                        PEL_Object* a_owner,
                                        PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataOnMeshingWriter:: create_replica_POST( 
                                    PEL_DataOnMeshingWriter const* result,
                                    PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
PEL_DataOnMeshingWriter:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "PEL_DataOnMeshingWriter descendant" ) ;
   return( result ) ;
}
