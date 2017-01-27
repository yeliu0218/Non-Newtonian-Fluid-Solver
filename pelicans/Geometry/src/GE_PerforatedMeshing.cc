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

#include <GE_PerforatedMeshing.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Int.hh>
#include <PEL_IntVector.hh>
#include <PEL_String.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>

#include <doubleArray2D.hh>
//#include <intVector.hh>
//#include <stringVector.hh>
//#include <size_t_vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>

#include <iostream>
#include <sstream>

GE_PerforatedMeshing const*
GE_PerforatedMeshing:: PROTOTYPE = new GE_PerforatedMeshing() ;

//-------------------------------------------------------------------------
GE_PerforatedMeshing:: GE_PerforatedMeshing( void )
//------------------------------------------------------------------------
   : GE_Meshing( "GE_PerforatedMeshing" )
   , INITIAL(0 )
   , HOLE_COLORS( 0 )
   , NB_CELLS( PEL::bad_index() )
   , REMOVED_CELLS( 0 )
   , NB_VERTS( PEL::bad_index() )
   , INITIAL_TO_LOCAL_VERTS( 0 )
   , VERTS_COL( 0 )
   , NB_FACES( PEL::bad_index() )
   , INITIAL_TO_LOCAL_FACES( 0 )
   , FACES_COL( 0 )
   , IDX_VERT( PEL::bad_index() )
   , IDX_FACE( PEL::bad_index() )
   , IDX_CELL( PEL::bad_index() )
   , LOCAL_VERTS(0)
   , LOCAL_CONNS(0)
{
}

//-------------------------------------------------------------------------
GE_PerforatedMeshing*
GE_PerforatedMeshing:: create_replica( PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp,
                                  size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;

   GE_PerforatedMeshing* result =
      new GE_PerforatedMeshing(  a_owner, exp, dim_space ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, dim_space ) ) ;
   return result ;
}

//-------------------------------------------------------------------------
GE_PerforatedMeshing:: GE_PerforatedMeshing( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp,
                                   size_t dim_space )
//------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , INITIAL(0 )
   , HOLE_COLORS( 0 )
   , NB_CELLS( PEL::bad_index() )
   , REMOVED_CELLS( 0 )
   , NB_VERTS( PEL::bad_index() )
   , INITIAL_TO_LOCAL_VERTS( 0 )
   , VERTS_COL( 0 )
   , NB_FACES( PEL::bad_index() )
   , INITIAL_TO_LOCAL_FACES( 0 )
   , FACES_COL( 0 )
   , IDX_VERT( PEL::bad_index() )
   , IDX_FACE( PEL::bad_index() )
   , IDX_CELL( PEL::bad_index() )
   , LOCAL_VERTS(0)
   , LOCAL_CONNS(0)
{
   PEL_LABEL( "GE_PerforatedMeshing:: GE_PerforatedMeshing" ) ;

   // Read initial meshing:
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "GE_Meshing" ) ;
      INITIAL = GE_Meshing::create( this, se, dim_space ) ;
      se->destroy() ; se = 0 ;
   }
   
   // Formula context:
   PEL_ContextSimple* ctx = PEL_ContextSimple::create( 0 ) ;
   PEL_DoubleVector* coords =
         PEL_DoubleVector::create( ctx, nb_space_dimensions() ) ;
   ctx->extend( PEL_Variable::object( "DV_X" ), coords ) ;
   std::vector<PEL_Data*> hole_formula( 0 ) ;      
   
   // Explore holes:
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "holes" ) ;
      for( se->start_entry_iterator() ;
           se->is_valid_entry() ;
           se->go_next_entry() )
      {
         GE_Color::extend( se->keyword() ) ;
         HOLE_COLORS.push_back( GE_Color::object( se->keyword() ) ) ;
         PEL_Data* formula = se->data( ctx, ctx ) ;
         hole_formula.push_back( formula ) ;
         if( formula->data_type() != PEL_Data::Bool )
         {
            PEL_Error::object()->raise_bad_data_type(
                               se, se->keyword(), PEL_Data::Bool ) ;
         }
         if( !formula->value_can_be_evaluated( 0 ) )
         {
            PEL_Error::object()->raise_not_evaluable(
               se, se->keyword(), formula->undefined_variables( 0 ) ) ;
         }
      }
      se->destroy() ; se = 0 ;
   }

   // Set new meshing:
   build_meshing( hole_formula, coords ) ;

   ctx->destroy() ; ctx = 0 ; coords = 0 ;
}

//------------------------------------------------------------------------
GE_PerforatedMeshing:: ~GE_PerforatedMeshing( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
size_t
GE_PerforatedMeshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: nb_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;
  
   return( NB_VERTS ) ;
}

//------------------------------------------------------------------------
size_t
GE_PerforatedMeshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: nb_cells" ) ;
   return( NB_CELLS ) ;
}

//------------------------------------------------------------------------
size_t
GE_PerforatedMeshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: nb_faces" ) ;
   return( NB_FACES ) ;
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   IDX_VERT = 0 ;
   INITIAL->start_vertex_iterator() ;
   find_next_vertex() ;
}

//------------------------------------------------------------------------
bool
GE_PerforatedMeshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( INITIAL->valid_vertex() ) ;   
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   IDX_VERT++ ;   
   INITIAL->go_next_vertex() ;
   find_next_vertex() ;
}

//------------------------------------------------------------------------
doubleVector const& 
GE_PerforatedMeshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;
   
   doubleVector const& result = INITIAL->vertex_coordinates() ;
   
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   IDX_CELL = 0 ;
   INITIAL->start_cell_iterator() ;
   find_next_cell() ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;
   
   IDX_CELL++ ;
   INITIAL->go_next_cell() ;
   find_next_cell() ;
}

//------------------------------------------------------------------------
bool
GE_PerforatedMeshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: valid_cell" ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = INITIAL->valid_cell() ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
GE_PerforatedMeshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = INITIAL->cell_polyhedron_name() ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_PerforatedMeshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = LOCAL_VERTS ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_PerforatedMeshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = LOCAL_CONNS ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;
   
   IDX_FACE = 0 ;
   INITIAL->start_face_iterator() ;
   find_next_face() ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   IDX_FACE++ ;
   INITIAL->go_next_face() ;
   find_next_face() ;
}

//------------------------------------------------------------------------
bool
GE_PerforatedMeshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: valid_face" ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = INITIAL->valid_face() ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
GE_PerforatedMeshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = INITIAL->face_polyhedron_name() ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_PerforatedMeshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = LOCAL_VERTS ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_PerforatedMeshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: default_vertex_color" ) ;
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const i_col = VERTS_COL( IDX_VERT ) ;  
   GE_Color const* result =
      ( i_col != PEL::bad_index() ? HOLE_COLORS[i_col] :
                                  INITIAL->vertex_color() ) ;
   
   PEL_CHECK( default_vertex_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_PerforatedMeshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: default_cell_color" ) ;
   PEL_CHECK( default_cell_color_PRE() ) ;
   
   GE_Color const* result = INITIAL->cell_color() ;
   
   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_PerforatedMeshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: default_face_color" ) ;
   PEL_CHECK( default_face_color_PRE() ) ;

   size_t const i_col = FACES_COL( IDX_FACE ) ;
   GE_Color const* result =
      ( i_col != PEL::bad_index() ? HOLE_COLORS[i_col] :
                                  INITIAL->face_color() ) ;
   
   PEL_CHECK( default_face_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Meshing::print( os, indent_width ) ;
   INITIAL->print( os, indent_width+3 ) ;
   std::string const s = std::string( indent_width+3, ' ' ) ;
   os << s << "holes: " ;
   for( size_t i=0 ; i<HOLE_COLORS.size() ; ++i )
   {
      os << "\"" << HOLE_COLORS[i]->name() << "\" " ;
   }
   os << std::endl ;
   size_t nb_cs = 0 ;
   for( size_t i=0 ; i<REMOVED_CELLS.size() ; ++i )
   {
      if( REMOVED_CELLS(i) ) ++nb_cs ;
   }
   os << s << "nb removed cells: " << nb_cs << std::endl ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: display( std::ostream& os, size_t indent_width )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: display" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   os << "INITIAL mesh : " << std::endl ;
   INITIAL->display( os, indent_width+2 ) ;
   os << std::endl ;
   os << "Perforated mesh :   " << std::endl ;
   GE_Meshing::display( os, indent_width+2  ) ;
   os << std::endl ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: find_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: find_next_vertex" ) ;
   for(  ; INITIAL->valid_vertex() ; INITIAL->go_next_vertex() )
   {
      if( INITIAL_TO_LOCAL_VERTS( IDX_VERT ) != PEL::bad_index() )
         break ;
      IDX_VERT++ ;
   }
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: find_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: find_next_cell" ) ;
   
   for(  ; INITIAL->valid_cell() ; INITIAL->go_next_cell() )
   {
      if( ! REMOVED_CELLS( IDX_CELL ) )
         break ;
      IDX_CELL++ ;
   }
   if( INITIAL->valid_cell() )
   {
      size_t_vector const& vert = INITIAL->cell_vertices() ;
      size_t_vector const& conn = INITIAL->cell_faces() ;
      LOCAL_VERTS.re_initialize( vert.size() ) ;
      for( size_t i=0 ; i<vert.size()  ; i++ )
         LOCAL_VERTS(i) = INITIAL_TO_LOCAL_VERTS( vert(i) ) ;
      LOCAL_CONNS.re_initialize( conn.size() ) ;
      for( size_t i=0 ; i<conn.size()  ; i++ )
         LOCAL_CONNS(i) = INITIAL_TO_LOCAL_FACES( conn(i) ) ;
   }   
}

//------------------------------------------------------------------------
void
GE_PerforatedMeshing:: find_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: find_next_face" ) ;
   
   for(  ; INITIAL->valid_face() ; INITIAL->go_next_face() )
   {
      if( INITIAL_TO_LOCAL_FACES( IDX_FACE ) != PEL::bad_index() )
         break ;
      IDX_FACE++ ;
   }
   if( INITIAL->valid_face() )
   {
      size_t_vector const& vert = INITIAL->face_vertices() ;

      LOCAL_VERTS.re_initialize( vert.size() ) ;
      for( size_t i=0 ; i<vert.size()  ; i++ )
      {
         LOCAL_VERTS(i) = INITIAL_TO_LOCAL_VERTS( vert(i) ) ;
      }
   }   
}

//----------------------------------------------------------------------
void
GE_PerforatedMeshing:: build_meshing(
   std::vector<PEL_Data*> const& hole_formula, PEL_DoubleVector* coords )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerforatedMeshing:: build_meshing" ) ;
   PEL_CHECK( hole_formula.size() == HOLE_COLORS.size() ) ;
   PEL_CHECK( coords != 0 ) ;

   REMOVED_CELLS.re_initialize( INITIAL->nb_cells() ) ;
   VERTS_COL.re_initialize( INITIAL->nb_vertices() ) ;
   VERTS_COL.set( PEL::bad_index() ) ;
   FACES_COL.re_initialize( INITIAL->nb_faces() ) ;
   FACES_COL.set( PEL::bad_index() ) ;

   size_t const nb_sp_dims = nb_space_dimensions() ;

   // Build vertices:
   doubleArray2D vert_coords( INITIAL->nb_vertices(), nb_sp_dims ) ;
   INITIAL->start_vertex_iterator() ;
   for( size_t i_vert = 0 ;
        INITIAL->valid_vertex() ;
        INITIAL->go_next_vertex(), ++i_vert )
   {
      doubleVector const& xx = INITIAL->vertex_coordinates() ;
      for( size_t i=0 ; i<nb_sp_dims ; ++i )
         vert_coords( i_vert, i ) = xx( i ) ;
   }

   doubleVector center( nb_sp_dims ) ;
   size_t_vector v( INITIAL->nb_vertices() ) ;
   v.set( 0 ) ;
   size_t_vector f( INITIAL->nb_faces() ) ;
   f.set( 0 ) ;

   NB_CELLS = 0 ;
   INITIAL->start_cell_iterator() ;
   for( size_t idx = 0 ;
        INITIAL->valid_cell() ;
        INITIAL->go_next_cell(), ++idx )
   {
      // Compute cell center coordinates:
      center.set(0.) ;
      size_t_vector const& verts = INITIAL->cell_vertices() ;
      double const alpha = 1./verts.size() ;
      for( size_t i=0 ; i<verts.size() ; ++i )
      {
         for( size_t j=0 ; j<nb_sp_dims ; ++j )
         {
            center(j) += alpha*vert_coords( verts(i), j ) ;
         }
      }
      coords->set( center ) ;

      // Loop on holes:
      GE_Color const* cell_col = 0 ;
      size_t i_hole = PEL::bad_index() ;
      for( size_t i=0 ; i<HOLE_COLORS.size() ; ++i )
      {
         if( hole_formula[i]->to_bool() )
         {
            if( cell_col != 0 )
            {
               std::ostringstream mesg ;
               mesg << "GE_PerforatedMeshing : " << std::endl ;
               mesg << "   cell of center: " << coords << std::endl ;
               mesg << "   can't be removed as a hole of color "
                    << cell_col->name() << std::endl ;
               mesg << "   and as a hole of color " << HOLE_COLORS[i]->name() ;
               PEL_Error::object()->raise_plain( mesg.str() ) ;
            }
            cell_col = HOLE_COLORS[i] ;
            i_hole = i ;
            REMOVED_CELLS(idx) = true ;
         }
      }

      // Vertices and sides colors:
      size_t_vector const& faces = INITIAL->cell_faces() ;
      if( ! REMOVED_CELLS(idx) )
      {
         ++NB_CELLS ;
         for( size_t iv=0 ; iv<verts.size() ; ++iv )
         {
            v( verts(iv) ) += 1 ;
         }
         for( size_t is=0 ; is<faces.size() ; ++is )
         {
            f( faces(is) ) += 1 ;
         }
      }
      else
      {
         for( size_t iv=0 ; iv<verts.size() ; ++iv )
         {
            VERTS_COL( verts(iv) ) = i_hole ;
         }
         for( size_t is=0 ; is<faces.size() ; ++is )
         {
            FACES_COL( faces(is) ) = i_hole ;
         }
      }
   }

   // No cells: error...
   if( NB_CELLS==0 )
   {
      PEL_Error::object()->raise_plain(
         "Unconsistant perforated meshing : empty submeshing" ) ;
   }

   // Vertices:
   NB_VERTS = 0 ;
   INITIAL_TO_LOCAL_VERTS.re_initialize( INITIAL->nb_vertices() ) ;
   INITIAL_TO_LOCAL_VERTS.set( PEL::bad_index() ) ;
   for( size_t iv=0 ; iv<v.size() ; ++iv )
   {
      if( v(iv) > 0 ) INITIAL_TO_LOCAL_VERTS(iv) = NB_VERTS++ ;
   }

   // Faces:
   NB_FACES = 0 ;
   INITIAL_TO_LOCAL_FACES.re_initialize( INITIAL->nb_faces() ) ;
   INITIAL_TO_LOCAL_FACES.set( PEL::bad_index() ) ;
   for( size_t is=0 ; is<f.size() ; ++is )
   {
      if( f(is) > 0 ) INITIAL_TO_LOCAL_FACES(is) = NB_FACES++ ;
   }
}

