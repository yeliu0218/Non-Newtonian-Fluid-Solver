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

#include <GE_EMC2Meshing.hh>

#include <PEL.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>

GE_EMC2Meshing const*
GE_EMC2Meshing::PROTOTYPE = new GE_EMC2Meshing() ;

//-------------------------------------------------------------------------
GE_EMC2Meshing:: GE_EMC2Meshing( void )
//------------------------------------------------------------------------
   : GE_Meshing( "GE_EMC2Meshing" )
   , INPUT_FILE()
   , INPUT_FORMAT()
   , VERT_POS()
   , CELL_POS()
   , NB_VERTS( PEL::bad_index() )
   , NB_FACES( PEL::bad_index() )
   , NB_TRIAS( PEL::bad_index() )
   , NB_QUADS( PEL::bad_index() )
   , FACE_POLY_NAMES( 0 )
   , CELL_POLY_NAMES( 0 )
   , VERT_IDX( 0 )
   , CELLS_2_FACES( 0, 0 )
   , FACES_2_VERTS( 0, 0 )
   , I_VERT( PEL::bad_index() )
   , I_FACE( PEL::bad_index() )
   , I_CELL( PEL::bad_index() )
   , DEF_COLOR( 0 )
   , VERT_COORDS( 0 )
   , MESH_2_VERTS( 0 )
   , CELL_2_FACES( 0 )
{
}

//-------------------------------------------------------------------------
GE_EMC2Meshing*
GE_EMC2Meshing:: create_replica( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp,
                                 size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;

   if( dim_space != 2 ) raise_invalid_nb_space_dimensions( dim_space, "2" ) ;
   
   GE_EMC2Meshing* result = new GE_EMC2Meshing( a_owner, exp ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return result ;
}

//-------------------------------------------------------------------------
GE_EMC2Meshing:: GE_EMC2Meshing( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp )
//------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, 2 )
   , INPUT_FILE(
      exp->string_data( "filename" ).c_str(), std::ios_base::binary )
   , INPUT_FORMAT(
      ( exp->string_data( "format" ) == "ftq" ) ? ftq : amdba )
   , VERT_POS()
   , CELL_POS()
   , NB_VERTS( PEL::bad_index() )
   , NB_FACES( PEL::bad_index() )
   , NB_TRIAS( PEL::bad_index() )
   , NB_QUADS( PEL::bad_index() )
   , FACE_POLY_NAMES( 0 )
   , CELL_POLY_NAMES( 0 )
   , VERT_IDX( 0 )
   , CELLS_2_FACES( 0, 4 )
   , FACES_2_VERTS( 0, 2 )
   , I_VERT( PEL::bad_index() )
   , I_FACE( PEL::bad_index() )
   , I_CELL( PEL::bad_index() )
   , DEF_COLOR( 0 )
   , VERT_COORDS( 2 )
   , MESH_2_VERTS( 0 )
   , CELL_2_FACES( 0 )
{
   PEL_LABEL( "GE_EMC2Meshing:: GE_EMC2Meshing" ) ;

   // Input file check:
   exp->test_file( "filename", "read" ) ;
   if( !INPUT_FILE )
   {
      PEL_Error::object()->raise_file_handling(
                   exp->string_data( "filename" ), "read" ) ;
   }
   exp->test_data_in( "format", "amdba,ftq" ) ;

   // Rounding:
   initialize_rounding_strategy( exp ) ;
   
   // Meshing sizes:
   if( INPUT_FORMAT==ftq )
   {
      size_t nb_tot = PEL::bad_index() ;
      INPUT_FILE >> NB_VERTS >> nb_tot >> NB_TRIAS >> NB_QUADS ;
      PEL_ASSERT( nb_tot == NB_TRIAS+NB_QUADS ) ;
   }
   else
   {
      NB_QUADS = 0 ;
      INPUT_FILE >> NB_VERTS >> NB_TRIAS ;
   }
   if( NB_TRIAS+NB_QUADS==0 )
   {
      PEL_Error::object()->raise_plain(
         "Bad or empty EMC2 file " +exp->string_data( "filename" ) ) ;
   }

   // Meshing polyhedra:
   FACE_POLY_NAMES.resize(3) ;
   if( exp->has_entry( "face_polyhedron" ) )
   {
      stringVector const& poly = exp->stringVector_data( "face_polyhedron" ) ;
      for( size_t i=0 ; i<poly.size() ; ++i )
      {
         GE_ReferencePolyhedron const* ref =
            GE_Mpolyhedron::reference_polyhedron( poly(i) ) ;
         if( ref->dimension() != nb_space_dimensions()-1 )
         {
            PEL_Error::object()->raise_data_error(
               exp, "face_polyhedron",
               "bad dimension for face polyhedron \""+poly(i)+"\"" ) ;
         }
         size_t const nb_verts = ref->nb_vertices() ;
         if( FACE_POLY_NAMES.size()<=nb_verts )
         {
            FACE_POLY_NAMES.resize( nb_verts+1 ) ;
         }
         if( !FACE_POLY_NAMES(nb_verts).empty() )
         {
            PEL_Error::object()->raise_data_error(
               exp, "face_polyhedron",
               "\""+FACE_POLY_NAMES(nb_verts)+"\" and \""+poly(i)+"\""
               " defined the same polyhedron type" ) ;
         }
         FACE_POLY_NAMES(nb_verts) = poly(i) ;
      }
   }
   if( FACE_POLY_NAMES(2).empty() ) FACE_POLY_NAMES(2) = "GE_Segment" ;
   
   CELL_POLY_NAMES.resize(5) ;
   if( exp->has_entry( "cell_polyhedron" ) )
   {
      stringVector const& poly = exp->stringVector_data( "cell_polyhedron" ) ;
      for( size_t i=0 ; i<poly.size() ; ++i )
      {
         GE_ReferencePolyhedron const* ref =
            GE_Mpolyhedron::reference_polyhedron( poly(i) ) ;
         if( ref->dimension() != nb_space_dimensions() )
         {
            PEL_Error::object()->raise_data_error(
               exp, "cell_polyhedron",
               "bad dimension for cell polyhedron \""+poly(i)+"\"" ) ;
         }
         size_t const nb_verts = ref->nb_vertices() ;
         if( CELL_POLY_NAMES.size()<=nb_verts )
         {
            CELL_POLY_NAMES.resize( nb_verts+1 ) ;
         }
         if( !CELL_POLY_NAMES(nb_verts).empty() )
         {
            PEL_Error::object()->raise_data_error(
               exp, "cell_polyhedron",
               "\""+CELL_POLY_NAMES(nb_verts)+"\" and \""+poly(i)+"\""
               " defined the same polyhedron type" ) ;
         }
         CELL_POLY_NAMES(nb_verts) = poly(i) ;
      }
   }
   if( CELL_POLY_NAMES(3).empty() ) CELL_POLY_NAMES(3) = "GE_Triangle" ;

   // NB : we assume vertex are read before meshes
   VERT_IDX.re_initialize( NB_VERTS, 0 ) ;   
   
   std::string tmp ;
   getline( INPUT_FILE, tmp ) ;
   if( INPUT_FORMAT == amdba )
   {
      VERT_POS = INPUT_FILE.tellg() ;
      for( size_t iLine=0 ; iLine<NB_VERTS ; iLine++ )
      {
         getline( INPUT_FILE, tmp ) ;
      }
   }
   CELL_POS = INPUT_FILE.tellg() ;

   // Building face connectivity
   NB_FACES = 0 ;
   CELLS_2_FACES.raise_first_index_bound( NB_TRIAS+NB_QUADS ) ;
   FACES_2_VERTS.raise_first_index_bound( 3*NB_TRIAS+4*NB_QUADS ) ;
   PEL_BalancedBinaryTree* side_tree =
                               PEL_BalancedBinaryTree:: create( 0 ) ;
   size_t_vector s_verts(2) ;
   for( size_t i=0 ; i<NB_TRIAS+NB_QUADS ; i++ )
   {
      PEL_ASSERT( !INPUT_FILE.eof() ) ;
      size_t idx, nb_cell_verts, i_col, s0, s1, s2 ;
      size_t_vector inod(4) ;
      
      if( INPUT_FORMAT==ftq )
      {
         INPUT_FILE >> nb_cell_verts ;
         if( nb_cell_verts != 3 && nb_cell_verts != 4 )
         {
            PEL_Error::object()->raise_plain(
               "*** GE_EMC2Meshing error:\n"
               "    reading file \""
                       +exp->string_data( "filename" )+"\":\n"
               "    only triangular or quadrangular cells are expected" ) ;
         }
         if( nb_cell_verts == 4 &&
             ( CELL_POLY_NAMES.size()<5 || CELL_POLY_NAMES(4).empty() ) )
         {
            PEL_Error::object()->raise_module_error(
               exp,
               "entry \"cell_polyhedron\" is mandatory with quadrangular cells\n" ) ;
         }
      }
      else 
      {
         INPUT_FILE >> idx ;
         nb_cell_verts = 3 ;
      }
      INPUT_FILE >> s0 ;
      
      // Fortran numbering
      s0-- ;
      s1 = s0 ;
      for( size_t j=0 ; j<nb_cell_verts ; ++j )
      {
         if( j==nb_cell_verts-1 )
         {
            s2 = s0 ;
         }
         else
         {
            INPUT_FILE >> s2 ;
            // Fortran numbering
            s2-- ;
         }
         inod(j)=s2;
         s_verts(0) = s1 ;
         s_verts(1) = s2 ;
         PEL_IndexSet* f = PEL_IndexSet::create( 0, s_verts, NB_FACES ) ;
         if( ! side_tree->has( f ) )
         {
            FACES_2_VERTS(NB_FACES,0) = ( s1<s2 ? s1 : s2 ) ;
            FACES_2_VERTS(NB_FACES,1) = ( s1<s2 ? s2 : s1 ) ;
            CELLS_2_FACES(i,j) = NB_FACES ;
            ++NB_FACES ;
            side_tree->extend( f ) ;
            f->set_owner( side_tree ) ;
         }
         else
         {
            PEL_IndexSet const* l =
               static_cast<PEL_IndexSet const*>( side_tree->item( f ) ) ;
            CELLS_2_FACES(i,j) = l->id() ;
            f->destroy() ; f = 0 ;
         }
         s1 = s2 ;
      }
      INPUT_FILE >> i_col ;
      if( i_col!=0 )
      {
         for( size_t n=0 ; n<nb_cell_verts ; ++n )
         {
            VERT_IDX(inod(n)) = i_col ;
         }
      }
   }
   side_tree->destroy() ; side_tree = 0 ;
   
   if( INPUT_FORMAT == ftq )
   {
      getline( INPUT_FILE, tmp ) ;
      VERT_POS = INPUT_FILE.tellg() ;
   }
   
   I_VERT = NB_VERTS+1 ;
   I_FACE = NB_FACES+1 ;
   I_CELL = NB_TRIAS+NB_QUADS+1 ;

}

//------------------------------------------------------------------------
GE_EMC2Meshing:: ~GE_EMC2Meshing( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
size_t
GE_EMC2Meshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   return( NB_VERTS ) ;
}

//------------------------------------------------------------------------
size_t
GE_EMC2Meshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   return( NB_TRIAS+NB_QUADS ) ;
}

//------------------------------------------------------------------------
size_t
GE_EMC2Meshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   return( NB_FACES ) ;
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;
   INPUT_FILE.seekg( VERT_POS ) ;
   I_VERT = 0 ;
   read_vertex() ;
}

//------------------------------------------------------------------------
bool
GE_EMC2Meshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( I_VERT<NB_VERTS ) ;   
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;
   I_VERT++ ;
   if( valid_vertex() ) read_vertex() ;
}

//------------------------------------------------------------------------
doubleVector const& 
GE_EMC2Meshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;
   doubleVector const& result = VERT_COORDS ;
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   I_CELL = 0 ;
   INPUT_FILE.seekg( CELL_POS ) ;

   read_cell() ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_EMC2Meshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: valid_cell" ) ;
   bool result = ( I_CELL<NB_TRIAS+NB_QUADS )  ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;

   I_CELL++ ;
   if( valid_cell() )
   {
      read_cell() ;
   }
}

//------------------------------------------------------------------------
std::string const&
GE_EMC2Meshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;
   
   size_t const nb_verts = MESH_2_VERTS.size() ;
   std::string const& result = CELL_POLY_NAMES( nb_verts ) ;
  
   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_EMC2Meshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = MESH_2_VERTS ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_EMC2Meshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = CELL_2_FACES ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;

   I_FACE = 0 ;
   MESH_2_VERTS.resize( 2 ) ;

   read_face() ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_EMC2Meshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: valid_face" ) ;

   bool result = ( I_FACE< NB_FACES ) ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   I_FACE++ ;
   if( valid_face() )
   {
      read_face() ;
   }
}

//------------------------------------------------------------------------
std::string const&
GE_EMC2Meshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   size_t const nb_verts = MESH_2_VERTS.size() ;
   std::string const& result = FACE_POLY_NAMES( nb_verts ) ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_EMC2Meshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = MESH_2_VERTS ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_EMC2Meshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string str = "r" ;
   str += ( '0' + VERT_IDX( I_VERT ) ) ;
   GE_Color::extend( str ) ;
   GE_Color const* result = GE_Color::object( str ) ;

   PEL_CHECK_POST( default_vertex_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_EMC2Meshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = DEF_COLOR ;

   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_EMC2Meshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = DEF_COLOR ;

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Meshing:: print( os, indent_width ) ;
   
   std::string const s( indent_width+3, ' ' ) ;
   if( NB_TRIAS>0 )
   {
      os << s << NB_TRIAS << " \"" << CELL_POLY_NAMES(3) << "\"" << std::endl ;
   }
   if( NB_QUADS>0 )
   {
      os << s << NB_QUADS << " \"" << CELL_POLY_NAMES(4) << "\"" << std::endl ;
   }
   os << s << NB_FACES << " \"" << FACE_POLY_NAMES(2) << "\"" << std::endl ;
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: read_vertex( void )
//------------------------------------------------------------------------
{
   PEL_CHECK( valid_vertex() ) ;   
   PEL_CHECK_INV( invariant() ) ;
   PEL_ASSERT( !INPUT_FILE.eof() ) ;
   double x0, x1 ;
   size_t col ;
   if( INPUT_FORMAT==ftq )
   {
      INPUT_FILE >> x0 >> x1  >> col ;
   }
   else
   {
      size_t idx ;
      INPUT_FILE >> idx >> x0 >> x1  >> col ;
   }
   
   VERT_COORDS(0) = roundoff( x0 ) ;
   VERT_COORDS(1) = roundoff( x1 ) ;
   if( col!=0 ) VERT_IDX(I_VERT) = col ;
   
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: read_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: read_cell" ) ;
   PEL_CHECK( valid_cell() ) ;
   
   size_t nb_cell_verts = 3 ;
   if( INPUT_FORMAT == ftq )
   {
      INPUT_FILE >> nb_cell_verts ;
   }
   else
   {
      size_t idx ;
      INPUT_FILE >> idx ;
   }
   PEL_ASSERT( nb_cell_verts == 3 || nb_cell_verts == 4 ) ;
   
   MESH_2_VERTS.resize( nb_cell_verts ) ;
   for( size_t i=0 ; i<nb_cell_verts ; i++ )
   {
      INPUT_FILE >> MESH_2_VERTS(i) ;
      // Fortran numbering
      MESH_2_VERTS(i)-- ;
   }
   size_t col ;
   INPUT_FILE >> col ;
   std::string str_col = "r" ;
   str_col += (char)( '0' + col ) ;
   GE_Color::extend( str_col ) ;
   DEF_COLOR = GE_Color::object( str_col ) ;

   PEL_ASSERT( nb_cell_verts<CELL_POLY_NAMES.size() ) ;
   std::string const& ref_name = CELL_POLY_NAMES( nb_cell_verts ) ;
   PEL_ASSERT( !ref_name.empty() ) ;
   GE_ReferencePolyhedron const* the_mesh_ref_polyhedron =
      GE_Mpolyhedron::reference_polyhedron( ref_name ) ;
   size_t nb_conn = the_mesh_ref_polyhedron->nb_faces() ;
   CELL_2_FACES.resize( nb_conn ) ;
   for( size_t i=0 ; i<nb_cell_verts ; i++ )
   {
      CELL_2_FACES(i) = CELLS_2_FACES(I_CELL,i) ;
   }
}

//------------------------------------------------------------------------
void
GE_EMC2Meshing:: read_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_EMC2Meshing:: read_face" ) ;
   PEL_CHECK( valid_face() ) ;
   
   MESH_2_VERTS(0) = FACES_2_VERTS(I_FACE,0) ;
   MESH_2_VERTS(1) = FACES_2_VERTS(I_FACE,1) ;
   // Find color
   size_t_vector indic_vector( MESH_2_VERTS.size() ) ;
   size_t n = 0 ;
   for( size_t i=0 ; i<MESH_2_VERTS.size() ; i++ )
   {
      size_t indic = VERT_IDX( MESH_2_VERTS(i) ) ;
      
      bool found = false ;
      for( size_t j=0 ; !found && j<n ; j++ )
      {
         if( indic_vector(j)>=indic )
         {
            found = true ;
            if( indic_vector(j)>indic )
            {
               for( size_t k=n-1 ; k>=j && k<n ; k-- )
               {
                  indic_vector(k+1) = indic_vector(k) ;
               }
               n++ ;
               indic_vector(j) = indic ;
            }
         }
      }
      if( !found )
      {
         indic_vector(n) = indic ;
         n++ ;
      }
   }
   std::string col ;
   for( size_t i=0 ; i<n ; i++ )
   {
      col += "r" ;
      col += ('0'+indic_vector(i)) ;
   }
   GE_Color::extend( col ) ;
   DEF_COLOR = GE_Color::object( col ) ;
}

//------------------------------------------------------------------------
bool
GE_EMC2Meshing:: invariant( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Meshing::invariant() ) ;
   return( true ) ;
}


   
