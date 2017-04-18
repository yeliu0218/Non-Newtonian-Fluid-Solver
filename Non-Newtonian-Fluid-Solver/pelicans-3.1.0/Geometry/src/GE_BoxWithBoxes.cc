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

#include <GE_BoxWithBoxes.hh>

#include <PEL_KeywordDataPair.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Data.hh>
#include <PEL_Error.hh>
#include <PEL_String.hh>
#include <PEL_Vector.hh>
#include <PEL.hh>
#include <intVector.hh>
#include <stringVector.hh>

#include <size_t_vector.hh>

#include <GE_Color.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Mpolyhedron.hh>

using std::string ;

//-------------------------------------------------------------------------
GE_BoxWithBoxes const* GE_BoxWithBoxes::PROTOTYPE = new GE_BoxWithBoxes() ;
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
GE_BoxWithBoxes:: GE_BoxWithBoxes( void )
//------------------------------------------------------------------------
  : GE_Meshing( "GE_BoxWithBoxes" )
  , FACE_POLY_NAME( "" )
  , FACE_POLY( 0 )
  , CELL_POLY_NAME( "" )
  , CELL_POLY( 0 )
  , IJK_COORDS( 0, 0 )
  , NB_VERTS( 0 )
  , NB_FACES( 0 )
  , IBD( 0 )
  , COLOR_TABLE( 0 )
  , GE_COLOR_TABLE( 0 )
  , I_POLY( 0 )
  , IBD_POLY( 0 )
  , POLY_TO_VERTS( 0 )
  , CELL_IT( false )
  , IVERT( PEL::bad_index() )
  , VCOORDS( 0 )
  , IFACE( PEL::bad_index() )
  , FACE_NORMAL( PEL::bad_index() )
  , CELL_TO_FACES( 0 )
{
}

//------------------------------------------------------------------------
GE_BoxWithBoxes:: ~GE_BoxWithBoxes( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
GE_BoxWithBoxes*
GE_BoxWithBoxes:: create_replica( PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp,
                                  size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;

   GE_BoxWithBoxes* result = new GE_BoxWithBoxes( a_owner, exp, dim_space ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dim_space ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
GE_BoxWithBoxes:: GE_BoxWithBoxes( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp,
                                   size_t dim_space )
//------------------------------------------------------------------------
  : GE_Meshing( a_owner, exp, dim_space )
  , FACE_POLY_NAME( "" )
  , FACE_POLY( 0 )
  , CELL_POLY_NAME( "" )
  , CELL_POLY( 0 )
  , IJK_COORDS( 0, dim_space )
  , NB_VERTS( dim_space )
  , NB_FACES( dim_space )
  , IBD( dim_space )
  , COLOR_TABLE( PEL_Vector::create( this, 2*dim_space ) )
  , GE_COLOR_TABLE( PEL_Vector::create( this, 2*dim_space ) )
  , I_POLY( dim_space )
  , IBD_POLY( dim_space )
  , POLY_TO_VERTS( 0 )
  , CELL_IT( false )
  , IVERT( PEL::bad_index() )
  , VCOORDS( 0 )
  , IFACE( PEL::bad_index() )
  , FACE_NORMAL( PEL::bad_index() )
  , CELL_TO_FACES( 0 )
{
   PEL_LABEL( "GE_BoxWithBoxes:: GE_BoxWithBoxes" ) ;
   
   read_mesh_polyhedron( exp, dim_space, FACE_POLY_NAME, CELL_POLY_NAME ) ;

   FACE_POLY = GE_Mpolyhedron::reference_polyhedron( FACE_POLY_NAME ) ;
   CELL_POLY = GE_Mpolyhedron::reference_polyhedron( CELL_POLY_NAME ) ;
   
   doubleVector const& left =
      exp->doubleVector_data( "vertices_coordinate_0" ) ;
   NB_VERTS(0) = left.size() ;
   size_t max_index = left.size() ;
   
   if( exp->has_entry( "vertices_coordinate_1" ) )
   {
      doubleVector const& bottom =
                     exp->doubleVector_data( "vertices_coordinate_1" ) ;
      max_index = PEL::max( max_index, bottom.size() ) ;
      NB_VERTS(1) = bottom.size() ;
   }
   
   if( exp->has_entry( "vertices_coordinate_2" ) )
   {
      if( dim_space != 3 )
      {
         raise_invalid_nb_space_dimensions(
            dim_space,
            "   3D is expexted together with\n"
            "   \"vertices_coordinate_2\" data" ) ;
      }
      doubleVector const& depth =
                     exp->doubleVector_data( "vertices_coordinate_2" ) ;
      max_index = PEL::max(  max_index, depth.size() ) ;
      NB_VERTS(2) = depth.size() ;
   }
   
   IJK_COORDS.raise_first_index_bound( max_index ) ;
   
   COLOR_TABLE->set_at( 0, PEL_String::create( this, "left" ) ) ;
   COLOR_TABLE->set_at( dim_space, PEL_String::create( this, "right" ) ) ;
   if( dim_space>1 )
   {
      COLOR_TABLE->set_at( 1, PEL_String::create( this, "bottom" ) ) ;
      COLOR_TABLE->set_at( dim_space+1, PEL_String::create( this, "top" ) ) ;
   }
   if( dim_space==3 )
   {
      COLOR_TABLE->set_at( 2, PEL_String::create( this, "behind" ) ) ;
      COLOR_TABLE->set_at( dim_space+2, PEL_String::create( this, "front" ) ) ;
   }
   
   IBD( 0 ) = left.size() ;
   for( size_t j = 0 ; j<left.size() ; j++ )
   {
      IJK_COORDS( j, 0 ) = left(j) ;
   }
   
   if( dim_space>1 )
   {
      doubleVector const& bottom = exp->doubleVector_data( "vertices_coordinate_1" ) ;
      IBD( 1 ) = bottom.size() ;
      for( size_t i = 0 ; i<bottom.size() ; i++ )
      {
         IJK_COORDS( i, 1 ) = bottom(i) ;
      }
   }
   
   
   if( dim_space==3 )
   {
      doubleVector const& depth = exp->doubleVector_data( "vertices_coordinate_2" ) ;
      IBD( 2 ) = depth.size() ;
      for( size_t j = 0 ; j<depth.size() ; j++ )
      {
         IJK_COORDS( j, 2 ) = depth(j) ;
      }
   }

   // create all colors
   for( size_t icol=0 ; icol<COLOR_TABLE->count() ; icol++ )
   {
      std::string const& n = static_cast<PEL_String const*>( COLOR_TABLE->at( icol ) )->to_string() ;
      GE_Color::extend( n ) ;
      GE_COLOR_TABLE->set_at( icol, const_cast<GE_Color*>(GE_Color::object( n )) ) ;
   }

   for( size_t n=0 ; n<dim_space ; n++ )
   {
      I_POLY(n) = 0 ;
   }
   do
   {
      // might create colors
      corner_color( I_POLY ) ;
      for( size_t n=dim_space-1 ; n<dim_space ; n-- )
      {
         I_POLY(n)++ ;
         if( I_POLY(n) < IBD(n) ) 
            break ;
         else if( n!=0 )
            I_POLY(n) = 0 ;
      }
   } while( I_POLY(0) < IBD(0) ) ;

   // set an invalid state for all iterators
   for( size_t n=0 ; n<dim_space ; n++ )
   {
      I_POLY(n) = IBD(n) ;
   }

   for( size_t n=0 ; n<nb_space_dimensions() ; n++ )
   {
      NB_FACES(n) = IBD(n) ;      
      for( size_t m=0 ; m<nb_space_dimensions() ; m++ )
      {
         if( m!=n )
         {
            NB_FACES(n) *= IBD(m)-1 ;
         }
      }
   }
}

//------------------------------------------------------------------------
size_t
GE_BoxWithBoxes:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: nb_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;
   size_t result = 1 ;
   for( size_t n = 0 ; n<nb_space_dimensions() ; n++ )
   {
      result *= IBD( n ) ;
   }
   return result ;
}

//------------------------------------------------------------------------
size_t
GE_BoxWithBoxes:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: nb_cells" ) ;
   size_t result = 1 ;
   for( size_t n = 0 ; n<nb_space_dimensions() ; n++ )
   {
      result *= IBD( n )-1 ;
   }
      
   return( result ) ;
}

//------------------------------------------------------------------------
size_t
GE_BoxWithBoxes:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: nb_faces" ) ;
   size_t result = 0 ;
   for( size_t n = 0 ; n<nb_space_dimensions() ; n++ )
   {
      result += NB_FACES(n) ;
   }
      
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_BoxWithBoxes:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t n=0 ; n<nb_space_dimensions() ; n++ )
   {
      I_POLY(n) = 0 ;
   }
   VCOORDS.re_initialize( nb_space_dimensions() ) ;
   IVERT = 0 ;
  
   PEL_CHECK_POST( IVERT==vertex_index_ijk( I_POLY ) ) ;
}

//------------------------------------------------------------------------
bool
GE_BoxWithBoxes:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( I_POLY( 0 ) < IBD( 0 ) ) ;
}

//------------------------------------------------------------------------
void
GE_BoxWithBoxes:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: go_next_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   
   for( size_t n=nb_space_dimensions()-1 ; n<nb_space_dimensions() ; n-- )
   {
      I_POLY(n)++ ;
      if( I_POLY(n) < IBD(n) )
      {
         break ;
      }
      else if( n>0 )
      {
         I_POLY(n) = 0 ;
      }
   }
   IVERT++ ;

   PEL_CHECK( IMPLIES( valid_vertex(),
                       IVERT==vertex_index_ijk( I_POLY ) ) ) ;
}

//------------------------------------------------------------------------
doubleVector const& 
GE_BoxWithBoxes:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;

   for( size_t i=0 ; i< nb_space_dimensions(); i++ )
   {
      VCOORDS( i ) = IJK_COORDS( I_POLY(i), i )  ;
   }
    
   doubleVector const& result = VCOORDS ;
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_BoxWithBoxes:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;
   for( size_t n=0 ; n<nb_space_dimensions() ; ++n )
   {
      I_POLY(n) = 0 ;
   }
   CELL_IT = true ;
   for( size_t n=0 ; n<nb_space_dimensions() ; n++ )
   {
      IBD_POLY(n) = IBD(n)-1 ;
   }

   POLY_TO_VERTS.re_initialize( CELL_POLY->nb_vertices() ) ;
   CELL_TO_FACES.re_initialize( CELL_POLY->nb_faces() ) ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_BoxWithBoxes:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: valid_cell" ) ;

   bool result = CELL_IT && ( I_POLY( 0 ) < IBD_POLY( 0 ) ) ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_BoxWithBoxes:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;
   for( size_t n=nb_space_dimensions()-1 ; n<nb_space_dimensions() ; n-- )
   {
      I_POLY(n)++ ;
      if( I_POLY(n) < IBD_POLY(n) )
      {
         break ;
      }
      else if( n>0 )
      {
         I_POLY(n) = 0 ;
      }
   }
}

//------------------------------------------------------------------------
std::string const&
GE_BoxWithBoxes:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = CELL_POLY_NAME ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_BoxWithBoxes:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;
   size_t n = 0 ;
   size_t_vector tmp( I_POLY ) ;
   if( nb_space_dimensions()==1 )
   {
      POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
      tmp( 0 )++ ;
      POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
      tmp( 0 )-- ;
   }
   else
   {
      for( size_t j=1 ; j<nb_space_dimensions() ; j++ )
      {
         if( j==2 )
         {
            tmp( 2 )++ ;
         }
         POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
         tmp( 0 )++ ;
         POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
         tmp( 1 )++ ;
         POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
         tmp( 0 )-- ;
         POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
         tmp( 1 )-- ;
      }
   }
   PEL_CHECK( n==POLY_TO_VERTS.size() ) ;
   size_t_vector const& result = POLY_TO_VERTS ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_BoxWithBoxes:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector tmp( I_POLY ) ;
   size_t n = 0 ;
   if( nb_space_dimensions()==1 )
   {
      CELL_TO_FACES( n++ ) = side_index_ijk( tmp, 0 ) ;
      tmp(0)++ ;
      CELL_TO_FACES( n++ ) = side_index_ijk( tmp, 0 ) ;
      tmp(0)-- ;
   }
   else
   {
      CELL_TO_FACES( n++ ) = side_index_ijk( tmp, 1 ) ;
      tmp(0)++ ;
      CELL_TO_FACES( n++ ) = side_index_ijk( tmp, 0 ) ;
      tmp(0)-- ;
      tmp(1)++ ;
      CELL_TO_FACES( n++ ) = side_index_ijk( tmp, 1 ) ;
      tmp(1)-- ;
      CELL_TO_FACES( n++ ) = side_index_ijk( tmp, 0 ) ;
      if( nb_space_dimensions()==3 )
      {
         CELL_TO_FACES( n++ ) = side_index_ijk( tmp, 2 ) ;
         tmp(2)++ ;
         CELL_TO_FACES( n++ ) = side_index_ijk( tmp, 2 ) ;      
      }
   }
   
   PEL_CHECK( n==CELL_TO_FACES.size() ) ;
   
   size_t_vector const& result = CELL_TO_FACES ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_BoxWithBoxes:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;
   for( size_t n=0 ; n<nb_space_dimensions() ; ++n )
   {
      I_POLY(n) = 0 ;
   }
   IFACE = 0 ;
   CELL_IT = false ;
   for( size_t n=0 ; n<nb_space_dimensions() ; n++ )
   {
      IBD_POLY(n) = IBD(n)-1 ;
   }
   FACE_NORMAL = 0 ;
   IBD_POLY( FACE_NORMAL ) = IBD(FACE_NORMAL) ;

   POLY_TO_VERTS.re_initialize( FACE_POLY->nb_vertices() ) ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
bool
GE_BoxWithBoxes:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: valid_face" ) ;
   
   bool result = ( !CELL_IT ) && ( I_POLY( 0 ) < IBD_POLY( 0 ) ) ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_BoxWithBoxes:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;
   
   for( size_t n=nb_space_dimensions()-1 ; n<nb_space_dimensions() ; n-- )
   {
      I_POLY(n)++ ;
      if( I_POLY(n) < IBD_POLY(n) )
      {
         break ;
      }
      else if( n>0 )
      {
         I_POLY(n) = 0 ;
      }
   }
   if( !valid_face() )
   {
      if( ++FACE_NORMAL < nb_space_dimensions() )
      {
         for( size_t n=0 ; n<nb_space_dimensions() ; n++ )
         {
            IBD_POLY(n) = IBD(n)-1 ;
            I_POLY(n) = 0 ;
         }
         IBD_POLY( FACE_NORMAL ) = IBD(FACE_NORMAL) ;
      }
   }
   if( valid_face() )
   {
      IFACE++ ;      

      PEL_CHECK( IMPLIES( valid_face(),
                 IFACE==side_index_ijk( I_POLY, FACE_NORMAL ) ) ) ;
   }
}

//------------------------------------------------------------------------
std::string const&
GE_BoxWithBoxes:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = FACE_POLY_NAME ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_BoxWithBoxes:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;
   size_t n = 0 ;
   size_t_vector tmp( I_POLY ) ;
   POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
   for( size_t j=0 ; j<nb_space_dimensions() ; j++ )
   {
      if( j!=FACE_NORMAL )
      {
         tmp(j)++ ;
         POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
         for( size_t k=0 ; k<nb_space_dimensions() ; k++ )
         {
            if( k!=FACE_NORMAL && k!=j )
            {
               tmp(k)++ ;
               POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
               tmp(j)-- ;
               POLY_TO_VERTS( n++ ) = vertex_index_ijk( tmp ) ;
            }
         }
         break ;
      }
   }
   
   PEL_CHECK( n==POLY_TO_VERTS.size() ) ;
   size_t_vector const& result = POLY_TO_VERTS ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_BoxWithBoxes:: face_reference_polyhedron( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: face_reference_polyhedron" ) ;
   PEL_CHECK_PRE( face_reference_polyhedron_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_ReferencePolyhedron const* result = FACE_POLY ;

   PEL_CHECK_POST( face_reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_BoxWithBoxes:: cell_reference_polyhedron( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: cell_reference_polyhedron" ) ;
   PEL_CHECK_PRE( cell_reference_polyhedron_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_ReferencePolyhedron const* result = CELL_POLY ;

   PEL_CHECK_POST( cell_reference_polyhedron_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_BoxWithBoxes:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Meshing:: print( os, indent_width ) ;
   
   std::string const s( indent_width+3, ' ' ) ;
   os << s << "nb cells along 0-axis: " << NB_VERTS(0)-1 << std::endl ;
   if( nb_space_dimensions()>1 )
      os << s << "nb cells along 1-axis: " << NB_VERTS(1)-1 << std::endl ;
   if( nb_space_dimensions()>2 )
      os << s << "nb cells along 2-axis: " << NB_VERTS(2)-1 << std::endl ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_BoxWithBoxes:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: default_vertex_color" ) ;
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result  = corner_color( I_POLY ) ;
   if( result==0 )
   {
      result = GE_Color::null_color() ;
   }

   PEL_CHECK( default_vertex_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_BoxWithBoxes:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: default_cell_color" ) ;
   PEL_CHECK( default_cell_color_PRE() ) ;

   GE_Color const* result = GE_Color::null_color() ;

   // Mesh herits color from vertices near boundaries.
   // Size of color name is used to fix prority
   size_t_vector tmp( I_POLY ) ;
   std::string current = "" ;
   for( size_t i=0 ; i<2 ; i++ )
   {
      tmp(0) += i ;
      if( nb_space_dimensions()==1 )
      {
         GE_Color const* col = corner_color( tmp ) ;
         if( col!=0 && col->name().length()>current.length() )
         {
            result = col ;
            current = result->name() ;
         }
      }
      else
      {
         for( size_t j=0 ; j<2 ; j++ )
         {
            tmp(1) += j ;
            GE_Color const* col = corner_color( tmp ) ;
            if( col!=0 && col->name().length()>current.length() )
            {
               result = col ;
               current = result->name() ;
            }
            if( nb_space_dimensions()==3 )
            {
               tmp(2)++ ;
               col = corner_color( tmp ) ;
               if( col!=0 && col->name().length()>current.length() )
               {
                  result = col ;
                  current = result->name() ;
               }
               tmp(2)-- ;
            }
            tmp(1) -= j ;
         }
      }
   }

   PEL_CHECK( default_cell_color_POST( result ) ) ;   
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_BoxWithBoxes:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: default_face_color" ) ;
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = GE_Color::null_color() ;
   if( I_POLY( FACE_NORMAL ) == 0 )
   {
      result =
	 static_cast<GE_Color const*>( GE_COLOR_TABLE->at( FACE_NORMAL ) ) ;
   }
   else if ( I_POLY( FACE_NORMAL ) == IBD_POLY( FACE_NORMAL )-1 )
   {
      result =
	 static_cast<GE_Color const*>(
            GE_COLOR_TABLE->at( nb_space_dimensions() + FACE_NORMAL ) ) ;
   } 
   
   PEL_CHECK( default_face_color_POST( result ) ) ;
   return( result)  ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_BoxWithBoxes:: corner_color( size_t_vector const& idx ) const
//------------------------------------------------------------------------
//????????? en fait :revoit la couleur des bords
{
   PEL_LABEL( "GE_BoxWithBoxes:: corner_color" ) ;
   
   GE_Color const* result = 0 ;
   bool prem = true ;
   string color_name ;
   for( size_t n = nb_space_dimensions()-1 ; n<nb_space_dimensions() ; n-- )
   {
      if( idx(n) == 0 )
      {
         if( !prem ) color_name += "_" ;
         color_name += static_cast<PEL_String const*>( 
                                       COLOR_TABLE->at( n ) )->to_string() ;
         prem = false ;
      }
      else if( idx(n) == IBD(n)-1  )
      {
         if( !prem ) color_name += "_" ;
         color_name += static_cast<PEL_String const*>(
                  COLOR_TABLE->at( nb_space_dimensions() + n ) )->to_string() ;
         prem = false ;
      }
   }
   if( !prem )
   {
      GE_Color::extend( color_name ) ;
      result = GE_Color::object( color_name ) ;
   }
   return( result ) ;
}

//------------------------------------------------------------------------
size_t
GE_BoxWithBoxes:: vertex_index_ijk( size_t_vector const& idx ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: vertex_index_ijk" ) ;
   PEL_CHECK( FORALL( (size_t i=0;i<nb_space_dimensions();i++),
                      idx(i)<IBD(i) ) ) ;
   
   // Global index du sommet de coordonnees discretes i,j,k
   size_t result = 0 ;
   size_t decal = 1 ;
   
   for( size_t n=nb_space_dimensions()-1 ; n<nb_space_dimensions() ; n-- )
   {
      result = result + idx(n)*decal ;
      decal *= IBD( n ) ;
   }
   return( result ) ;
}

//------------------------------------------------------------------------
size_t
GE_BoxWithBoxes:: side_index_ijk( size_t_vector const& idx,
                                  size_t normal ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: side_index_ijk" ) ;
   PEL_CHECK( normal<nb_space_dimensions() ) ;
   
   // Global index de la face de coin inferieur de coordonnees discretes i,j,k
   size_t result = 0 ;
   size_t decal = 1 ;
   
   for( size_t n=nb_space_dimensions()-1 ; n<nb_space_dimensions() ; n-- )
   {
      result += idx(n)*decal ;
      if( n==normal )
      {
         decal *= IBD( n ) ;
      }
      else
      {
         decal *= IBD( n )-1 ;
      }
      
   }
   for( size_t m=0 ; m<normal ; m++ )
   {
      result += NB_FACES( m ) ;
   }
   
   return( result ) ;
}

//-------------------------------------------------------------------------
void
GE_BoxWithBoxes:: check_mesh_polyhedron(
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_BoxWithBoxes:: check_mesh_polyhedron" ) ;
   PEL_CHECK( check_mesh_polyhedron_PRE(
                            dim_space, face_poly_name, cell_poly_name ) ) ;

   // Check the cell reference polyhedron:
   GE_ReferencePolyhedron const* cell_poly_ref = 
                   GE_Mpolyhedron::reference_polyhedron( cell_poly_name ) ;
   if( ! ( ( dim_space == 1 && cell_poly_ref->nb_vertices() == 2 ) || 
           ( dim_space == 2 && cell_poly_ref->nb_vertices() == 4 ) ||
           ( dim_space == 3 && cell_poly_ref->nb_vertices() == 8 ) ) )
   {
      raise_invalid_cell_polyhedron(
         cell_poly_name, "   inconsistent with box meshing" ) ;
   }

   // Check the side reference polyhedron:
   GE_ReferencePolyhedron const* side_poly_ref = 
                   GE_Mpolyhedron::reference_polyhedron( face_poly_name ) ;
   if( ! ( ( dim_space == 1 && side_poly_ref->nb_vertices() == 1 ) ||
           ( dim_space == 2 && side_poly_ref->nb_vertices() == 2 ) ||
           ( dim_space == 3 && side_poly_ref->nb_vertices() == 4 ) ) )
   {
      raise_invalid_face_polyhedron(
         face_poly_name, "   inconsistent with box meshing" ) ;
   }
}
