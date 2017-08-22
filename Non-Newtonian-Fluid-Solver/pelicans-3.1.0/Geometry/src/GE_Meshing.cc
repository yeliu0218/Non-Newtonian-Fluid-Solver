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

#include <GE_Meshing.hh>

#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Variable.hh>
#include <PEL.hh>

#include <GE_Color.hh>
#include <GE_Colorist.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>

#include <size_t_array2D.hh>
#include <stringVector.hh>

#include <ios>
#include <iostream>
#include <sstream>
#include <iomanip>

using std::cerr ;
using std::endl ;
using std::ios_base ;
using std::setprecision ; 
using std::setw ;
using std::string ;

//-------------------------------------------------------------------------
GE_Meshing*
GE_Meshing:: create( PEL_Object* a_owner,
                     PEL_ModuleExplorer const* exp,
                     size_t dim_space )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( dim_space==1 || dim_space==2 || dim_space==3 ) ;

   GE_Meshing* result = 0 ;

   string name = exp->string_data( "concrete_name" ) ;
   GE_Meshing const* proto =
      static_cast<GE_Meshing const*>( plugins_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
   
   result = proto->create_replica( a_owner, exp, dim_space ) ;
   
   if( result->COLORIST != 0 ) result->COLORIST->initialize( result ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_space_dimensions()==dim_space ) ;
   return result ;
}

//-------------------------------------------------------------------------
GE_Meshing:: GE_Meshing( PEL_Object* a_owner,
                         PEL_ModuleExplorer const* exp,
                         size_t dim_space  )
//------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , DIM( dim_space )
   , IS_PROTO( false )
   , COLORIST( 0 )
   , RDOFF( 0 )
   , XX( 0 )
{
   PEL_LABEL( "GE_Meshing:: GE_Meshing" ) ;
   if( exp->has_module( "GE_Colorist" ) )
   {
      PEL_ModuleExplorer const* sexp = 
                           exp->create_subexplorer( 0, "/GE_Colorist" ) ;
      COLORIST = GE_Colorist::create( this, sexp ) ;
      sexp->destroy() ;
   }
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//-------------------------------------------------------------------------
GE_Meshing:: GE_Meshing( std::string const& name )
//------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , DIM( 0 )
   , IS_PROTO( true )
   , COLORIST( 0 )
   , RDOFF( 0 )
   , XX( 0 )
{
   PEL_LABEL( "GE_Meshing:: GE_Meshing" ) ;
   
   plugins_map()->register_item( name, this ) ;
   
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//------------------------------------------------------------------------
GE_Meshing:: ~GE_Meshing( void )
//------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
GE_Meshing:: nb_space_dimensions( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: nb_space_dimensions" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return DIM ;
}

//-------------------------------------------------------------------------
GE_Color const*
GE_Meshing:: vertex_color( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: vertex_color" ) ;
   PEL_CHECK_PRE( valid_vertex() ) ;

   PEL_CHECK( !is_a_prototype() ) ;
   
   GE_Color const* result = 0 ;
   if( COLORIST != 0 )
   {
      result = COLORIST->vertex_color( vertex_coordinates() ) ;
   }

   if( result == 0 )
   {
      result = default_vertex_color() ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   return result ;
}

//-------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Meshing:: cell_reference_polyhedron( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: cell_reference_polyhedron" ) ;
   PEL_CHECK_PRE( cell_reference_polyhedron_PRE() ) ;

   GE_ReferencePolyhedron const* result = 
            GE_Mpolyhedron::reference_polyhedron( cell_polyhedron_name() ) ;

   PEL_CHECK_POST( cell_reference_polyhedron_POST( result ) ) ;
   return result ;
}

//-------------------------------------------------------------------------
GE_Color const*
GE_Meshing:: cell_color( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: cell_color" ) ;
   PEL_CHECK_PRE( valid_cell() ) ;

   PEL_CHECK( !is_a_prototype() ) ;
   
   GE_Color const* result = 0 ;
   if( COLORIST != 0 )
   {
      result = COLORIST->cell_color( cell_vertices() ) ;
   }

   if( result == 0 )
   {
      result = default_cell_color() ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   return result ;
}

//-------------------------------------------------------------------------
GE_ReferencePolyhedron const*
GE_Meshing:: face_reference_polyhedron( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: face_reference_polyhedron" ) ;
   PEL_CHECK_PRE( face_reference_polyhedron_PRE() ) ;

   GE_ReferencePolyhedron const* result = 
            GE_Mpolyhedron::reference_polyhedron( face_polyhedron_name() ) ;

   PEL_CHECK_POST( face_reference_polyhedron_POST( result ) ) ;
   return result ;
}

//-------------------------------------------------------------------------
GE_Color const*
GE_Meshing:: face_color( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: face_color" ) ;
   PEL_CHECK_PRE( valid_face() ) ;

   PEL_CHECK( !is_a_prototype() ) ;
   
   GE_Color const* result = 0 ;
   if( COLORIST != 0 )
   {
      result = COLORIST->face_color( face_vertices() ) ;
   }

   if( result == 0 )
   {
      result = default_face_color() ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   return result ;
}

//------------------------------------------------------------------------
void
GE_Meshing:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: print" ) ;

   string space( indent_width, ' ' ) ;
   os << space << type_name() << std::endl ;   
}

//------------------------------------------------------------------------
void
GE_Meshing:: display( std::ostream& os, size_t indent_width )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: display" ) ;

   string space( indent_width, ' ' ) ;

   os << space << "Space dimension " << nb_space_dimensions() << endl ;
   os << space << nb_vertices()  << " vertices : " << endl ;
   os << space << "n° : coordinates : color " << endl ;
   size_t idx = 0 ;
   for( start_vertex_iterator() ; valid_vertex() ; go_next_vertex() )
   {
      os << space << idx << " : " << vertex_coordinates() ;
      os << " : " << vertex_color()->name() << endl ;
      idx++ ;
   }
   
   os << space << nb_faces()  << " sides : " << endl ;
   os << space << "n° : vertices : color " << endl ;
   idx = 0 ;
   for( start_face_iterator() ; valid_face() ; go_next_face() )
   {
      os << space << idx << " : " << face_polyhedron_name() << " : " 
          << face_vertices() ;
      os << " : " << face_color()->name() << endl ;
      idx++ ;
   }
   
   os << space << nb_cells()  << " cells : " << endl ;
   os << space << "n° : vertices : connectivity : color" << endl ;
   idx = 0 ;
   for( start_cell_iterator() ; valid_cell() ; go_next_cell() )
   {
      os << space << idx << " : " << cell_polyhedron_name() << " : " 
          << cell_vertices()
          << " : " 
          << cell_faces() ;
      os << " : " << cell_color()->name() << endl ;
      idx++ ;
   }
}

//-----------------------------------------------------------------------------
void
GE_Meshing:: print_as_an_explicit_meshing( std::ostream& os )
//-----------------------------------------------------------------------------
{
   os << " concrete_name = \"GE_ExplicitMeshing\" " << endl ;

   // Vertices
   os << " MODULE vertices" << endl ;
   os << "    number = " << nb_vertices() << endl ;
   start_vertex_iterator() ;
   os << "    coordinates = <" << endl ;
   ios_base::fmtflags original_flags = os.flags() ;
   os.setf( ios_base::uppercase | ios_base::scientific ) ;
   std::streamsize p = os.precision() ;
   os << setprecision( 10 ) ;
   while( valid_vertex() )
   {
      for( size_t i=0; i<nb_space_dimensions(); i++ )
      {
	 os << setw( 23 ) << vertex_coordinates()( i ) << " " ;
      }
      os << endl ;
      go_next_vertex() ;
   }
   os << setprecision( p ) ;
   os.flags( original_flags ) ;
   os << "                  >" << endl ;
   size_t_array2D vert_color( (size_t)50, nb_vertices() ) ;
   stringVector colors( (size_t)0 ) ;
   start_vertex_iterator() ;
   size_t counter = 0 ;
   while( valid_vertex() )
   {
      std::string color_name = vertex_color()->name() ;
      if( !colors.has( color_name ) )
      {
	 colors.append( color_name ) ;
      }
      vert_color( colors.index_of( color_name ), counter ) = 1 ;
      go_next_vertex() ;
      counter++ ;
   }
   print_colors( os, colors, vert_color ) ;
   os << " END MODULE vertices " << endl ;

   // Sides
   stringVector side_names(0) ;
   for( start_face_iterator() ; valid_face() ; go_next_face() )
   {
      if( !side_names.has( face_polyhedron_name() ) )
      {
         side_names.append( face_polyhedron_name() ) ;
      }
   }
   os << " MODULE faces " << endl ;
   os << "    number = " << nb_faces() << endl ;
   os << "    MODULE polyhedra_and_connectivities " << endl ;
   for( size_t i_side_name=0 ; i_side_name<side_names.size() ; ++i_side_name )
   {
      os << "       " << side_names(i_side_name) << "  = < " << endl ;
      for( start_face_iterator() ; valid_face() ; go_next_face() )
      {
         if( face_polyhedron_name()==side_names(i_side_name) )
         {
            size_t_vector const& vert = face_vertices() ;
            os << "          " ;
            for( size_t i=0; i<vert.size(); i++ )
            {
               os << vert(i) << " " ;
            }
            os << endl ;
         }
      }
      os << "         > " << endl ;
   }
   os << "    END MODULE polyhedra_and_connectivities " << endl ;
   size_t_array2D sid_color( (size_t)50, nb_faces() ) ;
   colors.re_initialize( (size_t)0 ) ;
   start_face_iterator() ;
   counter = 0 ;
   while( valid_face() )
   {
      std::string color_name = face_color()->name() ;
      if( !colors.has( color_name ) )
      {
	 colors.append( color_name ) ;
      }
      sid_color( colors.index_of( color_name ), counter ) = 1 ;
      go_next_face() ;
      counter++ ;
   }
   print_colors( os, colors, sid_color ) ;
   os << " END MODULE faces " << endl ;

   // Cells
   stringVector cell_names(0) ;
   for( start_cell_iterator() ; valid_cell() ; go_next_cell() )
   {
      if( !cell_names.has( cell_polyhedron_name() ) )
      {
         cell_names.append( cell_polyhedron_name() ) ;
      }
      
   }
   os << " MODULE cells " << endl ;
   os << "    number = " << nb_cells() << endl ;
   os << "    MODULE polyhedra_and_connectivities " << endl ;
   for( size_t i_cell_name=0 ; i_cell_name<cell_names.size() ; ++i_cell_name )
   {
      os << "       " << cell_names(i_cell_name) << "  = < " << endl ;
      for( start_cell_iterator() ; valid_cell() ; go_next_cell() )
      {
         if( cell_polyhedron_name()==cell_names(i_cell_name) )
         {
            size_t_vector const& verti = cell_vertices() ;
            os << "          " ;
            for( size_t i=0; i<verti.size(); i++ )
            {
               os << verti(i) << " " ;
            }
            os << endl ;
            os << "          " ;
            size_t_vector const& sidei = cell_faces() ;
            for( size_t i=0; i<sidei.size(); i++ )
            {
               os << sidei(i) << " " ;
            }
            os << endl ;          
         }
      }
      os << "         > " << endl ;
   }
   os << "    END MODULE polyhedra_and_connectivities " << endl ;
   size_t_array2D cel_color( (size_t)50, nb_cells() ) ;
   colors.re_initialize( (size_t)0 ) ;
   start_cell_iterator() ;
   counter = 0 ;
   while( valid_cell() )
   {
      std::string color_name = cell_color()->name() ;
      if( !colors.has( color_name ) )
      {
	 colors.append( color_name ) ;
      }
      cel_color( colors.index_of( color_name ), counter ) = 1 ;
      go_next_cell() ;
      counter++ ;
   }
   print_colors( os, colors, cel_color ) ;
   os << " END MODULE cells" << endl ;
}

//---------------------------------------------------------------------------
void
GE_Meshing:: print_colors( std::ostream& os,
                           stringVector const& colors, 
                           size_t_array2D const& item_color )
//---------------------------------------------------------------------------
{
   size_t nb_items = item_color.index_bound( 1 ) ;

   if( colors.size()>0 )
   {
      bool col = false ;
      size_t halo_color = PEL::bad_index() ;
      size_t null_color = PEL::bad_index() ;
      for( size_t i=0; i<colors.size(); i++ )
      {
         if( colors(i)==GE_Color::halo_color()->name() )
         {
            PEL_ASSERT( halo_color==PEL::bad_index() ) ;
            halo_color = i ;
         }
         else if( colors(i)==GE_Color::null_color()->name() )
         {
            PEL_ASSERT( null_color==PEL::bad_index() ) ;
            null_color = i ;
         }
         else
         {
            col = true ;
         }
      }
      if( col )
      {
         os << "    MODULE colors" << endl ;
         for( size_t i=0; i<colors.size(); i++ )
         {
            if( i!=halo_color && i!=null_color )
            {
               os << "       " << colors(i) << " = < " ;
               for( size_t j=0 ; j<nb_items ; ++j )
               {
                  if( item_color( i, j )==1 ) os << j << " " ;
               }
               os << ">" << endl ;
            }
         }
         os << "    END MODULE colors" << endl ;
      }
      if( halo_color!=PEL::bad_index() || null_color!=PEL::bad_index() )
      {
         os << "    MODULE special_colors" << endl ;
         if( null_color!=PEL::bad_index() )
         {
            os << "       null = < " ;
            for( size_t j=0; j<nb_items; j++ )
            {
               if( item_color( null_color, j )==1 ) os << j << " " ;
            }
            os << ">" << endl ;
         }
         if( halo_color!=PEL::bad_index() )
         {
            os << "       halo = < " ;
            for( size_t j=0; j<nb_items; j++ )
            {
               if( item_color( halo_color, j )==1 ) os << j << " " ;
            }
            os << ">" << endl ;
         }
         os << "    END MODULE special_colors" << endl ;
      }
   }
}

//-------------------------------------------------------------------------
bool
GE_Meshing:: is_a_prototype( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: is_a_prototype" ) ;
   return IS_PROTO ;
}

//-------------------------------------------------------------------------
void
GE_Meshing:: read_mesh_polyhedron( PEL_ModuleExplorer const* exp,
                                   size_t dim_space,
                                   std::string& face_poly_name,
                                   std::string& cell_poly_name ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: read_mesh_polyhedron" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( dim_space == 1 || dim_space == 2 || dim_space == 3 ) ;
   
   stringVector const& vec = exp->stringVector_data( "mesh_polyhedron" ) ;
   if( vec.size() != 2 )
      raise_invalid_mesh_polyhedron( "   bad number of arguments" ) ;
   face_poly_name = vec( 0 ) ;
   cell_poly_name = vec( 1 ) ;
  
   GE_ReferencePolyhedron const* face_ref = 
                   GE_Mpolyhedron::reference_polyhedron( face_poly_name ) ;
   GE_ReferencePolyhedron const* cell_ref =
                   GE_Mpolyhedron::reference_polyhedron( cell_poly_name ) ;
   
   // Check the cell reference polyhedron:
   if( cell_ref->dimension() != dim_space )
   {
      raise_invalid_cell_polyhedron(
         cell_poly_name, "   inconsistent with space dimension" ) ;
   }

   // Check the face reference polyhedron:
   if( face_ref->dimension()+1 != dim_space )
   {
      raise_invalid_face_polyhedron(
         face_poly_name, "   inconsistent with space dimension" ) ;
   }

   // Check cell/face consistency:
   for( size_t i=0 ; i<cell_ref->nb_faces() ; ++i )
   {
      if( cell_ref->nb_face_vertices( i ) != face_ref->nb_vertices() )
      {
         raise_invalid_mesh_polyhedron( "   incompatible polyhedra" ) ;
      }
   }
   
   // Particular checks
   check_mesh_polyhedron( dim_space, face_poly_name, cell_poly_name ) ;

   PEL_CHECK_POST( !face_poly_name.empty() ) ;
   PEL_CHECK_POST( !cell_poly_name.empty() ) ;
}

//-------------------------------------------------------------------------
void
GE_Meshing:: check_mesh_polyhedron( size_t dim_space,
                                    std::string const& face_poly_name,
                                    std::string const& cell_poly_name ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: check_mesh_polyhedron" ) ;
   PEL_CHECK_PRE( check_mesh_polyhedron_PRE(
                            dim_space, face_poly_name, cell_poly_name ) ) ;
}

//-------------------------------------------------------------------------
void
GE_Meshing:: raise_invalid_nb_space_dimensions( size_t dim,
                                                std::string allowed_dims ) const
//-------------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << type_name() << " : " 
        << endl ;
   mesg << "   invalid number of space dimensions : " << dim 
        << endl ;
   mesg << "   allowed number of space dimensions : " << allowed_dims
	<< endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//-------------------------------------------------------------------------
void
GE_Meshing:: raise_invalid_mesh_polyhedron( std::string const& error ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: raise_invalid_mesh_polyhedron" ) ;
   PEL_CHECK_PRE( !error.empty() ) ;
   
   std::ostringstream mesg ;
   mesg << type_name() << " : " << endl ;
   mesg << "   \"mesh_polyhedron\" : " << error << endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//-------------------------------------------------------------------------
void
GE_Meshing:: raise_invalid_cell_polyhedron( 
             std::string const& poly_name, std::string const& error ) const
//-------------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << type_name() << " : " << endl ;
   mesg << "   invalid cell polyhedron : " << poly_name << endl ;
   if( !error.empty() ) mesg << error ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//-------------------------------------------------------------------------
void
GE_Meshing:: raise_invalid_face_polyhedron( 
             std::string const& poly_name, std::string const& error ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: raise_invalid_face_polyhedron" ) ;
   PEL_CHECK_PRE( !poly_name.empty() ) ;
   
   std::ostringstream mesg ;
   mesg << type_name() << " : " << endl ;
   mesg << "   invalid side polyhedron : " << poly_name << endl ;
   if( !error.empty() ) mesg << error ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//-------------------------------------------------------------------------
void
GE_Meshing:: initialize_rounding_strategy( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: initialize_rounding_strategy" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   if( exp->has_entry( "roundoff" ) )
   {
      PEL_ContextSimple* ct = PEL_ContextSimple::create( this ) ;
      XX = PEL_Double::create( ct, 0.0 ) ;
      ct->extend( PEL_Variable::object( "DS_x"), XX ) ;
      RDOFF = exp->abstract_data( this, "roundoff", ct ) ;

      if( !RDOFF->value_can_be_evaluated() )
      {
         PEL_Error::object()->raise_not_evaluable(
                     exp, "roundoff", RDOFF->undefined_variables() ) ;
      }
      if( RDOFF->data_type() != PEL_Data::Double )
      {
         PEL_Error::object()->raise_bad_data_type(
                                 exp, "roundoff", PEL_Data::Double ) ;
      }
   }
}

//-------------------------------------------------------------------------
double
GE_Meshing:: roundoff( double x ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing:: roundoff" ) ;

   double result = x ;

   if( RDOFF != 0 )
   {
      XX->set( x ) ;
      result = RDOFF->to_double() ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: go_next_vertex_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_vertex() ) ;
   return true ;
}

//-------------------------------------------------------------------------
bool
GE_Meshing:: vertex_coordinates_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( valid_vertex() ) ;
   return true ;
}

//-------------------------------------------------------------------------
bool
GE_Meshing:: vertex_coordinates_POST( doubleVector const& result ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result.size() == nb_space_dimensions() ) ;
   return true ;
}

//-------------------------------------------------------------------------
bool
GE_Meshing:: default_vertex_color_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( !is_a_prototype() ) ;
   PEL_ASSERT( valid_vertex() ) ;
   return true ;
}

//-------------------------------------------------------------------------
bool
GE_Meshing:: default_vertex_color_POST( GE_Color const* result ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: start_cell_iterator_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( !valid_face() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: start_cell_iterator_POST( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_cell() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: go_next_cell_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_cell() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: valid_cell_POST( bool result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result==true, valid_face()==false ) ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: cell_polyhedron_name_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_cell() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: cell_polyhedron_name_POST( std::string const& result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron( result )->dimension() <= 
               nb_space_dimensions() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: cell_reference_polyhedron_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_cell() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: cell_reference_polyhedron_POST(
                                GE_ReferencePolyhedron const* result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( result ==
          GE_Mpolyhedron::reference_polyhedron( cell_polyhedron_name() ) ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: default_cell_color_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( !is_a_prototype() ) ;
   PEL_ASSERT( valid_cell() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: default_cell_color_POST( GE_Color const* result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: cell_vertices_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_cell() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: cell_vertices_POST( size_t_vector const& result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( result.size() == GE_Mpolyhedron::reference_polyhedron( cell_polyhedron_name() )->nb_vertices() ) ;
   PEL_ASSERT( FORALL( ( size_t i=0 ; i<result.size() ; ++i ),
                       result(i) < nb_vertices() ) ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: cell_faces_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_cell() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: cell_faces_POST( size_t_vector const& result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( result.size() == GE_Mpolyhedron::reference_polyhedron( cell_polyhedron_name() )->nb_faces() ) ;
   PEL_ASSERT( FORALL( ( size_t i=0 ; i<result.size() ; ++i ),
                  result(i) < nb_faces() ) ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: start_face_iterator_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( !valid_cell() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: start_face_iterator_POST( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_face() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: go_next_face_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_face() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: valid_face_POST( bool result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result==true, valid_cell()==false ) ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: face_polyhedron_name_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_face() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: face_polyhedron_name_POST( std::string const& result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( GE_Mpolyhedron::reference_polyhedron( result )->dimension() <= 
               nb_space_dimensions()-1 ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: face_reference_polyhedron_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_face() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: face_reference_polyhedron_POST(
                                GE_ReferencePolyhedron const* result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( result ==
          GE_Mpolyhedron::reference_polyhedron( face_polyhedron_name() ) ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: default_face_color_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( !is_a_prototype() ) ;
   PEL_ASSERT( valid_face() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: default_face_color_POST( GE_Color const* result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: face_vertices_PRE( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( valid_face() ) ;
   return true ;
}

//----------------------------------------------------------------------------
bool
GE_Meshing:: face_vertices_POST( size_t_vector const& result ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( result.size() == 
               GE_Mpolyhedron::reference_polyhedron( 
                            face_polyhedron_name() )->nb_vertices() ) ;
   PEL_ASSERT( FORALL( ( size_t i=0 ; i<result.size() ; ++i ),
                       result(i) < nb_vertices() ) ) ;
   return true ;
}

//-------------------------------------------------------------------------
bool
GE_Meshing:: check_mesh_polyhedron_PRE( 
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const 
//------------------------------------------------------------------------
{
   PEL_ASSERT( dim_space == 1 || dim_space == 2 || dim_space == 3 ) ;
   PEL_ASSERT( !face_poly_name.empty() ) ;
   PEL_ASSERT( !cell_poly_name.empty() ) ;
   return( true ) ;   
}

//-------------------------------------------------------------------------
bool
GE_Meshing:: create_replica_PRE( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp,
                                 size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( exp!=0 ) ;
   PEL_ASSERT( dim_space==1 || dim_space==2 || dim_space==3 ) ;
   return( true ) ;   
}

//-------------------------------------------------------------------------
bool
GE_Meshing:: create_replica_POST( GE_Meshing const* result,
                                  PEL_Object* a_owner,
                                  size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   PEL_ASSERT( result->owner()==a_owner ) ;
   PEL_ASSERT( result->nb_space_dimensions()==dim_space ) ;
   return( true ) ;   
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
GE_Meshing:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
                 PEL_ObjectRegister::create( PEL_Root::object(),
                                             "GE_Meshing descendant" ) ;
   return( result ) ;
}
