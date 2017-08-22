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

#include <GE_SplitMeshing.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_GroupExp.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>

#include <iostream>

GE_SplitMeshing const* GE_SplitMeshing:: PROTOTYPE = new GE_SplitMeshing() ;

//-------------------------------------------------------------------------
GE_SplitMeshing:: GE_SplitMeshing( void )
//------------------------------------------------------------------------
   : GE_Meshing( "GE_SplitMeshing" )
   , INITIAL( 0 )
   , SPLIT_STRAT( 0 )
   , NB_VERTS( PEL::bad_index() )
   , LOCAL_VERT(0)
   , IDX_VERT( PEL::bad_index() )
   , NB_FACES( PEL::bad_index() )
   , LOCAL_FACE(0)
   , HALO_FACE(0)
   , IDX_FACE( PEL::bad_index() )
   , NB_CELLS( PEL::bad_index() )
   , LOCAL_CELL(0)
   , IDX_CELL( PEL::bad_index() )
   , ELM_2_VERTS(0)
   , CELL_2_FACES(0)
   , BAND( PEL::bad_index() )
   , HALO( 0 )
{
}

//-------------------------------------------------------------------------
GE_SplitMeshing*
GE_SplitMeshing:: create_replica( PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp,
                                  size_t dim_space ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, dim_space ) ) ;

   GE_SplitMeshing* result = new GE_SplitMeshing(  a_owner, exp, dim_space ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, dim_space ) ) ;
   return result ;
}

//-------------------------------------------------------------------------
GE_SplitMeshing:: GE_SplitMeshing( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp,
                                   size_t dim_space )
//------------------------------------------------------------------------
   : GE_Meshing( a_owner, exp, dim_space )
   , INITIAL(0 )
   , SPLIT_STRAT( 0 )
   , NB_VERTS( PEL::bad_index() )
   , LOCAL_VERT(0)
   , IDX_VERT( PEL::bad_index() )
   , NB_FACES( PEL::bad_index() )
   , LOCAL_FACE(0)
   , HALO_FACE(0)
   , IDX_FACE( PEL::bad_index() )
   , NB_CELLS( PEL::bad_index() )
   , LOCAL_CELL(0)
   , IDX_CELL( PEL::bad_index() )
   , ELM_2_VERTS(0)
   , CELL_2_FACES(0)
   , BAND( exp->int_data( "security_bandwidth" ) )
   , HALO( GE_Color::halo_color() )

{
   PEL_LABEL( "GE_SplitMeshing:: GE_SplitMeshing" ) ;

   // Test that GE_Colorist is not defined in order to not overload halo colors:
   if( exp->has_module( "GE_Colorist" ) )
   {
      PEL_Error::object()->raise_module_error(
         exp,
         "Sub-module \"GE_Colorist\" is avoided in \"GE_SplitMeshing\" module\n"
         "to not overload halo colors" ) ;
   }

   exp->test_data( "security_bandwidth", "security_bandwidth>=0" ) ;
   
   PEL_ModuleExplorer const* se = exp->create_subexplorer( 0, "GE_Meshing" ) ;
   INITIAL = GE_Meshing::create( this, se, dim_space ) ;
   se->destroy() ;

   se = exp->create_subexplorer( 0, "splitting_strategy" ) ;
   bool const has_nb_ranks = exp->has_entry( "nb_ranks" ) ;
   bool const has_rank =  exp->has_entry( "rank" ) ;
   if( has_nb_ranks || has_rank )
   {
      if( !has_nb_ranks || !has_rank )
      {
         PEL_Error::object()->raise_module_error(
            exp,
            "Entries \"nb_ranks\" and \"rank\" have to be specified both" ) ;
      }
      size_t const nb_ranks = exp->int_data( "nb_ranks" ) ;
      size_t const rank = exp->int_data( "rank" ) ;
      exp->test_data( "nb_ranks", "nb_ranks>0" ) ;
      exp->test_data( "rank", "rank>=0 && rank<nb_ranks" ) ;
      if( rank >= nb_ranks )
      {
         PEL_Error::object()->raise_bad_data_value( exp, "rank",
                                                    "0 to nb_ranks-1" ) ;
      }
      SPLIT_STRAT = GE_SplittingStrategy::create( this, se, INITIAL, 
                                                  nb_ranks, rank ) ;
   }
   else
   {
      SPLIT_STRAT = GE_SplittingStrategy::create( this, se, INITIAL, 
                                                  PEL_Exec::communicator() ) ;
   }
   se->destroy() ; se = 0 ;
   
   split_meshing() ;
}

//------------------------------------------------------------------------
GE_SplitMeshing:: ~GE_SplitMeshing( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
size_t
GE_SplitMeshing:: nb_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: nb_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;
  
   return( NB_VERTS ) ;
}

//------------------------------------------------------------------------
size_t
GE_SplitMeshing:: nb_cells( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: nb_cells" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( NB_CELLS ) ;
}

//------------------------------------------------------------------------
size_t
GE_SplitMeshing:: nb_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: nb_faces" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   return( NB_FACES ) ;
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: start_vertex_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: start_vertex_iterator" ) ;
   PEL_CHECK_INV( invariant() ) ;

   IDX_VERT=0 ;
   INITIAL->start_vertex_iterator() ;
   find_next_vertex() ;
}

//------------------------------------------------------------------------
bool
GE_SplitMeshing:: valid_vertex( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: valid_vertex" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( INITIAL->valid_vertex() ) ;
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: go_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: go_next_vertex" ) ;
   PEL_CHECK_PRE( go_next_vertex_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   IDX_VERT++ ;
   INITIAL->go_next_vertex() ;
   find_next_vertex() ;
}

//------------------------------------------------------------------------
doubleVector const& 
GE_SplitMeshing:: vertex_coordinates( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: vertex_coordinates" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_PRE( vertex_coordinates_PRE() ) ;

   doubleVector const& result = INITIAL->vertex_coordinates() ;
   
   PEL_CHECK_POST( vertex_coordinates_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: start_cell_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: start_cell_iterator" ) ;
   PEL_CHECK_PRE( start_cell_iterator_PRE() ) ;

   IDX_CELL=0 ;
   INITIAL->start_cell_iterator() ;
   find_next_cell() ;

   PEL_CHECK_POST( start_cell_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: go_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: go_next_cell" ) ;
   PEL_CHECK_PRE( go_next_cell_PRE() ) ;
   
   IDX_CELL++ ;
   INITIAL->go_next_cell() ;
   find_next_cell() ;
}

//------------------------------------------------------------------------
bool
GE_SplitMeshing:: valid_cell( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: valid_cell" ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = INITIAL->valid_cell() ;

   PEL_CHECK_POST( valid_cell_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
GE_SplitMeshing:: cell_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: cell_polyhedron_name" ) ;
   PEL_CHECK_PRE( cell_polyhedron_name_PRE() ) ;

   std::string const& result = INITIAL->cell_polyhedron_name() ;

   PEL_CHECK_POST( cell_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_SplitMeshing:: cell_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: cell_vertices" ) ;
   PEL_CHECK_PRE( cell_vertices_PRE() ) ;

   size_t_vector const& result = ELM_2_VERTS ;

   PEL_CHECK_POST( cell_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_SplitMeshing:: cell_faces( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: cell_faces" ) ;
   PEL_CHECK_PRE( cell_faces_PRE() ) ;

   size_t_vector const& result = CELL_2_FACES ;

   PEL_CHECK_POST( cell_faces_POST( result ) ) ;
   return( result) ;
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: start_face_iterator( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: start_face_iterator" ) ;
   PEL_CHECK_PRE( start_face_iterator_PRE() ) ;
   
   IDX_FACE = 0 ;
   INITIAL->start_face_iterator() ;
   find_next_face() ;

   PEL_CHECK_POST( start_face_iterator_POST() ) ;
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: go_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: go_next_face" ) ;
   PEL_CHECK_PRE( go_next_face_PRE() ) ;

   IDX_FACE++ ;
   INITIAL->go_next_face() ;
   find_next_face() ;
}

//------------------------------------------------------------------------
bool
GE_SplitMeshing:: valid_face( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: valid_face" ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = INITIAL->valid_face() ;

   PEL_CHECK_POST( valid_face_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
GE_SplitMeshing:: face_polyhedron_name( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: face_polyhedron_name" ) ;
   PEL_CHECK_PRE( face_polyhedron_name_PRE() ) ;

   std::string const& result = INITIAL->face_polyhedron_name() ;

   PEL_CHECK_POST( face_polyhedron_name_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
size_t_vector const& 
GE_SplitMeshing:: face_vertices( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: face_vertices" ) ;
   PEL_CHECK_PRE( face_vertices_PRE() ) ;

   size_t_vector const& result = ELM_2_VERTS ;

   PEL_CHECK_POST( face_vertices_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_SplitMeshing:: default_vertex_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: default_vertex_color" ) ;
   PEL_CHECK( default_vertex_color_PRE() ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Color const* result = INITIAL->vertex_color() ;
   
   PEL_CHECK( default_vertex_color_POST( result ) ) ;
   return result ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_SplitMeshing:: default_cell_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: default_cell_color" ) ;
   PEL_CHECK( default_cell_color_PRE() ) ;
   
   GE_Color const* result = INITIAL->cell_color() ;
   if( SPLIT_STRAT->cell_rank(IDX_CELL) != SPLIT_STRAT->rank() )
   {
      result = HALO ;
   }
   
   PEL_CHECK( default_cell_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
GE_Color const*
GE_SplitMeshing:: default_face_color( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: default_face_color" ) ;
   PEL_CHECK( default_face_color_PRE() ) ;

   GE_Color const* result = INITIAL->face_color() ;
   if( HALO_FACE(IDX_FACE) )
   {
      result = HALO ;
   }   

   PEL_CHECK( default_face_color_POST( result ) ) ;
   return( result ) ;
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: print( std::ostream& os, size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_Meshing:: print( os, indent_width ) ;
   
   std::string const s( indent_width+3, ' ' ) ;
   os << s << "Global meshing:" << std::endl ;
   INITIAL->print( os, indent_width+6 ) ;
   SPLIT_STRAT->print( os, indent_width+3 ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: display( std::ostream& os, size_t indent_width )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: display" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   os << "Initial mesh:" << std::endl ;
   INITIAL->display( os, indent_width+2 ) ;
   os << std::endl ;
   os << "Splitted mesh: rank=" << SPLIT_STRAT->rank() << "/" 
      << SPLIT_STRAT->nb_ranks() << " bandwidth=" << BAND << std::endl ;
   GE_Meshing::display( os, indent_width+2  ) ;
   os << std::endl ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: find_next_vertex( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: find_next_vertex" ) ;
   
   for( ; INITIAL->valid_vertex() ; INITIAL->go_next_vertex() )
   {
      if( LOCAL_VERT( IDX_VERT ) != PEL::bad_index() ) break ;
      IDX_VERT++ ;
   }
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: find_next_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: find_next_cell" ) ;
   
   for( ; INITIAL->valid_cell() ; INITIAL->go_next_cell() )
   {
      if( LOCAL_CELL( IDX_CELL ) != PEL::bad_index() ) break ;
      IDX_CELL++ ;
   }
   if( INITIAL->valid_cell() )
   {
      size_t_vector const& vert = INITIAL->cell_vertices() ;
      size_t_vector const& conn = INITIAL->cell_faces() ;

      ELM_2_VERTS.re_initialize( vert.size() ) ;
      for( size_t i=0 ; i<vert.size() ; i++ )
      {
         ELM_2_VERTS( i ) = LOCAL_VERT( vert( i ) ) ;
      }
      CELL_2_FACES.re_initialize( conn.size() ) ;
      for( size_t i=0 ; i<conn.size() ; i++ )
      {
         CELL_2_FACES( i ) = LOCAL_FACE( conn( i ) ) ;
      }
   }   
}

//------------------------------------------------------------------------
void
GE_SplitMeshing:: find_next_face( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: find_next_face" ) ;
   
   for( ; INITIAL->valid_face() ; INITIAL->go_next_face() )
   {
      if( LOCAL_FACE( IDX_FACE ) != PEL::bad_index() ) break ;
      IDX_FACE++ ;
   }
   if( INITIAL->valid_face() )
   {
      size_t_vector const& vert = INITIAL->face_vertices() ;
      ELM_2_VERTS.re_initialize( vert.size() ) ;
      for( size_t i=0 ; i<vert.size() ; i++ )
      {
         ELM_2_VERTS( i ) = LOCAL_VERT( vert( i ) ) ;
      }
   }   
}

//----------------------------------------------------------------------
void
GE_SplitMeshing:: split_meshing( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_SplitMeshing:: split_meshing" ) ;

   size_t const rank = SPLIT_STRAT->rank() ;
   
   // Cells and vertices:
   {
      NB_VERTS=0 ;
      LOCAL_VERT.resize( INITIAL->nb_vertices() ) ;
      LOCAL_VERT.set( PEL::bad_index() ) ;
      NB_CELLS=0 ;
      LOCAL_CELL.resize( INITIAL->nb_cells() ) ;
      LOCAL_CELL.set( PEL::bad_index() ) ;
      for( size_t b=0 ; b<1+BAND ; b++ )
      {
         size_t m = 0 ;
         INITIAL->start_cell_iterator() ;
         for( ; INITIAL->valid_cell() ; INITIAL->go_next_cell() )
         {
            if( LOCAL_CELL( m ) == PEL::bad_index() )
            {
               bool append_mesh = false ;
               if( b==0 )
               {
                  append_mesh = ( SPLIT_STRAT->cell_rank( m ) == rank ) ;
               }
               else
               {
                  size_t_vector const& verts = INITIAL->cell_vertices() ;
                  for( size_t j=0 ; j<verts.size() && !append_mesh ; j++ )
                  {
                     append_mesh |= ( LOCAL_VERT(verts(j))!=PEL::bad_index()) ;
                  }
               }
               if( append_mesh )
               {
                  LOCAL_CELL( m ) = NB_CELLS++ ;
               }
            }
            ++m ;
         }
         m = 0 ;
         INITIAL->start_cell_iterator() ;
         for( ; INITIAL->valid_cell() ; INITIAL->go_next_cell() )
         {
            if( LOCAL_CELL( m ) != PEL::bad_index() )
            {
               size_t_vector const& verts = INITIAL->cell_vertices() ;
               for( size_t j=0 ; j<verts.size() ; j++ )
               {
                  if( LOCAL_VERT( verts(j) ) == PEL::bad_index() )
                  {
                     LOCAL_VERT( verts(j) ) = NB_VERTS++ ;
                  }
               }
            }
            ++m ;
         }
      }
      if( NB_CELLS==0 )
      {
         PEL_Error::object()->raise_plain(
            "*** GE_SplitMeshing error\n"
            "    Unconsistant splitting : empty submeshing" ) ;
      }
   }

   // Face:
   {
      NB_FACES=0 ;
      LOCAL_FACE.resize( INITIAL->nb_faces() ) ;
      LOCAL_FACE.set( PEL::bad_index() ) ;
      HALO_FACE.resize( INITIAL->nb_faces() ) ;
      HALO_FACE.set( false ) ;
      
      size_t_vector face_2_cell_ini( INITIAL->nb_faces() ) ;
      face_2_cell_ini.set( (size_t) 0 ) ;
      size_t_vector face_2_cell( BAND == 0 ? INITIAL->nb_faces() : 0 ) ;
      face_2_cell.set( (size_t) 0 ) ;
      size_t m = 0 ;
      INITIAL->start_cell_iterator() ;
      for( ; INITIAL->valid_cell() ; INITIAL->go_next_cell() )
      {
         size_t_vector const& connec = INITIAL->cell_faces() ;
         for( size_t i=0 ; i<connec.size() ; i++ )
         {
            size_t j = connec( i ) ;
            face_2_cell_ini( j )++ ;
            if( LOCAL_CELL( m ) != PEL::bad_index() )
            {
               if( BAND == 0 ) face_2_cell( j )++ ;
               if( LOCAL_FACE( j ) == PEL::bad_index() )
               {
                  LOCAL_FACE( j ) = NB_FACES++ ;
                  HALO_FACE( j ) = ( SPLIT_STRAT->cell_rank( m ) != rank ) ;
               }
            }
         }
         ++m ;
      }

      for( size_t i=0 ; i<INITIAL->nb_faces() ; ++i )
      {
         if( face_2_cell_ini(i) == 2 )
         {
            if( BAND == 0 && face_2_cell(i) == 1 )
            {
               HALO_FACE(i) = true ; // ???
            }
         }
         else if( face_2_cell_ini(i) == 1 )
         {
            HALO_FACE(i) = false ; // bound
         }
         else
         {
            PEL_Error::object()->raise_plain(
                  "*** GE_SplitMeshing error\n"
                  "    Unconsistant face connectivity\n"
                  "    of the initial meshing" ) ;
         }
      }
   }

   // Reorder local vertices:
   {
      size_t nbv = 0 ;
      for( size_t idx=0 ; idx<LOCAL_VERT.size() ; ++idx )
      {
         if( LOCAL_VERT(idx) != PEL::bad_index()) LOCAL_VERT(idx) = nbv++ ;
      }
      PEL_ASSERT( nbv == NB_VERTS ) ;
   }

   // Reorder local faces:
   {
      size_t nbs=0 ;
      for( size_t idx=0 ; idx<LOCAL_FACE.size() ; ++idx )
      {
         if( LOCAL_FACE(idx) != PEL::bad_index()) LOCAL_FACE(idx) = nbs++ ;
      }
      PEL_ASSERT( nbs == NB_FACES ) ;
   }
}
