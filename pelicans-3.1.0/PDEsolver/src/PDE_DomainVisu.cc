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

#include <PDE_DomainVisu.hh>

#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_System.hh>
#include <boolVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_SetOfPoints.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <sstream>
#include <string>

using std::cout ; using std::endl ;
using std::ios_base ; using std::setprecision ; using std::setw ;

PDE_DomainVisu const* PDE_DomainVisu::PROTOTYPE = new PDE_DomainVisu() ;

struct PDE_DomainVisu_ERROR
{
   static void n0( std::string const& file_name ) ;
} ;

// Real names "halo color" and "null color" with a blank are unreadable by gmv
std::string const HALO_COLOR = "halo_color" ;
std::string const NULL_COLOR = "null_color" ;

//----------------------------------------------------------------------
PDE_DomainVisu:: PDE_DomainVisu( void )
//----------------------------------------------------------------------
   : PEL_Application( "PDE_DomainVisu" ) 
   , TRACE_CELLS( false )
   , TRACE_SIDES( false )
   , TRACE_BOUNDS( false )
   , CHECK_CONN(  false )
   , NB_DIMS( PEL::bad_index() )
   , VERTICES( 0 )
   , bFE( 0 )
   , cFE( 0 )
   , sFE( 0 )
   , SAVER( 0 )
   , NBS( PEL::bad_index() )
   , GMV_FILES( false )
{
}

//----------------------------------------------------------------------
PDE_DomainVisu*
PDE_DomainVisu:: create_replica( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PDE_DomainVisu* result = new PDE_DomainVisu( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DomainVisu*
PDE_DomainVisu:: create( PEL_Object* a_owner,
         	         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: create_replica" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   PDE_DomainVisu* result = new PDE_DomainVisu( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DomainVisu:: PDE_DomainVisu( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , TRACE_CELLS( false )
   , TRACE_SIDES( false )
   , TRACE_BOUNDS( false )
   , CHECK_CONN( exp->has_entry( "check_meshing_connexity" ) ?
                   exp->bool_data( "check_meshing_connexity" ) : false )
   , NB_DIMS( PEL::bad_index() )
   , VERTICES( 0 )
   , bFE( 0 )
   , cFE( 0 )
   , sFE( 0 )
   , SAVER( 0 )
   , NBS( PEL::bad_index() )
   , GMV_FILES( exp->has_entry( "build_GMV_files" ) ?
                   exp->bool_data( "build_GMV_files" ) : true )
{
   PEL_LABEL( "PDE_DomainVisu:: PDE_DomainVisu" ) ;

   size_t before = PEL_System::used_memory() ;
   PEL::out() << "Memory before meshing built: " << before << endl << endl ;
   PEL_ModuleExplorer* e =
      exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields const* dom = PDE_DomainAndFields::create( this, e ) ;
   e->destroy() ;
   size_t after = PEL_System::used_memory() ;
   PEL::out() << "Memory after meshing built: " << after << endl ;
   PEL::out() << "   difference: " << (after-before)/1024/1024 << " Mo" 
              << endl << endl ;
   
   SAVER = dom->result_saver() ;

   if( GMV_FILES )
   {
      NB_DIMS = dom->nb_space_dimensions() ;
      VERTICES = dom->set_of_vertices() ;
      bFE = dom->create_LocalFEbound( this ) ;
      bFE->include_color( GE_Color::halo_color() ) ;
      cFE = dom->create_LocalFEcell( this ) ;
      cFE->include_color( GE_Color::halo_color() ) ;
      sFE = dom->create_CursorFEside( this ) ;
      sFE->include_color( GE_Color::halo_color() ) ;

      NBS = 0 ;
      for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
      {
         NBS += ( sFE->is_periodic() ? 2 : 1 ) ;
      }
   
      if( exp->has_module( "traces" ) )
      {
         PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "traces" ) ;
         OFS.open( se->string_data( "output_file" ).c_str() ) ;
         if( se->has_entry( "trace_cells" ) )
            TRACE_CELLS = se->bool_data( "trace_cells" ) ;
         if( se->has_entry( "trace_sides" ) )
            TRACE_SIDES = se->bool_data( "trace_sides" ) ;
         if( se->has_entry( "trace_bounds" ) )
            TRACE_BOUNDS = se->bool_data( "trace_bounds" ) ;
         se->destroy() ;
      
         PDE_SetOfDiscreteFields* sdf = dom->set_of_discrete_fields() ;
         for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
         {
            PDE_DiscreteField const* ff = sdf->item() ;
            if( !dom->is_defined_on_boundary( ff->name() ) )
            {
               if( TRACE_CELLS ) cFE->require_field_calculation( ff, PDE_LocalFE::N ) ;
               if( TRACE_SIDES ) sFE->require_field_calculation( ff, PDE_LocalFE::N ) ;
            }
            if( TRACE_BOUNDS ) bFE->require_field_calculation( ff, PDE_LocalFE::N ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
PDE_DomainVisu:: ~PDE_DomainVisu( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: run( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: run" ) ;

   if( GMV_FILES )
   {
      PEL::out() << "*** Building gmv files" << endl ;
      build_file_with_bounds() ;
      build_file_with_faces() ;
      build_file_with_cells() ;
      if( CHECK_CONN ) check_connexity() ;
   }
   
   SAVER->start_cycle() ;
   SAVER->save_grid() ;
   SAVER->save_fields( 0 ) ;
   SAVER->terminate_cycle() ;
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: build_file_with_bounds( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: build_file_with_bounds" ) ;

   std::string fname = "bounds.gmv" ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   if( com->nb_ranks() > 1 )
   {
      std::ostringstream m ;
      m << "bounds_" << com->rank() << ".gmv" ;
      fname = m.str() ;
   }
   std::ofstream ofile( fname.c_str(),
                        std::ios::out | std::ios::trunc ) ;
   if( !ofile )
   {
      PDE_DomainVisu_ERROR::n0( fname ) ;
   }
   
   GE_SetOfPoints* skin_vertices = GE_SetOfPoints::create( 0, NB_DIMS ) ;
   size_t_vector glob_cell_2_skin_cell( cFE->nb_meshes() ) ;
   size_t nb_skin_cells ;
   build_skin( skin_vertices, nb_skin_cells, glob_cell_2_skin_cell ) ;

   ofile << "gmvinput ascii" << endl ;
   write_gmv_nodev( ofile, skin_vertices ) ;
   write_gmv_faces( ofile, 
                    skin_vertices, nb_skin_cells, glob_cell_2_skin_cell ) ;
   write_gmv_surface( ofile, skin_vertices, false ) ;
   write_gmv_surfvars( ofile, false ) ;

   ofile << "endgmv" << endl ;
   ofile.close() ;
   skin_vertices->destroy() ;
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: build_skin( GE_SetOfPoints* skin_vertices,
                             size_t& nb_skin_cells,
                             size_t_vector& glob_cell_2_skin_cell )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: build_skin" ) ;
   
   glob_cell_2_skin_cell.set( PEL::bad_index() ) ;
   
   nb_skin_cells = 0 ;
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      GE_Mpolyhedron const* poly = bFE->polyhedron() ;
      for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
      {
         skin_vertices->extend( poly->vertex( i ) ) ;
      }
      size_t ii = bFE->adjacent_cell_id() ;
      if( glob_cell_2_skin_cell( ii ) == PEL::bad_index() )
      {
         glob_cell_2_skin_cell( ii ) = nb_skin_cells ;
         ++nb_skin_cells ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: build_file_with_faces( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: build_file_with_faces" ) ;

   std::string fname = "faces.gmv" ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   if( com->nb_ranks() > 1 )
   {
      std::ostringstream m ;
      m << "faces_" << com->rank() << ".gmv" ;
      fname = m.str() ;
   }
   std::ofstream ofile( fname.c_str(),
                        std::ios::out | std::ios::trunc ) ;
   if( !ofile )
   {
      PDE_DomainVisu_ERROR::n0( fname ) ;
   }

   ofile << "gmvinput ascii" << endl ;
   write_gmv_nodev( ofile, VERTICES ) ;
   write_gmv_faces( ofile ) ;
   write_gmv_surface( ofile, VERTICES, true ) ;
   write_gmv_surfvars( ofile, true ) ;

   ofile << "endgmv" << endl ;
   ofile.close() ;
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: build_file_with_cells( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: build_file_with_cells" ) ;

   std::string fname = "cells.gmv" ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   if( com->nb_ranks() > 1 )
   {
      std::ostringstream m ;
      m << "cells_" << com->rank() << ".gmv" ;
      fname = m.str() ;
   }
   std::ofstream ofile( fname.c_str(),
                        std::ios::out | std::ios::trunc ) ;
   if( !ofile )
   {
      PDE_DomainVisu_ERROR::n0( fname ) ;
   }

   ofile << "gmvinput ascii" << endl ;
   write_gmv_nodev( ofile, VERTICES ) ;
   write_gmv_cells( ofile ) ;
   write_gmv_variable( ofile ) ;

   ofile << "endgmv" << endl ;
   ofile.close() ;
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: write_gmv_nodev( std::ostream& ofile,
                                  GE_SetOfPoints const* verts) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: write_gmv_nodev" ) ;

   ofile << "nodev " << verts->nb_points() << endl ;
   for( size_t iv=0 ; iv<verts->nb_points() ; ++iv )
   {
      GE_Point const* pt = verts->point( iv ) ;
      for( size_t d=0 ; d<NB_DIMS ; ++d )
      {
         ofile << pt->coordinate( d ) << " " ;
      }
      if( NB_DIMS == 2 ) ofile << "0" ;
      if( NB_DIMS == 1 ) ofile << "0 0" ;
      ofile << endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: write_gmv_faces( std::ostream& ofile )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: write_gmv_faces" ) ;

   ofile << "faces " << NBS + bFE->nb_meshes() << " "
         << cFE->nb_meshes() << endl ;
   PDE_LocalFEcell const* cfe_0 = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* cfe_1 = sFE->adjacent_localFEcell( 1 ) ;
   if( TRACE_SIDES )
   {
      OFS << "***************************************" << endl ;
      OFS << "PDE_CursorFEside : " << sFE->nb_meshes() << " sides" << endl ;
      OFS << "***************************************" << endl ;
   }
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      for( size_t ipoly=0 ; ipoly<( sFE->is_periodic() ? 2 : 1 ) ; ipoly++ )
      {
         GE_Mpolyhedron const* poly = ( sFE->is_periodic() ? sFE->polyhedron(ipoly) : sFE->polyhedron() ) ;
         ofile << poly->nb_vertices() << "  " ;
         for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
         {
            GE_Point const* pt = poly->vertex( i ) ;
            ofile << VERTICES->index( pt )+1 << "  " ;
         }
         ofile << cfe_0->mesh_id()+1 << " " << cfe_1->mesh_id()+1 << endl ;
         if( TRACE_SIDES )
         {
            OFS << "---------------------------------------" << endl ;
            OFS << "side" << endl ;
            sFE->print_current_mesh( OFS, 3 ) ;
         }
      }
   }
   
   if( TRACE_BOUNDS )
   {
      OFS << "***************************************" << endl ;
      OFS << "PDE_LocalFEbound : " << bFE->nb_meshes() << " bounds" << endl ;
      OFS << "***************************************" << endl ;
   }
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      GE_Mpolyhedron const* poly = bFE->polyhedron() ;
      ofile << poly->nb_vertices() << "  " ;
      for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
      {
         GE_Point const* pt = poly->vertex( i ) ;
         ofile << VERTICES->index( pt )+1 << "  " ;
      }
      ofile << bFE->adjacent_cell_id()+1 << " 0" << endl ;
      if( TRACE_BOUNDS )
      {
         OFS << "---------------------------------------" << endl ;
         OFS << "bound" << endl ;
         bFE->print_current_mesh( OFS, 3 ) ;
         OFS << "Local discretization of fields" << endl ;
         bFE->print_local_discretization_of_fields( OFS, 3 ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: write_gmv_faces( std::ostream& ofile,
                                  GE_SetOfPoints const* skin_vertices,
                                  size_t nb_skin_cells,
                                  size_t_vector const& glob_cell_2_skin_cell )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: write_gmv_faces" ) ;

   ofile << "faces " << bFE->nb_meshes() << " " << nb_skin_cells << endl ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      GE_Mpolyhedron const* poly = bFE->polyhedron() ;
      ofile << poly->nb_vertices() << "  " ;
      for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
      {
         GE_Point const* pt = poly->vertex( i ) ;
         ofile << skin_vertices->index( pt )+1 << "  " ;
      }
      ofile << glob_cell_2_skin_cell( bFE->adjacent_cell_id() )+1 
            << " 0" << endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: write_gmv_cells( std::ostream& ofile )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: write_gmv_cells" ) ;

   ofile << "cells " << cFE->nb_meshes() << endl ;
   if( TRACE_CELLS )
   {
      OFS << "***************************************" << endl ;
      OFS << "PDE_LocalFEcell : " << cFE->nb_meshes() << " cells" << endl ;
      OFS << "***************************************" << endl ;
   }
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      GE_Mpolyhedron const* poly = cFE->polyhedron() ;
      size_t nb_verts = poly->nb_vertices() ;
      std::string cell_type ;
      if( NB_DIMS==1 && nb_verts==2 )
      {
         cell_type = "line " ;
      }
      else if( NB_DIMS==2 && nb_verts==3 )
      {
         cell_type = "tri " ;
      }
      else if( NB_DIMS==2 && nb_verts==4 )
      {
         cell_type = "quad " ;
      }
      else if( NB_DIMS==3 && nb_verts==4 )
      {
         cell_type = "tet " ;
      }
      else if( NB_DIMS==3 && nb_verts==8 )
      {
         // Instead of standard hex we use the Patran numbered hexahedron 
         cell_type = "phex8 " ;
      }
      else 
      {
         PEL_Error::object()->raise_plain( "Invalid polyhedron type for GMV output " ) ;
      }

      ofile << cell_type << "  " << nb_verts << endl ;
      for( size_t i=0 ; i<nb_verts ; ++i )
      {
         GE_Point const* pt = poly->vertex( i ) ;
         ofile << VERTICES->index( pt )+1 << "  " ;
      }
      ofile << endl ;

      if( TRACE_CELLS )
      {
         OFS << "---------------------------------------" << endl ;
         OFS << "cell" << endl ;
         cFE->print_current_mesh( OFS, 3 ) ;
         OFS << "Local discretization of fields" << endl ;
         cFE->print_local_discretization_of_fields( OFS, 3 ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: write_gmv_surface( std::ostream& ofile,
                                    GE_SetOfPoints const* verts,
                                    bool with_sides ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: write_gmv_surface" ) ;

   if( with_sides )
   {
      ofile << "surface " << NBS + bFE->nb_meshes()<< endl ;
      for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
      {
         for( size_t ipoly=0 ; ipoly<( sFE->is_periodic() ? 2 : 1 ) ; ipoly++ )
         {
         
            GE_Mpolyhedron const* poly = ( sFE->is_periodic() ? sFE->polyhedron(ipoly) : sFE->polyhedron() ) ;
            ofile << poly->nb_vertices() << "  " ;
            for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
            {
               GE_Point const* pt = poly->vertex( i ) ;
               ofile << verts->index( pt )+1 << "  " ;
            }
            ofile << endl ;
         }
      }
   }
   else
   {
      ofile << "surface " << bFE->nb_meshes() << endl ;      
   }
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      GE_Mpolyhedron const* poly = bFE->polyhedron() ;
      ofile << poly->nb_vertices() << "  " ;
      for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
      {
         GE_Point const* pt = poly->vertex( i ) ;
         ofile << verts->index( pt )+1 << "  " ;
      }
      ofile << endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: write_gmv_surfvars( std::ostream& ofile,
                                     bool with_sides ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: write_gmv_surfvars" ) ;
   
   std::set< std::string > colors ;
   append_matching_color_names( bFE, colors ) ;
   append_matching_color_names( sFE, colors ) ;

   ofile << "surfvars" << endl ;

   std::ostringstream msg ;
   msg << "all" << colors.size() << "colors" ;
   ofile << msg.str() << endl ;
   if( with_sides )
   {
      for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
      {
         size_t nb_col = ( sFE->is_periodic() ? 2 : 1 ) ;
         for( size_t icol=0 ; icol<nb_col ; icol++ )
         {
            GE_Color const* fecol = ( sFE->is_periodic() ? 
                                         sFE->color( icol ) : sFE->color() ) ;
            ofile << fecol->identifier() << " " ;
         }
      }
   }
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      ofile << bFE->color()->identifier() << " " ;
   }
   ofile << endl ;

   std::set< std::string >::const_iterator it = colors.begin() ;
   for( ; it != colors.end() ; ++it )
   {
      std::string const& nn = *it ;
      GE_Color const* col = GE_Color::object( nn ) ;
      if( col == GE_Color::halo_color() )
      {
         ofile << HALO_COLOR << endl ;
      }
      else if( col == GE_Color::null_color() )
      {
         ofile << NULL_COLOR << endl ;         
      }
      else
      {
         ofile << nn << endl ;
      }
      size_t nb_s = 0 ;
      if( with_sides )
      {
         for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
         {
            size_t nb_col = ( sFE->is_periodic() ? 2 : 1 ) ;
            for( size_t icol=0 ; icol<nb_col ; icol++ )
            {
               GE_Color const* fecol = ( sFE->is_periodic() ? 
                                           sFE->color( icol ) : sFE->color() ) ;
               if( col->is_matching( fecol ) )
               {
                  ++nb_s ;
                  ofile << col->identifier() << " " ;
               }
               else
                  ofile << "-1 " ;
            }
         }
      }
      size_t nb_b = 0 ;
      for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
      {
         if( col->is_matching( bFE->color() ) )
         {
            ++nb_b ;
            ofile << col->identifier() << " " ;
         }
         else
            ofile << "-1 " ;
      }
      ofile << endl ;
      if( with_sides ) display_color_info( col->name(), nb_s, nb_b ) ;
   }
   ofile << "endsvars" << endl ;
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: write_gmv_variable( std::ostream& ofile ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: write_gmv_variable" ) ;

   std::set< std::string > colors ;
   append_matching_color_names( cFE, colors ) ;

   ofile << "variable" << endl ;

   std::ostringstream msg ;
   msg << "all" << colors.size() << "colors" ;
   ofile << msg.str() << " 0 " << endl ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      ofile << cFE->color()->identifier() << " " ;
   }
   ofile << endl ;

   std::set< std::string >::const_iterator it = colors.begin() ;
   for( ; it != colors.end() ; ++it )
   {
      std::string const& nn = *it ;
      GE_Color const* col = GE_Color::object( nn ) ;
      if( col == GE_Color::halo_color() )
      {
         ofile << HALO_COLOR << " 0 " << endl ;
      }
      else if( col == GE_Color::null_color() )
      {
         ofile << NULL_COLOR << " 0 " << endl ;
      }
      else
      {
         ofile << nn << " 0 " << endl ;
      }
      size_t nb_c = 0 ;
      for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
      {
         if( col->is_matching( cFE->color() ) )
         {
            ++nb_c ;
            ofile << col->identifier() << " " ;
         }
         else
            ofile << "-1 " ;
      }
      ofile << endl ;
      display_color_info( col->name(), nb_c ) ;
   }
   ofile << "endvars" << endl ;
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: append_matching_color_names( 
                                       PDE_LocalFE* fe,
                                       std::set< std::string >& colors ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: append_matching_color_names" ) ;
   
   stringVector const& coltable = GE_Color::color_table() ;
   
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      // if( fe->color() != GE_Color::null_color() )
      colors.insert( fe->color()->name() ) ;
      for( size_t j=0 ; j<coltable.size() ; ++j )
      {
         GE_Color const* col = GE_Color::object( coltable( j ) ) ;
         if( col->is_composite() && col->has( fe->color()->name() ) )
         {
            colors.insert( coltable( j ) ) ;
         }
      }
   }
}

// exactly the same as above... 
// but PDE_CursorFEside is not derived from PDE_LocalFE
//----------------------------------------------------------------------
void
PDE_DomainVisu:: append_matching_color_names( 
   PDE_CursorFEside* fe,
   std::set< std::string >& colors ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: append_matching_color_names" ) ;
   
   stringVector const& coltable = GE_Color::color_table() ;
   
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      size_t nb_col = ( fe->is_periodic() ? 2 : 1 ) ;
      for( size_t icol=0 ; icol<nb_col ; icol++ )
      {
         GE_Color const* fecol = ( fe->is_periodic() ? fe->color(icol) : fe->color() ) ;
         
         if( fecol != GE_Color::null_color() )
         {
            colors.insert( fecol->name() ) ;
         }
         for( size_t j=0 ; j<coltable.size() ; ++j )
         {
            GE_Color const* col = GE_Color::object( coltable( j ) ) ;
            if( col->is_composite() && col->has( fecol->name() ) )
            {
               colors.insert( coltable( j ) ) ;
            }
         }
      }
   }
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: display_color_info( std::string const& col_name,
                                     size_t nb_s, size_t nb_b ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: display_color_info" ) ;

   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   if( col_name.size() > 32 )
   {
      PEL::out() << endl << "*** PDE_DomainVisu: " << endl ;
      PEL::out() << "    the color \"" << col_name << "\" has " 
                 << col_name.size() << " characters and" << endl
                 << "    gmv limits the number of characters of " << endl
                 << "    surface field names to 32" << endl ;
   }

   static size_t nn = 0 ;
   if( nn == 0 )
   {
      PEL::out() << setw( 32 ) << "color"
                 << setw( 16 ) << "nb_sides"
                 << setw( 16 ) << "nb_bounds" << endl ;
      PEL::out() << setw( 32 ) << "-----"
                 << setw( 16 ) << "--------"
                 << setw( 16 ) << "---------" << endl ;
      nn = 1 ;
   }
   PEL::out() << setw( 32 ) << col_name
              << setw( 16 ) << nb_s
              << setw( 16 ) << nb_b << endl ;
   

   PEL::out().flags( original_flags ) ;
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: display_color_info( std::string const& col_name,
                                     size_t nb_c ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: display_color_info" ) ;

   ios_base::fmtflags original_flags = PEL::out().flags() ;
   PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   if( col_name.size() > 32 )
   {
      PEL::out() << endl << "*** PDE_DomainVisu: " << endl ;
      PEL::out() << "    the color \"" << col_name << "\" has " 
                 << col_name.size() << " characters and" << endl
                 << "    gmv limits the number of characters of " << endl
                 << "    variable names to 32" << endl ;
   }

   static size_t nn = 0 ;
   if( nn == 0 )
   {
      PEL::out() << endl ;
      PEL::out() << setw( 32 ) << "color"
                 << setw( 16 ) << "nb_cells" << endl ;
      PEL::out() << setw( 32 ) << "-----"
                 << setw( 16 ) << "---------" << endl ;
      nn = 1 ;
   }
   PEL::out() << setw( 32 ) << col_name
              << setw( 16 ) << nb_c << endl ;
   

   PEL::out().flags( original_flags ) ;
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: check_connexity( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: check_connexity" ) ;

   PEL::out() << endl << "*** Checking mesh connexity" << endl ;
   boolVector ok_cell( cFE->nb_meshes() ) ;
   ok_cell.set( false ) ;

   append_neigh( 0, ok_cell ) ;

   bool ok = true ;
   size_t nb_connected_cells = 0 ;
   for( size_t i=0 ; i<cFE->nb_meshes() ; ++i )
   {
      if( !ok_cell( i ) )
      {
         ok = false ;
      }
      else
      {
         ++nb_connected_cells ;
      }
   }

   if( ok )
   {
      PEL::out() << "       successful " << endl ;
   }
   else
   {
      PEL::out() << "       FAILED" << endl ;
   }
   PEL::out() << "       " << nb_connected_cells 
              << " cells traversed starting from cell 0" << endl << endl ;   
}

//----------------------------------------------------------------------
void
PDE_DomainVisu:: append_neigh( size_t icell, boolVector& ok_cell )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainVisu:: append_neigh" ) ;

   PDE_LocalFEcell const* fe_1 = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* fe_2 = sFE->adjacent_localFEcell( 1 ) ;

   cFE->go_i_th( icell ) ;
   PEL_ASSERT( cFE->is_valid() ) ;
   
   ok_cell( icell ) = true ;

   //???? probablement pas efficace pour le temps CPU
   size_t_vector sides = cFE->adjacent_side_ids() ;
   for( size_t i=0 ; i<sides.size() ; ++i )
   {
      sFE->go_i_th( sides( i ) ) ;
      size_t icell_1 = fe_1->mesh_id() ;
      size_t icell_2 = fe_2->mesh_id() ;
      size_t icell_other = PEL::bad_index() ;
      if( icell_1 == icell ) 
      {
         icell_other = icell_2 ;
      }
      else if( icell_2 == icell )
      {
         icell_other = icell_1 ;
      }
      PEL_ASSERT( icell_other != PEL::bad_index() ) ;
      if( !ok_cell( icell_other ) )
         append_neigh( icell_other, ok_cell ) ;
   }
}

//internal--------------------------------------------------------------
void
PDE_DomainVisu_ERROR:: n0( std::string const& file_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl
       << "*** PDE_DomainVisu: saving failure" << std::endl << std::endl ;
   msg << "    Unable to open file \""
       << file_name << "\" for writing" << std::endl ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
