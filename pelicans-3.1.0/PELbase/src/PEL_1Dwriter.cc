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

#include <PEL_1Dwriter.hh>

#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Map.hh>
#include <PEL_String.hh>

#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::endl ; using std::cout ;
using std::ostringstream ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

PEL_1Dwriter const* PEL_1Dwriter::PROTOTYPE = new PEL_1Dwriter() ;

//----------------------------------------------------------------------
PEL_1Dwriter:: PEL_1Dwriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "PEL_1Dwriter" )
   , CELL_VARS_NAMES( 0 )
   , VERT_VARS_NAMES( 0 )
   , CELL_MATRIX( 0, 0 )
   , VERT_MATRIX( 0, 0 )
{
}

//----------------------------------------------------------------------
PEL_1Dwriter*
PEL_1Dwriter:: create_replica( PEL_Object* a_owner,
				PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_1Dwriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_1Dwriter* result = new PEL_1Dwriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PEL_1Dwriter:: PEL_1Dwriter( PEL_Object* a_owner,
			       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , HAS_CELL_VARS( false )
   , CELL_VARS_NAMES( 0 )
   , HAS_VERT_VARS( false )
   , VERT_VARS_NAMES( 0 )
   , DIM( PEL::bad_index() )
   , CELL_MATRIX( 0, 0 )
   , VERT_MATRIX( 0, 0 )
{
   FILEBASENAME = exp->string_data( "files_basename" ) ;
   FILEXTENSION = ".1d" ;
   OPENMODE =  std::ios::out ;

   FORMAT = "ascii" ;

   if( exp->string_data( "writing_mode" ) == "binary" )
   {
      ostringstream mesg ;
      mesg << "*** PEL_1Dwriter :" << endl ;
      mesg << "    the only \"writing_mode\" that is handled is : \"text\"" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}


//----------------------------------------------------------------------
PEL_1Dwriter:: ~PEL_1Dwriter( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_1Dwriter:: determine_conditions( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_1Dwriter:: determine_conditions" ) ;

   if( DIM == PEL::bad_index() )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "meshing" ) ;
      DIM = se->int_data( "nb_sp_dims" ) ;
      se->destroy() ; se=0 ;

      if( DIM != 1 )
      {
         ostringstream mesg ;
         mesg << "*** PEL_1Dwriter :" << endl ;
         mesg << "    the geometrical domain should be 1D" << endl ;
         mesg << "    (it is currently " << DIM << "D )" ;
         PEL_Error::object()->raise_plain( mesg.str() ) ;
      }

      if( exp->has_module( "fields" ) )
      {
         se = exp->create_subexplorer( 0, "fields" ) ;
         se->start_module_iterator() ;
         for( ; se->is_valid_module() ; se->go_next_module() )
         {
            PEL_ModuleExplorer const* sse = se->create_subexplorer( 0 ) ;
            std::string const& name = sse->string_data( "name" ) ;
            std::string const& location = sse->string_data( "location" ) ;
            cout << "sse->name()=" << sse->name() << endl ;
            doubleArray2D const& vals = sse->doubleArray2D_data( "value" )  ;
            cout << "lu" << endl ;
            if( location=="at_cell_centers" ) 
            {
               HAS_CELL_VARS = true ;
               size_t nb_comps = vals.index_bound( 0 ) ;
               for( size_t ic=0 ; ic<nb_comps ; ++ic )
               {
                  ostringstream nn ;
                  nn << name ;
                  if( nb_comps != 1 ) nn << "(" << ic << ")" ;
                  CELL_VARS_NAMES.append( nn.str() ) ;
               }
            }
            else if( location=="at_vertices" ) 
            {
               HAS_VERT_VARS = true ;
               size_t nb_comps = vals.index_bound( 0 ) ;
               for( size_t ic=0 ; ic<nb_comps ; ++ic )
               {
                  ostringstream nn ;
                  nn << name ;
                  if( nb_comps != 1 ) nn << "(" << ic << ")" ;
                  VERT_VARS_NAMES.append( nn.str() ) ;
               }
            }
            else raise_field_location_error(
                           name, location, 
                           "at_vertices,at_cell_centers" ) ;
            sse->destroy() ; sse=0 ;
         }
         se->destroy() ; se=0 ;
      }
   }
   else
   {
      // *** peut etre à faire : vérifs de consistance.
   }
}

//----------------------------------------------------------------------
void
PEL_1Dwriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_1Dwriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;

   determine_conditions( exp ) ;

   // Current time
   size_t icycle = exp->int_data( "cycle_number" ) ;

   // Open file and test
   if( HAS_CELL_VARS )
   {
      std::string fname = output_file_name( icycle, "_c" ) ;
      FILE_FOR_CELL_VARS.open( fname.c_str(), OPENMODE ) ;
      if( FILE_FOR_CELL_VARS.fail() || !FILE_FOR_CELL_VARS.is_open() )
      {
         PEL_Error::object()->raise_plain(
            "unable to create the 1D output file : " + fname ) ;
      }
   }
   if( HAS_VERT_VARS )
   {
      std::string fname = output_file_name( icycle, "_v" ) ;
      FILE_FOR_VERT_VARS.open( fname.c_str(), OPENMODE ) ;
      if( FILE_FOR_VERT_VARS.fail() || !FILE_FOR_VERT_VARS.is_open() )
      {
         PEL_Error::object()->raise_plain(
            "unable to create the GMV output file : " + fname ) ;
      }
   }

   // Does a new grid must be saved?
   if( exp->has_module( "meshing" ) )
   {
      PEL_ModuleExplorer const* se = exp->create_subexplorer( 0, "meshing" ) ;
      doubleArray2D const& vertices = se->doubleArray2D_data( "vertices" ) ;
      intVector const& cell_nb_verts = 
                                     se->intVector_data( "cell_nb_vertices" ) ;
      intArray2D const& cell2vertex = se->intArray2D_data( "cell2vertex" ) ;
      if( HAS_CELL_VARS )
      {
         resize_and_add_cell_coords( CELL_MATRIX, vertices, cell_nb_verts,
                                     cell2vertex, CELL_VARS_NAMES.size()+1 ) ;
      }
      if( HAS_VERT_VARS )
      {
         resize_and_add_vert_coords( VERT_MATRIX, vertices, 
                                     VERT_VARS_NAMES.size()+1 ) ;
      }
      se->destroy() ; se = 0 ;
   }

   TIME = PEL::bad_double() ;
   if( exp->has_module( "variables" ) )
   {
      PEL_ModuleExplorer const* se = exp->create_subexplorer( 0, "variables" );
      if( se->has_entry( "TIME" ) )
      {
         TIME = se->double_data( "TIME" ) ;
      }
      se->destroy() ; se=0 ;
   }

   // Saving of fields
   size_t current_vert_column = 1 ;
   size_t current_cell_column = 1 ;
   if( exp->has_module( "fields" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "fields" ) ;
      se->start_module_iterator() ;
      for( ; se->is_valid_module() ; se->go_next_module() )
      {
         PEL_ModuleExplorer const* sse = se->create_subexplorer( 0 ) ;
         doubleArray2D const& vals = sse->doubleArray2D_data( "value" ) ;
         std::string const& location = sse->string_data( "location" ) ;
         if( location == "at_cell_centers" )
         {
            add_one_field( CELL_MATRIX, current_vert_column, vals ) ;
         }
         else if( location=="at_vertices" )
         {
            add_one_field( VERT_MATRIX, current_cell_column, vals ) ;
         }

         sse->destroy() ; sse=0 ;
      }
      se->destroy() ; se = 0 ;
   }

   if( HAS_CELL_VARS )
   {
      write_matrix( CELL_MATRIX, TIME, CELL_VARS_NAMES, FILE_FOR_CELL_VARS ) ;
   }
   if( HAS_VERT_VARS )
   {
      write_matrix( VERT_MATRIX, TIME, VERT_VARS_NAMES, FILE_FOR_VERT_VARS ) ;
   }

   if( HAS_CELL_VARS ) FILE_FOR_CELL_VARS.close() ;
   if( HAS_VERT_VARS ) FILE_FOR_VERT_VARS.close() ;
}

//----------------------------------------------------------------------
void
PEL_1Dwriter:: resize_and_add_cell_coords( doubleArray2D& cell_matrix,
                                           doubleArray2D const& vertices,
                                           intVector const& cell_nb_verts,
                                           intArray2D const& cell2vert,
                                           size_t nb_columns )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_1Dwriter:: resize_and_add_cell_coords" ) ;
   PEL_ASSERT( vertices.index_bound( 0 ) == 1 ) ;

   size_t nb_cells = cell2vert.index_bound( 1 ) ;
   cell_matrix.re_initialize( nb_cells, nb_columns ) ;
   for( size_t i=0 ; i<nb_cells ; ++i )
   {
      size_t ii = (size_t)i ;
      PEL_ASSERT( cell_nb_verts( ii ) == 2 ) ;
      cell_matrix( ii, 0 ) = ( vertices( 0, cell2vert( 0, ii ) ) + 
                               vertices( 0, cell2vert( 1, ii ) ) ) / 2.0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_1Dwriter:: resize_and_add_vert_coords( doubleArray2D& vert_matrix,
                                           doubleArray2D const& vertices,
                                           size_t nb_columns )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_1Dwriter:: resize_and_add_vert_coords" ) ;
   PEL_ASSERT( vertices.index_bound( 0 ) == 1 ) ;

   size_t nb_verts = vertices.index_bound( 1 ) ;
   vert_matrix.re_initialize( nb_verts, nb_columns ) ;
   for( size_t i=0 ; i<nb_verts ; ++i )
   {
      vert_matrix( i, 0 ) = vertices( 0, i ) ;
   }
}

//-----------------------------------------------------------------------
void
PEL_1Dwriter:: add_one_field( doubleArray2D& matrix,
                              size_t& current_column,
                              doubleArray2D const& vals )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_1Dwriter:: add_one_field" ) ;

   size_t nb_comps = vals.index_bound( 0 ) ;
   size_t nb_lines = vals.index_bound( 1 ) ;
   PEL_ASSERT( matrix.index_bound( 0 ) == nb_lines  ) ;

   for( size_t ic=0 ; ic<nb_comps ; ++ic )
   {
      for( size_t i=0 ; i<nb_lines ; ++i )
      {
         matrix( i, current_column ) = vals( ic, i ) ;
      }
      ++current_column ;
   }
}

//-----------------------------------------------------------------------
void
PEL_1Dwriter:: write_matrix( doubleArray2D const& matrix,
                             double time,
                             stringVector const& var_names,
                             std::ofstream& os )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_1Dwriter:: write_matrix" ) ;
   PEL_ASSERT( var_names.size()+1 == matrix.index_bound( 1 ) ) ;

   ios_base::fmtflags original_flags = os.flags() ;
   os.setf( ios_base::uppercase | ios_base::scientific ) ;

   os << "#" << endl ;
   if( time != PEL::bad_double() )
   {
      os << "#   time : " << time << endl ;
      os << "#" << endl ;
   }
   os << "#" << setw( 19 ) << "coordinate" ;
   size_t nb_vars = var_names.size() ;
   for( size_t i=0 ; i<nb_vars ; ++i )
   {
      os << setw( 20 ) << var_names( i ) ;
   }
   os << endl << "#" << endl ;
   size_t nb_lines = matrix.index_bound( 0 ) ;
   size_t nb_cols = matrix.index_bound( 1 ) ;
   for( size_t i=0 ; i<nb_lines ; ++i )
   {
      for( size_t j=0 ; j<nb_cols ; ++j )
      {
         os << setprecision( 10 ) << setw( 20 ) <<  matrix( i, j ) ;
      }
      os << endl ;
   }

   os.flags( original_flags ) ;
}

//-----------------------------------------------------------------------
std::string
PEL_1Dwriter:: output_file_name( size_t nb, std::string add_string )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_1Dwriter:: output_file_name" ) ;
   PEL_CHECK( nb<99999 ) ;

   std::ostringstream tmp ;
   tmp << nb ;
   std::string nb_string = tmp.str() ;
   std::string result = FILEBASENAME+add_string+FILEXTENSION+".00000";
   result.replace( result.length()-nb_string.length(), 
                   nb_string.length(),
                   nb_string ) ;

   PEL_CHECK( result.length()==( FILEBASENAME.length()+add_string.length()
                                +FILEXTENSION.length()+6) ) ;

   return( result ) ;
}
