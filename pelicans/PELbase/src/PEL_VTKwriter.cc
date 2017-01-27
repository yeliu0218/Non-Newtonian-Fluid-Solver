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

#include <PEL_VTKwriter.hh>

#include <PEL_Communicator.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Int.hh>
#include <PEL_IntVector.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Map.hh>
#include <PEL_Module.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>

#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifdef ZLIB
#include <zlib.h>
#endif

using std::endl ; 

PEL_VTKwriter const* PEL_VTKwriter::PROTOTYPE = new PEL_VTKwriter() ;
std::string PEL_VTKwriter::UNAMED = "_unamed_" ;

static size_t sizeof_Float32 = 4 ;
static size_t sizeof_Int32 = 4 ;

//----------------------------------------------------------------------
bool
PEL_VTKwriter:: is_parallel_writer( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_VTKwriter:: PEL_VTKwriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "PEL_VTKwriter" )
   , SPACE_DIM( PEL::bad_index() )
   , PVD_STRINGS( 0 )
   , BINARY( false )
   , CYCLE_NUMBER( 0 )
   , VERTICES( 0, 0 )
   , LINEAR_CONNECTIVITY( 0 ) 
   , OFFSET_CONNECTIVITY( 0 )
   , CELL_TYPES( 0 )
{
}

//----------------------------------------------------------------------
PEL_VTKwriter*
PEL_VTKwriter:: create_replica( PEL_Object* a_owner,
				PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_VTKwriter* result = new PEL_VTKwriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PEL_VTKwriter:: PEL_VTKwriter( PEL_Object* a_owner,
			       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , SPACE_DIM( PEL::bad_index() )
   , BASE_FILENAME( exp->string_data( "files_basename" ) )
   , PVD_FILENAME( exp->string_data( "files_basename" ) + ".pvd" )
   , PVD_STRINGS( 0 )
   , BINARY( false )
   , CYCLE_NUMBER( 0 )
   , VERTICES( 0, 0 )
   , LINEAR_CONNECTIVITY( 0 )
   , OFFSET_CONNECTIVITY( 0 )
   , CELL_TYPES( 0 )
   , BUFFER( 0 )
   , ALLOCATED( 0 )
   , OFFSET( 0 )
   , ENCODING( PEL_System::big_endian_encoding() ?
                                       "BigEndian" : "LittleEndian" ) 
   , COMPRESS( false )
   , COM( PEL_Exec::communicator() )
{
   if( exp->string_data( "writing_mode" )=="binary" )
   {
      BINARY = true ;
      if( exp->has_entry( "compress" ) )
      {
         COMPRESS = exp->bool_data( "compress" ) ;
      }
#ifdef ZLIB
      else COMPRESS=true ;
#else
      if( COMPRESS ) 
      {
		  PEL_Error::object()->raise_plain(
           "No compression allowed without ZLIB library linking" ) ;
      }
#endif
   }
   PEL_CHECK_INV( invariant() ) ;
}


//----------------------------------------------------------------------
PEL_VTKwriter:: ~PEL_VTKwriter( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}


//----------------------------------------------------------------------
void
PEL_VTKwriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;

   // Current time
   CYCLE_NUMBER = exp->int_data( "cycle_number" ) ;

   // Open files and test
   std::string fname = output_file_name( CYCLE_NUMBER, false, COM->rank() ) ;
   std::ofstream file( fname.c_str() ) ;
   if( !file )
   {
      PEL_Error::object()->raise_plain(
	 "unable to open the VTK output file : " + fname ) ;
   }
   if( COM->nb_ranks()>1 && COM->rank()==0 )
   {
      write_pvd_file( exp, output_file_name( CYCLE_NUMBER, true, 0 ) ) ;
   }
   else if( COM->rank()==0 )
   {
      write_pvd_file( exp, output_file_name( CYCLE_NUMBER, false, 0 ) ) ;
   }
   PEL_Module * vtk = PEL_Module::create( 0, "VTKFile" ) ;
   
   if( exp->has_module( "meshing" ) )
   {
      PEL_ModuleExplorer const* se = exp->create_subexplorer( 0,
                                                              "meshing" ) ;
      store_meshing( se ) ;
      se->destroy() ;
   }
                                 
   build_vtu( vtk, exp, false ) ;
   PEL_ModuleExplorer* mexp = PEL_ModuleExplorer::create( vtk, vtk ) ;
   
   write_vtk( mexp, file, 0, false ) ;
                                 
   file.close() ;
   vtk->destroy() ;

   if( COM->nb_ranks() > 1 && COM->rank()==0 ) 
   {
      vtk = PEL_Module::create( 0, "VTKFile" ) ;
      build_vtu( vtk, exp, true ) ;
      PEL_Module* grid = vtk->module( "PUnstructuredGrid" ) ;
      grid->add_entry( "GhostLevel", PEL_Int::create( grid, 0 ) ) ;
      for( size_t i=0 ; i<COM->nb_ranks() ; i++ ) 
      {
         std::ostringstream npiece ;
         npiece << "Piece#" << i ;
         PEL_Module* piece = PEL_Module::create( grid, npiece.str() ) ;
         std::string file_name =
            PEL_System::basename( output_file_name(CYCLE_NUMBER, false, i ) ) ;
         piece->add_entry( "Source", PEL_String::create( piece, file_name ) ) ;
         grid->add_module( piece ) ;
      }
      
      std::string pfname = output_file_name( CYCLE_NUMBER, true, 0 ) ;
      
      std::ofstream pfile( pfname.c_str() ) ;
      if( !pfile )
      {
         PEL_Error::object()->raise_plain(
            "unable to open the VTK output file : " + pfname ) ;
      }
      mexp = PEL_ModuleExplorer::create( vtk, vtk ) ;
   
      write_vtk( mexp, pfile, 0, true ) ;
      pfile.close() ;
      vtk->destroy() ;
   }
}

//----------------------------------------------------------------------
void
PEL_VTKwriter:: store_meshing( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL("PEL_VTKwriter:: store_meshing") ;
   PEL_CHECK( exp!=0 ) ;

   // VERTICES
   VERTICES = exp->doubleArray2D_data( "vertices" ) ;
   
   // Cells
   intArray2D const& connectivity = exp->intArray2D_data( "cell2vertex" ) ;
   intVector const& cell_nb_vertices = exp->intVector_data( "cell_nb_vertices" ) ;
   
   NB_CELLS = connectivity.index_bound(1) ;
   LINEAR_CONNECTIVITY.resize( 0 ) ;
   OFFSET_CONNECTIVITY.resize( NB_CELLS ) ;
   CELL_TYPES.resize( NB_CELLS ) ;
   SPACE_DIM = VERTICES.index_bound( 0 ) ;
   
   size_t n = 0 ;
   for( size_t i=0 ; i<NB_CELLS ; ++i )
   {
      size_t const vert_by_mesh = cell_nb_vertices( i ) ;
      VTK_TYPE elm_type = VTK_UNDEF ;
      if( SPACE_DIM==1 && vert_by_mesh==2 )
      {
         elm_type = VTK_LINE ;
      }
      else if( SPACE_DIM==2 && vert_by_mesh==3 )
      {
         elm_type = VTK_TRIANGLE ;
      }
      else if( SPACE_DIM==2 && vert_by_mesh==4 )
      {
         elm_type = VTK_QUAD ;
      }
      else if( SPACE_DIM==3 && vert_by_mesh==4 )
      {
         elm_type = VTK_TETRA ;
      }
      else if( SPACE_DIM==3 && vert_by_mesh==8 )
      {
         elm_type = VTK_HEXAHEDRON ;
      }
      else 
      {
         PEL_Error::object()->raise_plain(
            "Invalid polyhedron type for VTK output " ) ;
      }
      for( size_t j=0 ; j<vert_by_mesh ; j++ )
      {
         LINEAR_CONNECTIVITY.append( connectivity(j,i) ) ;
      }
      n += vert_by_mesh ;
      OFFSET_CONNECTIVITY(i) = n ;
      CELL_TYPES(i) = (int) elm_type ;
   }   
}

//-----------------------------------------------------------------------
void
PEL_VTKwriter:: write_pvd_file( PEL_ModuleExplorer const* exp,
                                std::string const& vtu_filename )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: write_pvd_file" ) ;

   // First step : add a new line for new cycle if any
   std::string new_line = "<DataSet timestep=\"" ;
   std::ostringstream oss;
   //oss << exp->int_data( "cycle_number" ) ;
   oss << exp->double_data( "Time" ) ;

   new_line += oss.str() ;
   new_line += "\" group=\"\" part=\"0\" file=\"" ;
   new_line += PEL_System::basename( vtu_filename ).c_str() ;
   new_line += "\"/>" ;

   PVD_STRINGS.append( new_line ) ;
   
   // Second step : build pvd file with header, one line by cycle and end. 
   std::ofstream file( PVD_FILENAME.c_str(), std::ios::trunc ) ;
   if( !file )
   {
      PEL_Error::object()->raise_plain(
	 "unable to open the FTK output file : " + PVD_FILENAME ) ;
   }

   file << "<?xml version=\"1.0\"?>" << endl ;
   file << "<VTKFile type=\"Collection\" version=\"0.1\""
        << " byte_order=\"" << ENCODING << "\"" ;
   if( COMPRESS )
   {
      file << " compressor=\"vtkZLibDataCompressor\"" ;
   }
   file << ">" << endl ;
   file << "<Collection>" << endl ;

   for( size_t i=0; i<PVD_STRINGS.size(); ++i )
   {
      file << PVD_STRINGS( i ) << endl ;
   }

   file << "</Collection>" << endl ;
   file << "</VTKFile>" << endl ;
   file.close() ;
}

//----------------------------------------------------------------------
void
PEL_VTKwriter:: build_vtu( PEL_Module* vtk,
                           PEL_ModuleExplorer const* exp,
                           bool parallel ) 
//----------------------------------------------------------------------
{
   PEL_LABEL("PEL_VTKwriter:: build_vtu") ;
   PEL_CHECK( exp!=0 ) ;

   std::string pre = ( parallel ? "P" : "" ) ;
   
   vtk->add_entry( "type", PEL_String::create( vtk,
                                               pre+"UnstructuredGrid" ) ) ;

   if( BINARY )
   {
      vtk->add_entry( "byte_order",
		      PEL_String::create( vtk, ENCODING ) ) ;
   }
   
   if( COMPRESS )
      vtk->add_entry( "compressor", PEL_String::create( vtk, "vtkZLibDataCompressor" ) ) ;
   
   PEL_Module* grid = PEL_Module::create( vtk, pre+"UnstructuredGrid" ) ;
   PEL_Module* piece = PEL_Module::create( grid, "Piece" ) ;
   PEL_Module* base = ( parallel ? grid : piece ) ;
   
   write_grid(base, parallel) ;
   
   if( exp->has_module( "fields" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "fields" ) ;
      write_fields( base, se,  parallel) ;
      se->destroy() ; se = 0 ;
      
   }
      
   if( !parallel ) grid->add_module(piece) ;
   vtk->add_module(grid) ;
   
}

//----------------------------------------------------------------------
void
PEL_VTKwriter:: write_grid( PEL_Module* base,
                            bool parallel ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: write_grid" ) ;
   std::string pre = ( parallel ? "P" : "" ) ;

   if( !parallel ) 
   {  
      base->add_entry( "NumberOfPoints",
                       PEL_Int::create( base,
                                        VERTICES.index_bound(1) ) ) ;
      
      base->add_entry( "NumberOfCells",
                       PEL_Int::create( base, NB_CELLS ) ) ;
   }
   
   PEL_Module* Points = PEL_Module::create( base, pre+"Points" ) ;
   Points->add_entry( UNAMED, PEL_DoubleArray2D::create(Points,VERTICES) ) ;
   
   PEL_Module* PointData = PEL_Module::create( base, pre+"PointData" ) ;
   
   PEL_Module* Cells = PEL_Module::create( base, pre+"Cells" ) ;
   Cells->add_entry( "connectivity",
                     PEL_IntVector::create(Cells,LINEAR_CONNECTIVITY) ) ;
   Cells->add_entry( "offsets",
                     PEL_IntVector::create(Cells,OFFSET_CONNECTIVITY) ) ;
   Cells->add_entry( "types",
                     PEL_IntVector::create(Cells,CELL_TYPES) ) ;

   PEL_Module* CellData = PEL_Module::create( base, pre+"CellData" ) ;
   
   base->add_module(Cells) ;
   base->add_module(Points) ;
   base->add_module(PointData) ;
   base->add_module( CellData ) ;
}


//----------------------------------------------------------------------
void
PEL_VTKwriter:: write_fields( PEL_Module* base,
                              PEL_ModuleExplorer* exp,
                              bool parallel ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: write_fields" ) ;
   PEL_CHECK( exp!=0 ) ;

   std::string pre = ( parallel ? "P" : "" ) ;
   
   PEL_Module*point_data = base->module(pre+"PointData") ;
   PEL_Module*cell_data = base->module(pre+"CellData") ;
   for( exp->start_module_iterator()  ;
        exp->is_valid_module() ;
        exp->go_next_module() )
   {
      PEL_ModuleExplorer const* fexp = exp->create_subexplorer( 0 ) ;
      write_one_field( fexp, point_data, cell_data ) ;
      fexp->destroy() ; fexp = 0 ;
   }
   
}


//----------------------------------------------------------------------
void
PEL_VTKwriter:: write_one_field( PEL_ModuleExplorer const* fexp,
                                 PEL_Module* point_data,
                                 PEL_Module* cell_data ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: write_one_field" ) ;
   PEL_CHECK( fexp!=0 ) ;
   PEL_CHECK( point_data!=0 ) ;
   PEL_CHECK( cell_data!=0 ) ;

   std::string const& name = fexp->string_data( "name" ) ;
   doubleArray2D const& X = fexp->doubleArray2D_data( "value" ) ;
   std::string const& location = fexp->string_data( "location" ) ;
   PEL_Module* target = cell_data ;
   
   if( location == "at_vertices" ) 
   {
      target = point_data ;
   }
   else if( location != "at_cell_centers" )
   {
      raise_field_location_error(
         name, location, "at_vertices,at_cell_centers" ) ;
   }

   target->add_entry( name, PEL_DoubleArray2D::create( target, X ) ) ;
   bool scalar = X.index_bound(0)==1 ;
   if(scalar&&!target->has_entry("Scalars")) 
   {
      target->add_entry( "Scalars", PEL_String::create( target, name ) ) ;
   }
   else if( !scalar&&!target->has_entry("Vectors")) 
   {
      target->add_entry( "Vectors", PEL_String::create( target, name ) ) ;
   }
}

//-----------------------------------------------------------------------
std::string
PEL_VTKwriter:: output_file_name( size_t nb,
                                  bool parallel,
                                  size_t rank )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: output_file_name" ) ;
   PEL_CHECK( nb<99999 ) ;
   std::string file_extension = ( parallel ? ".pvtu" : ".vtu" ) ;
   
   std::ostringstream tmp ;
   tmp << BASE_FILENAME ;
   tmp << "T" << nb ;
   if( !parallel && COM->nb_ranks() > 1 )
   {
      tmp << "_" << rank ;
   }
   
   std::string result = tmp.str() + file_extension ;

   return( result ) ;
}


//-----------------------------------------------------------------------
void
PEL_VTKwriter:: write_vtk( PEL_ModuleExplorer* vtk,
                           std::ofstream& file,
                           size_t level,
                           bool parallel )  
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: write_vtk" ) ;

   std::string data_array = ( parallel ? "PDataArray" : "DataArray" ) ;
   
   if( level==0 )
      file << "<?xml version=\"1.0\"?>" << endl ;
   std::string bl( 2*level, ' ' ) ;
   bl = "\n" + bl ;
   std::string name =  vtk->name() ;
   if( name.find( "#" )<name.length() )
      name = name.substr( 0, name.find( "#" ) ) ;
   
   file << bl << "<" << name ;
   for( vtk->start_entry_iterator() ;
        vtk->is_valid_entry() ;
        vtk->go_next_entry() ) 
   {
      PEL_Data* data = vtk->data( 0 ) ;
      PEL_Data::Type dt = data->data_type() ;
      
      if( dt==PEL_Data::String  ) 
      {
         file << " " << vtk->keyword() << "=" ;
         data->print(file,0) ;
         file<<" " ;
      }
      else if(  dt==PEL_Data::Int ) 
      {
         file << " " << vtk->keyword() << "=\"" ;
         data->print(file,0) ;
         file<<"\" " ;
      }
      data->destroy() ;
   }
   file << ">" ;
   for( vtk->start_entry_iterator() ;
        vtk->is_valid_entry() ;
        vtk->go_next_entry() ) 
   {
      PEL_Data* data = vtk->data( 0 ) ;
      PEL_Data::Type dt = data->data_type() ;
      if( dt==PEL_Data::DoubleArray2D || dt==PEL_Data::IntVector ) 
      {
         file << "<" << data_array << " " ;
         
         if( vtk->keyword()!=UNAMED ) 
         {
            file << "Name=\"" << vtk->keyword() << "\" " ;
         }
         std::string a_type = "Int32" ;
         
         std::string format = "ascii" ;
         
         if( BINARY && !parallel )
         {
            format="appended" ;
            file << "offset=\""<<OFFSET<<"\" " ;
         }
         
         if( dt==PEL_Data::DoubleArray2D )
         {
            int nbc = data->to_double_array2D().index_bound(0) ;
            if( nbc>1 )
               file << "NumberOfComponents=\"3\" " ;
            a_type="Float32" ;
         }
         
         file << "type=\"" << a_type << "\" " ;
         
         file << "format=\"" << format << "\" " ;

         file << ">" << endl ;

         if( !parallel )
         {
            
            if( dt==PEL_Data::DoubleArray2D ) 
            {
               doubleArray2D const& tab = data->to_double_array2D() ;
               size_t nbc = data->to_double_array2D().index_bound(0) ;
               size_t ncmax = ( nbc==1 ? 1 : 3 ) ;
            
               start_output( sizeof_Float32, ncmax*tab.index_bound(1) ) ;
               for( size_t i=0 ; i<tab.index_bound(1) ; i++ )
               {
                  for( size_t j=0 ; j<ncmax ; j++ )
                     write_double( file, ( j<nbc ? tab(j,i) : 0.0 ) ) ;
                  
               }
               flush(file) ;
               
            }
            else if( dt==PEL_Data::IntVector  )
            {
               intVector const& tab = data->to_int_vector() ;
               start_output( sizeof_Int32, tab.size() ) ;
               for( size_t i=0 ; i<tab.size() ; i++ )
                  write_int( file, tab(i) ) ;
               flush(file) ;
            }
         }
         
         file << "</" << data_array << ">" << endl ;
      }
      else if( ! ( dt==PEL_Data::String || dt==PEL_Data::Int ) ) 
      {
         PEL_Error::object()->raise_internal(
            " Bad type "+data->data_type() ) ;
      }
      data->destroy() ;
   }
   for( vtk->start_module_iterator() ;
        vtk->is_valid_module() ;
        vtk->go_next_module() ) 
   {
      PEL_ModuleExplorer* sexp = vtk->create_subexplorer( 0 ) ;
      write_vtk(sexp,file,level+1,parallel) ;
      sexp->destroy() ;
   }
   if( level==0 && BINARY && !parallel ) 
   {
      file << bl << "<AppendedData encoding=\"raw\" > " << endl << "    _" ;
      file.write( BUFFER, OFFSET ) ;
      file << bl << "</AppendedData>" << endl ;
      delete [] BUFFER ; BUFFER = 0 ;
      ALLOCATED = 0 ;
      OFFSET = 0 ;
   }
   
   file << bl << "</" << name <<">" ;
  
   
}



//-----------------------------------------------------------------------
void
PEL_VTKwriter:: start_output( size_t size, size_t number )  
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: start_output" ) ;
   if(BINARY)
   {
      int current_output_size = size*number ;
      unsigned long ncomp = current_output_size + (current_output_size+999)/1000 + 12 + sizeof_Int32 ;
      check_allocated( ncomp ) ;
      CURRENT_LENGTH = store_int(current_output_size) ;
   }
}


//-----------------------------------------------------------------------
void
PEL_VTKwriter:: write_double( std::ofstream& file, double val )  
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: write_double" ) ;
   PEL_CHECK( sizeof_Float32==sizeof(float) ) ;
   
   if(BINARY)
   {
      PEL_CHECK( OFFSET + sizeof_Float32 <= ALLOCATED ) ;
      *((float*)&(BUFFER[OFFSET])) = (float)val ;
      OFFSET += sizeof_Float32  ;
   }
   else
      file << " " << (float)val ;
}


//-----------------------------------------------------------------------
void
PEL_VTKwriter:: write_int( std::ofstream& file, int val )  
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: write_int" ) ;
   PEL_CHECK( sizeof_Int32==sizeof(int) ) ;
   
   if(BINARY)
   {
      store_int(val) ;
   }
   else
      file << " " << val ;
}



//-----------------------------------------------------------------------
size_t
PEL_VTKwriter:: store_int( int val )  
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: store_int" ) ;
   PEL_CHECK( OFFSET + sizeof_Int32 <= ALLOCATED ) ;
   
   size_t result = OFFSET ;
   *((int*)&(BUFFER[OFFSET])) = val ;
   OFFSET += sizeof_Int32  ;
   return result ;
}


//-----------------------------------------------------------------------
void
PEL_VTKwriter:: check_allocated( size_t size )  
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: check_allocated" ) ;
   if(OFFSET+size>=ALLOCATED) 
   {
      size_t new_size = PEL::max( 2*ALLOCATED, (size_t)1024 ) ;
      new_size = PEL::max( new_size, 2*(OFFSET+size) ) ;
      new_size = 4 * ( new_size/4 +1 ) ; // allignement sur 4 bytes
      
      char * new_buffer = new char [ new_size ] ;
      for( size_t i=0 ;i<OFFSET ;i++ )
         new_buffer[i] = BUFFER[i] ;
      if( BUFFER!=0 ) delete [] BUFFER ;
      BUFFER = new_buffer ;
      ALLOCATED = new_size ;
      
   }
   PEL_CHECK_POST( OFFSET+size<ALLOCATED ) ;
}

//-----------------------------------------------------------------------
void
PEL_VTKwriter:: flush( std::ofstream& file )  
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: flush" ) ;
   if(COMPRESS)
   {
      compress_segment(CURRENT_LENGTH) ;         
   }
   file << endl ;
}



//-----------------------------------------------------------------------
void
PEL_VTKwriter:: compress_segment( size_t seg )  
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VTKwriter:: compress_segment" ) ;
   static size_t BlockSize = 32768 ;
   size_t size = (size_t)(*((int*)&BUFFER[seg])) ;
   
   size_t numFullBlocks = size / BlockSize;
   size_t lastBlockSize = size % BlockSize;
   size_t numBlocks = numFullBlocks + (lastBlockSize?1:0);

   size_t headerLength = numBlocks+3;

   int * CompressionHeader = new int[headerLength];
   CompressionHeader[0] = numBlocks;
   CompressionHeader[1] = BlockSize;
   CompressionHeader[2] = lastBlockSize;

   unsigned long encoded_buff_size = PEL::max(BlockSize,size)  ;
   unsigned char* encoded_buff = new unsigned char [ encoded_buff_size ] ;
   size_t encoded_offset = 0 ;
   for( size_t block=0 ; block<numBlocks ; block++ )
   {
      size_t buffer_start = seg + sizeof_Int32 + block*BlockSize ;
      size_t length = ( block+1<numBlocks || !lastBlockSize ? BlockSize : lastBlockSize ) ;
      unsigned char* to_encode = (unsigned char *)(&BUFFER[buffer_start]) ;
      unsigned char* encoded = &encoded_buff[encoded_offset] ;
      unsigned long ncomp = encoded_buff_size - encoded_offset ;
#ifdef ZLIB
      if(compress2((Bytef*)encoded,
                   &ncomp,
                   (const Bytef*)to_encode,
                   length,
                   Z_DEFAULT_COMPRESSION) != Z_OK)
#endif
      {
         PEL_Error::object()->raise_plain(
            "Zlib error while compressing data.");
      }
      CompressionHeader[3+block] = ncomp ;
      encoded_offset += ncomp ;
   }
   
   OFFSET = seg ;
   check_allocated( headerLength*sizeof_Int32 + encoded_offset ) ;
   
   for(size_t i=0 ; i<headerLength ; i++ )
      store_int(CompressionHeader[i]) ;     

   for(size_t i=0 ; i<encoded_offset ; i++ )
      BUFFER[OFFSET++] = encoded_buff[i] ;

   if( OFFSET%4 != 0 )
      OFFSET = 4*( OFFSET/4 +1 ) ; // Re-allignement
   
   delete [] CompressionHeader ;
   delete [] encoded_buff ;
}

