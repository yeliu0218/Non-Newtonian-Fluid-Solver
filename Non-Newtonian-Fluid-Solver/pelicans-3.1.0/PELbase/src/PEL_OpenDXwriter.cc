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

#include <PEL_OpenDXwriter.hh>

#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Map.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>

#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

PEL_OpenDXwriter const* PEL_OpenDXwriter::PROTOTYPE = new PEL_OpenDXwriter() ;

//----------------------------------------------------------------------
PEL_OpenDXwriter:: PEL_OpenDXwriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "PEL_OpenDXwriter" )
   , file()
   , object_number( 0 )
   , conn_object( 0 )
   , pos_object( 0 )
   , obj_map( 0 )
   , t( 0 )
   , icycle( 0 )
{
}

//----------------------------------------------------------------------
PEL_OpenDXwriter*
PEL_OpenDXwriter:: create_replica( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_OpenDXwriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_OpenDXwriter* result = new PEL_OpenDXwriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PEL_OpenDXwriter:: PEL_OpenDXwriter( PEL_Object* a_owner,
			 	     PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , object_number(0)
   , t(0)
   , icycle(0)
   , binary( false )
{
   filename = exp->string_data( "files_basename" ) ;
   filename += ".dx" ;
   file_pos = 0 ;

   std::string const& format = exp->string_data( "writing_mode" ) ;
   if( format=="binary" )
   {
      binary = true ;
      bin_filename = filename+".bin" ;
      binfile_pos = 0 ;
   }
      
   obj_map = PEL_Map::create( this ) ;

   PEL_CHECK_INV( invariant() ) ;
}


//----------------------------------------------------------------------
PEL_OpenDXwriter:: ~PEL_OpenDXwriter( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_OpenDXwriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_OpenDXwriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;

   icycle = exp->int_data( "cycle_number" ) ;

   if( t.size()<icycle+1 ) t.resize( icycle+1 ) ;   
   t( icycle ) = (double)icycle ;

   // default mode = create in first cycle
   std::ios_base::openmode default_mode = std::ios::out ; 
   if ( file_pos != (std::ofstream::pos_type) 0 )
      default_mode |= std::ios::in ; // except for next cycles

   file.open( filename.c_str(), default_mode ) ;
   if( !file )
   {
      PEL_Error::object()->raise_plain(
         "unable to create the OpenDX output file : " + filename ) ;
   }
   file.seekp( file_pos );
   file << "########## CYCLE " << icycle << " at " << file_pos << std::endl ;

   if( binary )
   {
      binfile.open( bin_filename.c_str(), default_mode | std::ios::binary ) ;
      if( !binfile )
      {
         PEL_Error::object()->raise_plain(
            "unable to create the OpenDX binary output file : " + filename ) ;
      }
      binfile.seekp( binfile_pos );
   }


   if( exp->has_module( "meshing" ) )
   {
      PEL_ModuleExplorer const* se = exp->create_subexplorer( 0, "meshing" ) ;
      write_grid( se ) ;
      se->destroy() ;
   }

   if( exp->has_module( "fields" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "fields" ) ;
      se->start_module_iterator() ;
      for( ; se->is_valid_module() ; se->go_next_module() )
      {
	 PEL_ModuleExplorer* sse = se->create_subexplorer( 0 ) ;
	 write_field( sse ) ;
	 sse->destroy() ;
      }
      se->destroy() ;
   }

   if( exp->has_module( "variables" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "variables" ) ;
      se->start_entry_iterator() ;
      for( ; se->is_valid_entry() ; se->go_next_entry() )
      {
         write_one_variable( se->keyword(), se->data( se ) ) ;
      }
      se->destroy() ;
   }

   if( exp->has_module( "integration_domain" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "integration_domain" ) ;
      se->start_entry_iterator() ;
      for( ; se->is_valid_entry() ; se->go_next_entry() )
      {
         write_one_variable( se->keyword(), se->data( se ) ) ;
      }
      se->destroy() ;
   }

   file_pos = file.flush().tellp() ;

   finish() ;

   file.close() ;

   if( binary )
   {
      binfile_pos  = binfile.tellp() ;
      binfile.close() ;
   }
}

//----------------------------------------------------------------------
void
PEL_OpenDXwriter:: write_grid( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL("PEL_OpenDXwriter::write_grid");

   doubleArray2D const& vertices = exp->doubleArray2D_data( "vertices" ) ;
   pos_object = save_array_with_vectors( vertices ) ;   
   size_t nb_dims = vertices.index_bound( 0 ) ;

   intArray2D const& connectivity = exp->intArray2D_data( "cell2vertex" ) ;
   size_t nvert = connectivity.index_bound(0) ;
   size_t_vector permut( nvert ) ;
   for( size_t i=0 ; i<nvert ; i++ ) permut(i)=i ;
   std::string elm_type ;
   if( nb_dims==3 && nvert==8 )
   {
      elm_type = "cubes" ;
      permut(1)=3 ;
      permut(2)=1 ;
      permut(3)=2 ;
      permut(5)=7 ;
      permut(6)=5 ;
      permut(7)=6 ;
   } 
   else if( nb_dims==3 && nvert==4 )
   {
      elm_type = "tetrahedra" ;
   } 
   else if( nb_dims==1 && nvert==2 )
   {
      elm_type = "lines" ;
   } 
   else if( nb_dims==2 && nvert==4 )
   {
      elm_type = "quads" ;
      permut(1)=3 ;
      permut(2)=1 ;
      permut(3)=2 ;
   } 
   else if( nb_dims==2 && nvert==3 )
   {
      elm_type = "triangles" ;
   } 
   else 
   {
      PEL_Error::object()->raise_plain( "Invalid polyhedron for openDx output " ) ;
   }
   intArray2D I( connectivity ) ;
   for( size_t m = 0 ; m<connectivity.index_bound(1) ; m++ )
   {
      for( size_t l=0 ; l<nvert ; ++l )
      {
         I( l, m ) = connectivity( permut(l), m ) ;
      }
   }

   conn_object = save_array_with_vectors( I ) ;
   file << "attribute \"element type\" string \""
        << elm_type << "\"" << std::endl ;
   file << "attribute \"ref\" string \"positions\"" << std::endl ;
   
}

//----------------------------------------------------------------------
void
PEL_OpenDXwriter:: write_one_variable( std::string const& name,
                                       PEL_Data const* val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_OpenDXwriter:: write_one_variable" ) ;
   
//??????????????? tres tres douteux ?????????????????

   if( name == "TIME" ) 
   {
      t( icycle ) = val->to_double() ;
   }
   else
   {
      size_t obj = 0 ;
      if( val->data_type()==PEL_Data::Double )
      {
         obj = save_scalar( val->to_double() ) ;
      }
      else if( val->data_type()==PEL_Data::DoubleVector )
      {
         obj = save_vector( val->to_double_vector() ) ;
      }
      else if( val->data_type()==PEL_Data::DoubleArray2D )
      {
         obj = save_array( val->to_double_array2D() ) ;
      }
      else if( val->data_type()==PEL_Data::Int )
      {
         obj = save_scalar( val->to_int() ) ;
      }
      else if( val->data_type()==PEL_Data::IntVector )
      {
         obj = save_vector( val->to_int_vector() ) ;
      }
      else if( val->data_type()==PEL_Data::IntArray2D )
      {
         obj = save_array( val->to_int_array2D() ) ;
      }
      else
      {
         PEL_Error::object()->raise_plain(
            "Unsupported data type "+PEL_Data::type_name( val->data_type() )
            +" in PEL_OpenDXwriter:: write_one_variable method" ) ;
      }      
      size_t object = new_object() ;
      file << "object " << object << " class field " << std::endl ;
      file << "component \"data\" value " << obj << std::endl ;
      file << "attribute \"name\" string \"" << name << "\"" << std::endl ;
      file << "attribute \"series position\" number " << icycle << std::endl ;
   
      add_serie( name, object, icycle ) ;

   }
   PEL_CHECK_INV( invariant() ) ;   
}


//----------------------------------------------------------------------
void
PEL_OpenDXwriter:: write_field( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_OpenDXwriter:: write_field" ) ;
   
   doubleArray2D const& X = exp->doubleArray2D_data( "value" ) ;
   std::string const& name = exp->string_data( "name" ) ;
   std::string const& location = exp->string_data( "location" ) ;

   std::string dep ;
   
   if( location == "at_cell_centers" ) 
   {
      dep = "connections" ;
   } 
   else if( location == "at_vertices" ) 
   {
      dep = "positions" ;
   } 
   else
   { 
      raise_field_location_error(
         name, location, "at_vertices,at_cell_centers" ) ;
   }
   size_t obj = save_array_with_vectors( X ) ;   
   file << "attribute \"dep\" string \"" << dep << "\"" << std::endl ;
   size_t field_object = new_object() ;
   file << "object " << field_object << " class field " << std::endl ;
   file << "component \"data\" value " << obj << std::endl ;
   file << "component \"positions\" value " << pos_object << std::endl ;
   file << "component \"connections\" value " << conn_object << std::endl ;
   file << "attribute \"name\" string \"" << name << "\"" << std::endl ;
   file << "attribute \"series position\" number " << icycle << std::endl ;
   
   add_serie( name, field_object, icycle ) ;

}

//-----------------------------------------------------------------------
size_t
PEL_OpenDXwriter:: new_object( void ) 
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_OpenDXwriter:: new_object" ) ;
   
   file << "########## OBJECT " << object_number + 1 << std::endl ;
   return ++object_number ;
}



//-----------------------------------------------------------------------
void
PEL_OpenDXwriter:: write_data_header( void ) 
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_OpenDXwriter:: write_data_header" ) ;
   
   if( binary )
   {
      if( PEL_System::big_endian_encoding() )
      {
         file << " msb " ;
      }
      else
      {
         file << " lsb " ;
      }
      file << " binary data file " << PEL_System::basename( bin_filename )
           << "," << binfile.tellp() << std::endl ;
   }
   else
   {
      file << " data follows " << std::endl ;
   }
}


//-----------------------------------------------------------------------
void
PEL_OpenDXwriter:: write_data_endl( void ) 
//-----------------------------------------------------------------------
{
   if( !binary )
   {
      file << std::endl ;
   }
}

//-----------------------------------------------------------------------
void
PEL_OpenDXwriter:: write_data( int val ) 
//-----------------------------------------------------------------------
{
   static const size_t s = sizeof(int) ;
   
   if( !binary )
   {
      file << " " << val ;
   }
   else
   {
      binfile.write((const char*)&val,s) ;
   }   
}


//-----------------------------------------------------------------------
void
PEL_OpenDXwriter:: write_data( double val ) 
//-----------------------------------------------------------------------
{
   static const size_t s = sizeof(float) ;
   float fval = (float)val ;
   
   if( !binary )
   {
      file << " " << val ;
   }
   else
   {
      binfile.write((char*)&fval,s) ;
   }   
}


//-----------------------------------------------------------------------
size_t
PEL_OpenDXwriter:: save_array( intArray2D const& I ) 
//-----------------------------------------------------------------------
{
   size_t obj = new_object() ;
   size_t dim = I.index_bound(0) ;
   size_t n = I.index_bound(1) ;
   
   file << "object " << obj << " array type int " 
	<< "rank 2 shape " << n << " " << dim 
	<< " items 1 " ;
   
   write_data_header() ;
   
   for( size_t i= 0 ; i<n ; i++ ) {
      for( size_t j= 0 ; j<dim ; j++ ) {
         write_data( I(j, i) ) ;
      }
      write_data_endl() ;
   }
   return obj ;
}



//-----------------------------------------------------------------------
size_t
PEL_OpenDXwriter:: save_vector( intVector const& I ) 
//-----------------------------------------------------------------------
{
   size_t obj = new_object() ;
   size_t n = I.size() ;
   
   file << "object " << obj << " array type int " 
        << "rank 1 shape " << n 
	<< " items 1 " ;
   write_data_header() ;
   for( size_t i= 0 ; i<n ; i++ ) {
      write_data( I(i) ) ;
   }
   write_data_endl() ;
   return obj ;
}



//-----------------------------------------------------------------------
size_t
PEL_OpenDXwriter:: save_scalar( int val ) 
//-----------------------------------------------------------------------
{
   size_t obj = new_object() ;
   
   file << "object " << obj << " array type int " 
        << "rank 0 items 1 data follows " << val << std::endl ;
   return obj ;
}



//-----------------------------------------------------------------------
size_t
PEL_OpenDXwriter:: save_scalar( double val ) 
//-----------------------------------------------------------------------
{
   size_t obj = new_object() ;
   
   file << "object " << obj << " array type float " 
        << "rank 0 items 1 data follows " << val << std::endl ;
   return obj ;
}



//-----------------------------------------------------------------------
size_t
PEL_OpenDXwriter:: save_vector( doubleVector const& D ) 
//-----------------------------------------------------------------------
{
   size_t obj = new_object() ;
   size_t n = D.size();

   file << "object " << obj << " array type float " 
        << "rank 1 shape " << n 
	<< " items 1 " ;
   write_data_header() ;

   for( size_t i=0 ; i<n ; i++ ) {
	 write_data( D(i) ) ;
      }
   write_data_endl() ;

   return obj ;
}

//-----------------------------------------------------------------------
size_t
PEL_OpenDXwriter:: save_array( doubleArray2D const& D ) 
//-----------------------------------------------------------------------
{
   size_t obj = new_object() ;
   size_t dim = D.index_bound(0);
   size_t n = D.index_bound(1);

   file << "object " << obj << " array type float " 
	<< "rank 2 shape " << n << " " << dim 
	<< " items 1 " ;
   write_data_header() ;
   
   for( size_t i=0 ; i<n ; i++ ) {
      for( size_t j=0 ; j<dim ; j++ ) {
	 write_data( D(j,i) ) ;
      }
      write_data_endl() ;
   }

   return obj ;
}


//-----------------------------------------------------------------------
size_t
PEL_OpenDXwriter:: save_array_with_vectors( intArray2D const& I ) 
//-----------------------------------------------------------------------
{
   size_t obj = new_object() ;
   size_t dim = I.index_bound(0);
   size_t n = I.index_bound(1);

   file << "object " << obj << " array type int "
	<< "rank 1 shape " << dim 
	<< " items " << n ;
   write_data_header() ;
   
   for( size_t i=0 ; i<n ; i++ ) {
      for( size_t j=0 ; j<dim ; j++ ) {
	 write_data( I(j,i) ) ;
      }
      write_data_endl() ;
   }

   return obj ;
}


//-----------------------------------------------------------------------
size_t
PEL_OpenDXwriter:: save_array_with_vectors( doubleArray2D const& D ) 
//-----------------------------------------------------------------------
{
   size_t obj = new_object() ;
   size_t dim = D.index_bound(0);
   size_t n = D.index_bound(1);

   file << "object " << obj << " array type float " ;
   if( dim==1 ) {
      file << "rank 0" ;
   } else {
      file << "rank 1 shape " << dim ;
   }
   file << " items " << n ;
   write_data_header() ;

   for( size_t i=0 ; i<n ; i++ ) {
      for( size_t j=0 ; j<dim ; j++ ) {
	 write_data( D(j,i) ) ;
      }
      write_data_endl() ;
   }

   return obj ;
}



//----------------------------------------------------------------------------
void
PEL_OpenDXwriter:: add_serie( std::string const& name,
                             size_t obj,
                             int serie )
//----------------------------------------------------------------------------
{
   PEL_String * str = PEL_String::create( 0, name ) ;
   if( !obj_map->has_key( str ) )
   {
      obj_map->set_item_at( str->create_clone( obj_map ),
                            PEL_List::create( obj_map ) ) ;
   }
   PEL_List*lst = static_cast<PEL_List*>( obj_map->item_at( str ) ) ;
   lst->append( PEL_Int::create( lst, obj ) ) ;
   lst->append( PEL_Int::create( lst, serie ) ) ;
   str->destroy() ;
}



//----------------------------------------------------------------------------
void
PEL_OpenDXwriter:: finish( void )
//----------------------------------------------------------------------------
{
   PEL_Map*serie_map = PEL_Map::create( this ) ;
   
   PEL_MapIterator* it = obj_map->create_iterator( this ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PEL_String const* str = dynamic_cast<PEL_String const*>( it->key() ) ;
      PEL_ASSERT( str!=0 ) ;
      PEL_List const* lst = dynamic_cast<PEL_List const*>( it->item() ) ;
      PEL_ASSERT( lst!=0 ) ;
      size_t serie_object = new_object() ;
      file << "object " << serie_object << " class series" << std::endl ;
      PEL_Iterator * lit = lst->create_iterator( this ) ;
      size_t idx=0 ;
      for( lit->start() ; lit->is_valid() ; lit->go_next() )
      {
         PEL_Int const* obj = dynamic_cast<PEL_Int const* >( lit->item() ) ;
         PEL_ASSERT( obj!=0 ) ;
         lit->go_next() ;
         PEL_Int const* pos =  dynamic_cast<PEL_Int const* >( lit->item() ) ;
         PEL_ASSERT( pos!=0 ) ;
         file << "member " << idx
              << " position " << t( pos->to_int() )
              << " value " << obj->to_int() << std::endl ;
         idx++ ;
      }
      serie_map->set_item_at( str->create_clone( this ),
                              PEL_Int::create( this, serie_object ) ) ;
   }
   
   file << "object \"default\" class group" << std::endl ;
   it = serie_map->create_iterator( this ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PEL_String const* str = dynamic_cast<PEL_String const*>( it->key() ) ;
      PEL_ASSERT( str!=0 ) ;
      PEL_Int const* n = dynamic_cast<PEL_Int const*>( it->item() ) ;
      PEL_ASSERT( n!=0 ) ;
      file << "member \"" << str->to_string() << "\" value "
           << n->to_int() << std::endl ;
   }
   file << "end" << std::endl ;  
}

