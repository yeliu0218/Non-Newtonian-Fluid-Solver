/*
 *  Copyright :
 *    "Institut de Radioprotection et de Sret�Nucl�ire - IRSN" (1995-2008)
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

#include <PEL_MATLABwriter.hh>

#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <PEL_DoubleVector.hh>

#include <doubleArray2D.hh>
//#include <doubleVector.hh>
#include <intArray2D.hh>
#include <intVector.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>

PEL_MATLABwriter const* PEL_MATLABwriter::PROTOTYPE = new PEL_MATLABwriter() ;

//----------------------------------------------------------------------
PEL_MATLABwriter:: PEL_MATLABwriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "PEL_MATLABwriter" )
   , ICYCLE( 0 )
   , COORDS(0)
{
}

//----------------------------------------------------------------------
PEL_MATLABwriter*
PEL_MATLABwriter:: create_replica( PEL_Object* a_owner,
   			        PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABwriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_MATLABwriter* result = new PEL_MATLABwriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_MATLABwriter:: PEL_MATLABwriter( PEL_Object* a_owner,
			       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , FILEBASENAME( exp->string_data( "files_basename" ) )
   , FILEXTENSION(".dat")
   , FORMAT( PEL_MATLABio::Unspecified )
   , CONTEXT( PEL_ContextSimple::create( this ) )
   , COORDS( PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) )
   {

   std::string const& format = exp->string_data( "writing_mode" ) ;
   if( format == "text" )
   {
      FORMAT = PEL_MATLABio::Text ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value(
         exp,
         "writing_mode",
         "   - \"text\"\n" ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_MATLABwriter:: ~PEL_MATLABwriter( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_MATLABwriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABwriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;

   // Current time
   ICYCLE = exp->int_data( "cycle_number" ) ;

   // Open files and test
   std::string fname = output_file_name( ICYCLE, "T" ) ;
   std::string fnameU = output_file_name( ICYCLE, "U" ) ;
   std::string fnameSTRESS = output_file_name( ICYCLE, "S" ) ;
   FILE_NAME=fname;
   FILE_NAME_U = fnameU; 
   FILE_NAME_STRESS = fnameSTRESS;
   PEL_MATLABio::create_file( fname, FORMAT ) ;
   PEL_MATLABio::create_file( fnameU, FORMAT ) ;
   PEL_MATLABio::create_file( fnameSTRESS, FORMAT ) ;

   //int nbv = nb_variables( exp ) ;
   //PEL_MATLABio::write_new_cycle( fname, FORMAT, nbv, ICYCLE ) ;

   if( exp->has_module( "meshing" ) )
   {
	  PEL_ModuleExplorer const* se = exp->create_subexplorer( 0, "meshing" ) ;
      write_grid( se ) ;
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
   
   if( exp->has_module( "integration_domain" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "integration_domain" ) ;
      write_integration_domain( se ) ;
      se->destroy() ;
   }
   
}

//----------------------------------------------------------------------
int
PEL_MATLABwriter:: nb_variables( PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   int result = 0 ;

   result += exp->int_data( "nb_variables" ) ;
   result += exp->int_data( "nb_fields" ) ;
   if( exp->has_module( "meshing" ) )
   {
      result += 2 + exp->int_data( "meshing/nb_sp_dims" ) ;
      intVector const& cell_nb_vertices =
         exp->intVector_data( "meshing/cell_nb_vertices" ) ;
      bool ok = true ;
      int n = cell_nb_vertices(0) ;
      for( size_t i=1 ; ok && i<cell_nb_vertices.size() ; ++i )
      {
         ok = ( n == cell_nb_vertices(i) ) ;
      }
      if( ok ) ++result ;
   }
   if( exp->has_module( "integration_domain" ) )
   {
      result += 4 ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_MATLABwriter:: write_grid( PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   doubleArray2D const& vertices = exp->doubleArray2D_data( "vertices" ) ;
   size_t const dimension = vertices.index_bound(0) ;
   size_t const nbvertices = vertices.index_bound(1) ;

   intArray2D const& connectivity = exp->intArray2D_data( "cell2vertex" ) ;
   size_t const nb_verts_per_cell = connectivity.index_bound(0) ;
   size_t const my_nb_cells = connectivity.index_bound(1) ;


	   doubleArray2D  cell_centers(dimension,my_nb_cells);

	   for( size_t f=0 ; f<my_nb_cells ; f++ ){
		   cell_centers( 0, f )=0;
		   cell_centers( 1, f )=0;
	         for( size_t i=0 ; i<nb_verts_per_cell ; i++ ){
	        	 cell_centers( 0, f ) += vertices(0,connectivity( i, f )) ;
	        	 cell_centers( 1, f ) += vertices(1,connectivity( i, f )) ;
//	        	 PEL::out() << i << ","<< f <<"  "<< vertices(0,connectivity( i, f ))<< std::endl;
	        	 }
//	         PEL::out() << "------------"<< std::endl;
	         cell_centers( 0, f )=cell_centers( 0, f )/nb_verts_per_cell;
	         cell_centers( 1, f )=cell_centers( 1, f )/nb_verts_per_cell;
	   }

       doubleVector coordsCellX( my_nb_cells ) ;
	   doubleVector coordsCellY( my_nb_cells ) ;
	   doubleVector coordsCellZ( my_nb_cells ) ;
	   for( size_t i=0 ; i<my_nb_cells ; i++ ) {
			   coordsCellX( i ) = cell_centers( 0, i );
	   }

	   if( dimension >= 2 )
	   {
		   for( size_t j=0 ; j<my_nb_cells ; j++ ) coordsCellY( j ) = cell_centers( 1, j ) ;
	   }
	   if( dimension >= 3 )
	   {
		   for( size_t k=0 ; k<my_nb_cells ; k++ ) coordsCellZ( k ) = cell_centers( 2, k );
	   }
	   PEL_MATLABio::write_gridXY_variable( FILEBASENAME, FORMAT, dimension,
	                                          "_cell" ,coordsCellX, coordsCellY, coordsCellZ ) ;

	   doubleVector coordsX( nbvertices ) ;
	   doubleVector coordsY( nbvertices ) ;
	   doubleVector coordsZ( nbvertices ) ;
	   for( size_t i=0 ; i<nbvertices ; i++ ) coordsX( i ) = vertices( 0, i ) ;

	   if( dimension >= 2 )
	   {
		   for( size_t j=0 ; j<nbvertices ; j++ ) coordsY( j ) = vertices( 1, j ) ;
	   }
	   if( dimension >= 3 )
	   {
		   for( size_t k=0 ; k<nbvertices ; k++ ) coordsZ( k ) = vertices( 2, k );
	   }
	   PEL_MATLABio::write_gridXY_variable( FILEBASENAME, FORMAT, dimension,
	                                          "_vert" ,coordsX, coordsY, coordsZ ) ;


}


//----------------------------------------------------------------------
void
PEL_MATLABwriter:: write_field( PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   doubleArray2D const& X = exp->doubleArray2D_data( "value" ) ;
   std::string const& name = exp->string_data( "name" ) ;

   std::string const& location = exp->string_data( "location" ) ;

   if( location == "at_vertices" )
   {
      doubleVector X_linear( X.index_bound(0)*X.index_bound(1) ) ;
      for( size_t i=0 ; i<X.index_bound(0) ; i++ )
      {
	  for( size_t j=0 ; j<X.index_bound(1) ; j++ )
	  {
	    X_linear( i + j*X.index_bound(0) ) = X( i, j ) ;
	  }
      }
      if( name == "STRESS" || name == "GAMMA" || name == "GAMMADOT")
	PEL_MATLABio::write_doubleVector_variableSTRESS(FILE_NAME_STRESS , FORMAT,
							name, X_linear ) ;
      else
	PEL_MATLABio::write_doubleVector_variableUU(FILE_NAME_U , FORMAT,
						    name, X_linear ) ;
   }	
   else if( location == "at_cell_centers" )
   {
      doubleVector X_linear( X.index_bound(0)*X.index_bound(1) ) ;
      for( size_t i=0 ; i<X.index_bound(0) ; i++ )
      {
         for( size_t j=0 ; j<X.index_bound(1) ; j++ )
         {
            X_linear( i + j*X.index_bound(0) ) = X( i, j ) ;
         }
      }
	if( name == "ExtForce")
	  PEL_MATLABio::write_doubleVector_variableUU( FILE_NAME_U, FORMAT,
		                               name, X_linear ) ;
	else
      PEL_MATLABio::write_doubleVector_variable( FILE_NAME, FORMAT,
                                           name, X_linear ) ;

   }
   else
   {
      raise_field_location_error(
         name, location, "at_vertices,at_cell_centers" ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_MATLABwriter:: write_integration_domain(
                                   PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   doubleArray2D const& inner_boundary =
                            exp->doubleArray2D_data( "inner_boundary" ) ;
   doubleVector xFS( inner_boundary.index_bound(1) ) ;
   doubleVector zFS( inner_boundary.index_bound(1) ) ;
   inner_boundary.extract_section( 0, 0, xFS ) ;
   inner_boundary.extract_section( 0, 1, zFS ) ;

   doubleArray2D const& polygon = exp->doubleArray2D_data( "polygon" ) ;
   doubleVector xDOM( polygon.index_bound(1) ) ;
   doubleVector zDOM( polygon.index_bound(1) ) ;
   polygon.extract_section( 0, 0, xDOM ) ;
   polygon.extract_section( 0, 1, zDOM ) ;

   PEL_MATLABio::write_doubleVector_variable( FILE_NAME, FORMAT, "XFS", xFS ) ;
   PEL_MATLABio::write_doubleVector_variable( FILE_NAME, FORMAT, "ZFS", zFS ) ;
   PEL_MATLABio::write_doubleVector_variable( FILE_NAME, FORMAT, "XDOM", xDOM ) ;
   PEL_MATLABio::write_doubleVector_variable( FILE_NAME, FORMAT, "ZDOM", zDOM ) ;
}

//----------------------------------------------------------------------
void
PEL_MATLABwriter:: write_one_variable( std::string const& name,
                                    PEL_Data const* val ) const
//----------------------------------------------------------------------
{
   if( val->data_type()==PEL_Data::Int )
   {
      PEL_MATLABio::write_int_variable( FILE_NAME, FORMAT,
                                     name, val->to_int() ) ;
   }
   else if( val->data_type()==PEL_Data::Double )
   {
      PEL_MATLABio::write_double_variable( FILE_NAME, FORMAT,
                                        name, val->to_double() ) ;
   }
   else if( val->data_type()==PEL_Data::DoubleVector )
   {
      PEL_MATLABio::write_doubleVector_variable( FILE_NAME, FORMAT,
                                              name, val->to_double_vector() ) ;
   }
   else if( val->data_type()==PEL_Data::IntVector )
   {
      PEL_MATLABio::write_intVector_variable( FILE_NAME, FORMAT,
                                           name, val->to_int_vector() ) ;
   }
   else if( val->data_type()==PEL_Data::DoubleArray2D )
   {
      PEL_MATLABio::write_doubleArray2D_variable( FILE_NAME, FORMAT,
                                               name, val->to_double_array2D() ) ;
   }
   else if( val->data_type()==PEL_Data::IntArray2D )
   {
      PEL_MATLABio::write_intArray2D_variable( FILE_NAME, FORMAT,
                                            name, val->to_int_array2D() ) ;
   }
   else
   {
      std::string mesg = "*** PEL_MATLABwriter error :\n" ;
      mesg += "   cannot write variables of type \"" ;
      mesg += PEL_Data::type_name( val->data_type() ) ;
      mesg += "\"" ;
      PEL_Error::object()->raise_plain( mesg ) ;
   }
}

//-----------------------------------------------------------------------
std::string
PEL_MATLABwriter:: output_file_name( size_t nb, std::string add_string )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GMVwriter:: output_file_name" ) ;
   PEL_CHECK( nb<9999 ) ;

   std::ostringstream tmp ;
   tmp << nb ;
   std::string nb_string = tmp.str() ;
   std::string result = FILEBASENAME+add_string+"0000"+FILEXTENSION;
   result.replace( result.length()-FILEXTENSION.length()-nb_string.length(),
                   nb_string.length(),
                   nb_string ) ;

   PEL_CHECK( result.length()==( FILEBASENAME.length()+add_string.length()
                                +FILEXTENSION.length()+6) ) ;

   return( result ) ;
}
