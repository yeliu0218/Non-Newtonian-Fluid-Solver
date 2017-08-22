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

#include <PEL_TICwriter.hh>

#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <intArray2D.hh>
#include <intVector.hh>
#include <stringVector.hh>

PEL_TICwriter const* PEL_TICwriter::PROTOTYPE = new PEL_TICwriter() ;

//----------------------------------------------------------------------
PEL_TICwriter:: PEL_TICwriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "PEL_TICwriter" )
{
}

//----------------------------------------------------------------------
PEL_TICwriter*
PEL_TICwriter:: create_replica( PEL_Object* a_owner,
   			        PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICwriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_TICwriter* result = new PEL_TICwriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_TICwriter:: PEL_TICwriter( PEL_Object* a_owner,
			       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , FILE_NAME( exp->string_data( "files_basename" )+".gene" )
   , FORMAT( PEL_TICio::Unspecified )
   , APPEND_MODE( false )
{
   // Binary save :
   std::string const& format = exp->string_data( "writing_mode" ) ;
   if( format == "text" )
   {
      FORMAT = PEL_TICio::Text ;
   }
   else if( format == "binary" )
   {
      FORMAT = PEL_TICio::Binary ;
   }
   else if ( format == "binary_no_local" )
   {
      FORMAT = PEL_TICio::Binary_no_local ;
   }
   else if( format == "Cbinary" )
   {
      FORMAT = PEL_TICio::CBinary ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value(
         exp,
         "writing_mode",
         "   - \"text\"\n"
         "   - \"binary\"\n"
         "   - \"binary_no_local\"\n"
         "   - \"Cbinary\"" ) ;
   }
  
   if( exp->has_entry( "append_mode" ) )
   {
      APPEND_MODE = exp->bool_data( "append_mode" ) ;
   }
   if( !APPEND_MODE )
   {
      PEL_TICio::create_file( FILE_NAME, FORMAT ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_TICwriter:: ~PEL_TICwriter( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PEL_TICwriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICwriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;

   static int cycle_number = 0 ;

   if( APPEND_MODE )
   {
      cycle_number++ ;
   }
   else
   {
      cycle_number = exp->int_data( "cycle_number" ) ;
   }

   int nbv = nb_variables( exp ) ;
   PEL_TICio::write_new_cycle( FILE_NAME, FORMAT, nbv, cycle_number ) ;

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
PEL_TICwriter:: nb_variables( PEL_ModuleExplorer const* exp ) const
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
PEL_TICwriter:: write_grid( PEL_ModuleExplorer const* exp ) const 
//----------------------------------------------------------------------
{
   doubleArray2D const& vertices = exp->doubleArray2D_data( "vertices" ) ;
   size_t const dimension = vertices.index_bound(0) ;
   size_t const nbvertices = vertices.index_bound(1) ;

   intArray2D const& connectivity = exp->intArray2D_data( "cell2vertex" ) ;
   size_t const nbmeshes = connectivity.index_bound(1) ;

   intVector const& cell_nb_vertices =
                             exp->intVector_data( "cell_nb_vertices" ) ;
   
   intVector linearConnectivity( 0 ) ;
   for( size_t i=0 ; i<nbmeshes ; i++ ) 
   {
      size_t meshsize = cell_nb_vertices(i) ;
      for (size_t j=0 ; j<meshsize ; j++ ) 
      {
	 linearConnectivity.append( connectivity( j, i ) + 1 ) ;
      }
   }

   intVector ngr(0) ;
   int cell_nb_verts = cell_nb_vertices(0) ;
   int n = 1 ;
   for( size_t i=1 ; i<cell_nb_vertices.size() ; ++i )
   {
      if( cell_nb_vertices(i) != cell_nb_verts )
      {
         ngr.append( n ) ;
         ngr.append( cell_nb_verts ) ;
         n = 1 ;
         cell_nb_verts = cell_nb_vertices(i) ;
      }
      else
      {
         ++n ;
      }
   }
   ngr.append( n ) ;
   ngr.append( cell_nb_verts ) ;

   if( n == (int) nbmeshes )
   {
      PEL_TICio::write_int_variable( FILE_NAME, FORMAT,
                                     "NVER", cell_nb_verts ) ;
   }
   PEL_TICio::write_intVector_variable( FILE_NAME, FORMAT,
                                        "NGR", ngr ) ;
   PEL_TICio::write_intVector_variable( FILE_NAME, FORMAT,
                                        "NODE", linearConnectivity ) ;

   doubleVector coords( nbvertices ) ;

   for( size_t i=0 ; i<nbvertices ; i++ ) coords( i ) = vertices( 0, i ) ;
   PEL_TICio::write_doubleVector_variable( FILE_NAME, FORMAT,
                                           "XNOD", coords ) ;

   if( dimension >= 2 )
   {
      for( size_t i=0 ; i<nbvertices ; i++ ) coords( i ) = vertices( 1, i ) ;
      PEL_TICio::write_doubleVector_variable( FILE_NAME, FORMAT,
                                              "ZNOD", coords ) ;
   }

   if( dimension >= 3 )
   {
      for( size_t i=0 ; i<nbvertices ; i++ ) coords( i ) = vertices( 2, i ) ;
      PEL_TICio::write_doubleVector_variable( FILE_NAME, FORMAT,
                                              "YNOD", coords ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_TICwriter:: write_field( PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   doubleArray2D const& X = exp->doubleArray2D_data( "value" ) ;
   std::string const& name = exp->string_data( "name" ) ;

   std::string const& location = exp->string_data( "location" ) ;
   if( location == "at_vertices" || location == "at_cell_centers" )
   {
      doubleVector X_linear( X.index_bound(0)*X.index_bound(1) ) ;
      for( size_t i=0 ; i<X.index_bound(0) ; i++ )
      {
         for( size_t j=0 ; j<X.index_bound(1) ; j++ )
         {
            X_linear( i + j*X.index_bound(0) ) = X( i, j ) ;
         }
      }
      PEL_TICio::write_doubleVector_variable( FILE_NAME, FORMAT,
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
PEL_TICwriter:: write_integration_domain(
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
   
   PEL_TICio::write_doubleVector_variable( FILE_NAME, FORMAT, "XFS", xFS ) ;
   PEL_TICio::write_doubleVector_variable( FILE_NAME, FORMAT, "ZFS", zFS ) ;
   PEL_TICio::write_doubleVector_variable( FILE_NAME, FORMAT, "XDOM", xDOM ) ;
   PEL_TICio::write_doubleVector_variable( FILE_NAME, FORMAT, "ZDOM", zDOM ) ;
}

//----------------------------------------------------------------------
void
PEL_TICwriter:: write_one_variable( std::string const& name,
                                    PEL_Data const* val ) const
//----------------------------------------------------------------------
{
   if( val->data_type()==PEL_Data::Int )
   {
      PEL_TICio::write_int_variable( FILE_NAME, FORMAT,
                                     name, val->to_int() ) ;
   }
   else if( val->data_type()==PEL_Data::Double )
   {
      PEL_TICio::write_double_variable( FILE_NAME, FORMAT,
                                        name, val->to_double() ) ;
   }
   else if( val->data_type()==PEL_Data::DoubleVector )
   {
      PEL_TICio::write_doubleVector_variable( FILE_NAME, FORMAT,
                                              name, val->to_double_vector() ) ;
   }
   else if( val->data_type()==PEL_Data::IntVector )
   {
      PEL_TICio::write_intVector_variable( FILE_NAME, FORMAT,
                                           name, val->to_int_vector() ) ;
   }
   else if( val->data_type()==PEL_Data::DoubleArray2D )
   {
      PEL_TICio::write_doubleArray2D_variable( FILE_NAME, FORMAT,
                                               name, val->to_double_array2D() ) ;
   }
   else if( val->data_type()==PEL_Data::IntArray2D )
   {
      PEL_TICio::write_intArray2D_variable( FILE_NAME, FORMAT,
                                            name, val->to_int_array2D() ) ;
   }
   else 
   {
      std::string mesg = "*** PEL_TICwriter error :\n" ;
      mesg += "   cannot write variables of type \"" ;
      mesg += PEL_Data::type_name( val->data_type() ) ;
      mesg += "\"" ;
      PEL_Error::object()->raise_plain( mesg ) ;
   }
}
