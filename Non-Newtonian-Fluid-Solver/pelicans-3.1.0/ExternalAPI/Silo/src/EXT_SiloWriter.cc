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

#include <EXT_SiloWriter.hh>

#include <silo.h>

#include <PEL_Bool.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_assertions.hh>
#include <doubleArray2D.hh>
#include <intArray2D.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>
#include <iomanip>

EXT_SiloWriter const* EXT_SiloWriter::PROTOTYPE = new EXT_SiloWriter() ;

//----------------------------------------------------------------------
EXT_SiloWriter:: EXT_SiloWriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "EXT_SiloWriter" )
   , dbfile( 0 )
{
   PEL_Bool* val = PEL_Bool::create( 0, true ) ;
   PEL_Exec::add_variable_to_execution_context(
                         PEL_Variable::object( "BS_with_SILO" ), val ) ;
}

//----------------------------------------------------------------------
EXT_SiloWriter*
EXT_SiloWriter:: create_replica( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SiloWriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   EXT_SiloWriter* result = new EXT_SiloWriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
EXT_SiloWriter:: EXT_SiloWriter( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , dbfile( 0 )
{

   std::string filename = exp->string_data( "files_basename" ) ;
   filename += ".silo" ;
   char* pathname = const_cast<char*>( filename.c_str() ) ;

   dbfile = DBCreate( pathname, DB_CLOBBER, DB_LOCAL, NULL, DB_PDB ) ;
   if( dbfile == 0 ) {
      PEL_Error::object()->raise_plain( "unable to create the Silo database" ) ;
   }
   dir_name = "" ;
   meshing_name = "" ;
   
   PEL_CHECK_INV( invariant() ) ;
}


//----------------------------------------------------------------------
EXT_SiloWriter:: ~EXT_SiloWriter( void )
//----------------------------------------------------------------------
{
   if( !is_a_prototype() )
   {
      int error = DBClose( dbfile ) ;
      if( error != 0 )
         PEL_Error::object()->raise_plain( "unable to close the Silo database" ) ;
   }
   else
   {
      PROTOTYPE = 0  ;
   }
}

//----------------------------------------------------------------------
void
EXT_SiloWriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SiloWriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;

   int cycle_number = exp->int_data( "cycle_number" ) ;
   std::ostringstream os ;
   os <<  "/cycle" << std::setw(4) << std::setfill('0') << cycle_number ;
   dir_name = os.str() ;
   DBMkDir( dbfile, const_cast<char*>( dir_name.c_str() ) ) ;
   DBSetDir( dbfile, const_cast<char*>( dir_name.c_str() ) ) ;

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
         write_field( cycle_number, sse ) ;
	 sse->destroy() ;
      }
      se->destroy() ;
   }
}

//----------------------------------------------------------------------
void
EXT_SiloWriter:: write_grid( PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SiloWriter:: write_grid" ) ;
   
   doubleArray2D const& vertices = exp->doubleArray2D_data( "vertices" ) ;
   size_t const ndims = vertices.index_bound(0) ;
   size_t const nnodes = vertices.index_bound(1) ;

   intArray2D const& connectivity = exp->intArray2D_data( "cell2vertex" ) ;
   size_t const nzones = connectivity.index_bound(1) ;
   size_t const nvert = connectivity.index_bound(0) ;

   int nzshapes = 1 ;
   int zshapecnt  [] = { nzones } ; 
   int zshapesize [] = { nvert } ; 
      
   int lznodelist = nzones*zshapesize[0] ;
   int* znodelist = new int[ lznodelist ] ;

   for (size_t m=0; m<nzones; m++) 
   {
      for (size_t n=0; n<nvert; n++) 
      {
	 znodelist[m*nvert+n] = connectivity(n, m);
      }
   }
   std::ostringstream os ;
   os <<  dir_name << "//mesh" ;
   
   meshing_name = os.str() ;

   char zonelname [] = "zonelist" ;

   int error ;
   error = DBPutZonelist( dbfile,
                          zonelname,
                          nzones,
                          ndims,
                          znodelist,
                          lznodelist,
                          0,
                          zshapesize,
                          zshapecnt,
                          nzshapes ) ;

   if( error != 0 ) 
   {
      PEL_Error::object()->raise_plain( "DBPutZoneList failure" ) ;
   }

   double** coords = new double* [ndims] ;
   for( size_t i=0 ; i<ndims ; ++i ) 
   {
      coords[i] = new double[ nnodes ] ;
   }

   for( size_t j=0 ; j<nnodes ; ++j ) 
   {
      for( size_t ic=0 ; ic<ndims ; ++ic )
	 coords[ic][j] = vertices(ic, j) ;
   }

   error = DBPutUcdmesh( dbfile, 
                         const_cast<char*>(meshing_name.c_str()),
                         ndims, 
			 NULL,
			 (float**)coords,
                         nnodes,
                         nzones,
                         zonelname,
                         NULL, //facelname, 
                         DB_DOUBLE,
                         NULL ) ;
   if( error != 0 ) 
   {
      PEL_Error::object()->raise_plain( "DBPutUcdmesh failure" ) ;
   }

    delete [] znodelist ;
    for( size_t i=0 ; i<ndims ; ++i ) 
    {
       delete [] coords[i] ;
    }
    delete[] coords ;
}

//----------------------------------------------------------------------
void
EXT_SiloWriter:: write_field( int cycle_number,
                              PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SiloWriter:: write_field" ) ;
   
   doubleArray2D const& X = exp->doubleArray2D_data( "value" ) ;
   std::string const& name = exp->string_data( "name" ) ;
   std::string const& location = exp->string_data( "location" ) ;

   int error = 0 ;

   DBoptlist* optlist = DBMakeOptlist( 10 ) ;
   int opt = cycle_number ;
   DBAddOption( optlist, DBOPT_CYCLE, &opt ) ;

   size_t nels = X.index_bound( 1 ) ;

   int centering = 0 ;
   if( location == "at_cell_centers" ) 
   {
      centering = DB_ZONECENT ;
   } 
   else if( location == "at_vertices" ) 
   {
      centering = DB_NODECENT ;
   } 
   else 
   {
      raise_field_location_error(
         name, location, "at_vertices,at_cell_centers" ) ;
   }

   if( X.index_bound(0) == 1 ) 
   {
      error = -1 ;
      double* var = new double[ nels ] ;
      for( size_t i=0 ; i<nels ; ++i ) var[i] = X(0,i) ;

      error = DBPutUcdvar1( dbfile,
                            const_cast<char*>(name.c_str()),
              	            const_cast<char*>(meshing_name.c_str()),
			    (float*)var,
                            nels,
                            NULL,
                            0,
                            DB_DOUBLE,
                            centering,
                            optlist ) ;

      delete [] var ;
    } 
    else 
    {
       size_t nvars = X.index_bound(0) ;
       char** varnames = new char* [ nvars ] ;
       double** vars = new double* [ nvars ] ;
       for( size_t i=0 ; i<nvars ; i++ ) 
       {
          varnames[i] = new char[ 512 ] ;
          std::ostringstream os ;
          os << name << "_" << std::setw(2) << std::setfill('0') << i ;
          for( size_t c=0 ; c<=os.str().length() ; c++ )
             varnames[i][c] = os.str()[c] ;
          
          vars[i] = new double [ nels ] ;
          for( size_t j=0 ; j<nels ; ++j ) vars[i][j]=X(i,j) ;
       }
       error = DBPutUcdvar( dbfile,
                            const_cast<char*>(name.c_str()),
                            const_cast<char*>(meshing_name.c_str()),
                            nvars,
                            varnames,
                            (float**)vars,
                            nels,
                            NULL,
                            0,
                            DB_DOUBLE,
                            centering,
                            optlist ) ;

        for( size_t i=0 ; i<nvars ; ++i ) 
        {
           delete [] varnames[i] ;
           delete [] vars[i] ;
        } 
        delete [] varnames ;
        delete [] vars ;
     }

     DBFreeOptlist( optlist ) ;

     if( error != 0 ) 
     {
         PEL_Error::object()->raise_plain( " unable to write variables into a Silo file" ) ;
     }

}
