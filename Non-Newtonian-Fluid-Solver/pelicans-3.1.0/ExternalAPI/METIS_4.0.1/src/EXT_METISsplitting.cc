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

#include <EXT_METISsplitting.hh>

#include <PEL_Bool.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Timer.hh>
#include <PEL_Variable.hh>
#include <PEL.hh>

#include <intVector.hh>

#include <GE_Meshing.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_ReferenceTriangle.hh>
#include <GE_ReferenceTetrahedron.hh>
#include <GE_ReferenceCube.hh>
#include <GE_ReferenceSquare.hh>

#include <iostream>
#include <string>
#include <sstream>

extern "C" {
#include <metis.h>
}

EXT_METISsplitting const* 
EXT_METISsplitting::PROTOTYPE = new EXT_METISsplitting() ;

struct EXT_METISsplitting_ERROR
{
   static void n0( void ) ;
   static void n1( std::string const& s1, std::string const& s2 ) ;
   static void n2( int rk, size_t nb_rks ) ;
   static void n3( std::string const& cell_poly_name ) ;
   static void n4( int rk ) ;
} ;

//-------------------------------------------------------------------------
EXT_METISsplitting:: EXT_METISsplitting( void )
//------------------------------------------------------------------------
   : GE_SplittingStrategy( "EXT_METISsplitting" )
   , CELL_RANK( 0 )
   , TIMER( 0 )
   , RANK_CELLS( 0 )
{
   PEL_Bool* val = PEL_Bool::create( 0, true ) ;
   PEL_Exec::add_variable_to_execution_context(
                         PEL_Variable::object( "BS_with_METIS" ), val ) ;
}

//-------------------------------------------------------------------------
EXT_METISsplitting*
EXT_METISsplitting:: create_replica( PEL_Object* a_owner,
				     PEL_ModuleExplorer const* exp,
				     GE_Meshing* meshing,
				     PEL_Communicator const* com ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_METISsplitting:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, meshing, com ) ) ;

   EXT_METISsplitting* result =
      new EXT_METISsplitting( a_owner, exp, meshing, com ) ;

   PEL_CHECK( create_replica_POST( a_owner, exp, meshing, com, result ) ) ;
   return ( result ) ;
}

//-------------------------------------------------------------------------
EXT_METISsplitting:: EXT_METISsplitting( PEL_Object* a_owner,
					 PEL_ModuleExplorer const* exp,
					 GE_Meshing* meshing,
					 PEL_Communicator const* com )
//------------------------------------------------------------------------
   : GE_SplittingStrategy( a_owner, exp, meshing, com )
   , CELL_RANK( meshing->nb_cells() )
   , TIMER( 0 )
   , RANK_CELLS( 0 )
{
   PEL_LABEL( "EXT_METISsplitting:: EXT_METISsplitting" ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const nb_rks = com->nb_ranks() ;
   size_t const rk = com->rank() ;

   CELL_RANK.set( 0 ) ;
   if( nb_rks>1 && rk == 0 )
   {
      TIMER = PEL_Timer::create( this ) ;
      
      // CALL METIS
      METIS_balancing( meshing ) ;
   }

   com->broadcast( CELL_RANK, 0 ) ;
}

//-------------------------------------------------------------------------
EXT_METISsplitting*
EXT_METISsplitting:: create_replica( PEL_Object* a_owner,
				     PEL_ModuleExplorer const* exp,
				     GE_Meshing* meshing,
				     size_t nb_rks, size_t rk ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_METISsplitting:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, meshing, nb_rks, rk ) ) ;

   EXT_METISsplitting* result =
      new EXT_METISsplitting( a_owner, exp, meshing, nb_rks, rk ) ;

   PEL_CHECK( create_replica_POST( a_owner, exp, meshing, nb_rks, rk, 
                                   result ) ) ;
   return ( result ) ;
}

//-------------------------------------------------------------------------
EXT_METISsplitting:: EXT_METISsplitting( PEL_Object* a_owner,
					 PEL_ModuleExplorer const* exp,
					 GE_Meshing* meshing,
					 size_t nb_rks, size_t rk )
//------------------------------------------------------------------------
   : GE_SplittingStrategy( a_owner, exp, meshing, nb_rks, rk )
   , CELL_RANK( meshing->nb_cells() )
   , TIMER( 0 )
   , RANK_CELLS( 0 )
{
   PEL_LABEL( "EXT_METISsplitting:: EXT_METISsplitting" ) ;
   PEL_CHECK_INV( invariant() ) ;

   CELL_RANK.set( 0 ) ;
   if( nb_rks>1 )
   {
      TIMER = PEL_Timer::create( this ) ;
      
      // CALL METIS
      METIS_balancing( meshing ) ;
   }
}

//-------------------------------------------------------------------------
EXT_METISsplitting:: ~EXT_METISsplitting( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
size_t 
EXT_METISsplitting:: cell_rank( size_t mesh_id ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_METISsplitting:: cell_rank" ) ;
   PEL_CHECK_PRE( cell_rank_PRE( mesh_id ) ) ;

   size_t result = (size_t) CELL_RANK( mesh_id ) ;
 
   PEL_CHECK_POST( cell_rank_POST( mesh_id, result ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
EXT_METISsplitting:: print( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_METISsplitting:: print" ) ;
   
   if( TIMER != 0 )
   {
      std::string const s( indent_width, ' ' ) ;
      os << s << "METIS partitioning time: " ;
      PEL_Timer::print_time( TIMER->elapsed_time(), os, 0 ) ;
      os << std::endl ;
      int min_val = (int) nb_cells() ;
      int max_val = 0 ;
      for( size_t i=0 ; i<nb_ranks() ; ++i )
      {
         min_val = PEL::min( min_val, RANK_CELLS(i) ) ;
         max_val = PEL::max( max_val, RANK_CELLS(i) ) ;
      }
      os << s << "   maximal number of cells per processes: "
         << max_val << std::endl ;
      os << s << "   minimal number of cells per processes: "
         << min_val << std::endl ;
   }
}

//------------------------------------------------------------------------
void 
EXT_METISsplitting:: METIS_balancing( GE_Meshing* meshing )
//------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_METISsplitting:: METIS_balancing" ) ;
   PEL_CHECK( meshing != 0 ) ;
   PEL_CHECK( meshing->nb_cells() == nb_cells() ) ;
   PEL_CHECK( nb_ranks()>1 ) ;
   PEL_CHECK( TIMER != 0 ) ;

   size_t const nb_sp_dims = meshing->nb_space_dimensions() ;
   if( nb_sp_dims != 2 && nb_sp_dims != 3 )
   {
      EXT_METISsplitting_ERROR::n0() ;
   }
   
   int ne = meshing->nb_cells() ;
   int nn = meshing->nb_vertices() ;
   int numflag = 0 ;
   int nparts = nb_ranks() ;
   int edgecut = PEL::bad_index() ;
   idxtype* epart = new idxtype[ne] ;
   idxtype* npart = new idxtype[nn] ;

   meshing->start_cell_iterator() ;
   std::string cell_poly = meshing->cell_polyhedron_name() ;
   int etype = METIS_cell_type( cell_poly ) ;
   size_t nb_vertices_per_cell 
                = meshing->cell_reference_polyhedron()->nb_vertices() ;
   idxtype* elmnts = new idxtype[ne*nb_vertices_per_cell] ;

   size_t cell_counter = 0 ;
   for( ; meshing->valid_cell() ; meshing->go_next_cell(), ++cell_counter )
   {
      if( meshing->cell_polyhedron_name() != cell_poly )
      {
         EXT_METISsplitting_ERROR::n1(
                            cell_poly, meshing->cell_polyhedron_name() ) ;
      } 
      size_t_vector const& cv = meshing->cell_vertices() ;
      PEL_CHECK( cv.size()==nb_vertices_per_cell ) ;
      size_t i_shift = cell_counter*nb_vertices_per_cell ;
      if( etype == 4 ) // quadrilateral
      {
         elmnts[i_shift+0]=(idxtype)cv( 0 ) ;
         elmnts[i_shift+1]=(idxtype)cv( 3 ) ;
         elmnts[i_shift+2]=(idxtype)cv( 2 ) ;
         elmnts[i_shift+3]=(idxtype)cv( 1 ) ;
      }
      else // triangles, tetrahedra, hexahedra
      {
         for( size_t i=0 ; i<cv.size() ; ++i )
         {
            elmnts[i_shift+i]=(idxtype)cv( i ) ;
         }
      }
   }
  
   // Call to METIS
   TIMER->start() ;

   // --- ancienne implémentation ---------------------------
   //
   // METIS_PartMeshDual( &ne, &nn, elmnts, &etype, &numflag, 
   //                     &nparts, &edgecut, epart, npart ) ;
   // --- nouvelle implémentation----------------------------

   int esize, esizes[] = {-1, 3, 4, 8, 4};
   esize = esizes[etype];
   idxtype* xadj   = new idxtype[ne+1] ;
   idxtype* adjncy = new idxtype[esize*(ne)] ;
   METIS_MeshToDual( &ne, &nn, elmnts, &etype, &numflag, xadj, adjncy ) ;

   int options[10] ;
   options[0] = 0 ;
   int wgtflag = 0 ;
   METIS_PartGraphKway( &ne, xadj, adjncy, NULL, NULL, &wgtflag, &numflag, 
                        &nparts, options, &edgecut, epart ) ;

   delete xadj ; xadj = 0 ;
   delete adjncy ; adjncy = 0 ;
   
   // -------------------------------------------------------

   TIMER->stop() ;

   // Fill local data structure
   RANK_CELLS.re_initialize( nb_ranks() ) ;
   RANK_CELLS.set( 0 ) ;
   for( size_t idx=0 ; idx<nb_cells() ; ++idx )
   {
      int const cell_owner = (int) epart[idx] ;
      if( cell_owner<0 || cell_owner>=(int) nb_ranks() )
      {
         EXT_METISsplitting_ERROR::n2( cell_owner, nb_ranks() ) ;
      }
      CELL_RANK( idx ) = cell_owner ;
      RANK_CELLS( (size_t) cell_owner )++ ;
   }
   for( size_t i=0 ; i<nb_ranks() ; ++i )
   {
      if( RANK_CELLS(i) == 0 ) EXT_METISsplitting_ERROR::n4( i ) ;
   }

   delete[] epart ;
   delete[] npart ;
   delete[] elmnts ;
}

//-------------------------------------------------------------------------
int
EXT_METISsplitting:: METIS_cell_type( std::string const& cell_poly_name ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_METISsplitting:: METIS_cell_type" ) ;
   PEL_CHECK( !cell_poly_name.empty() ) ;

   GE_ReferencePolyhedron const* ref_poly =
                    GE_Mpolyhedron::reference_polyhedron( cell_poly_name ) ;
   
   int result = PEL::bad_int() ;
   if( ref_poly == GE_ReferenceTriangle::object() )
   {
      result = 1 ;
   }
   else if( ref_poly == GE_ReferenceTetrahedron::object() )
   {
      result = 2 ;
   }
   else if( ref_poly == GE_ReferenceCube::object() )
   {
      result = 3 ;
   }
   else if( ref_poly == GE_ReferenceSquare::object() )
   {
      result = 4 ;
   }
   else
   {
      EXT_METISsplitting_ERROR::n3( cell_poly_name ) ;
   }

   PEL_CHECK_POST( result==1 || result==2 || result==3 || result==4 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
EXT_METISsplitting:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( GE_SplittingStrategy::invariant() ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
EXT_METISsplitting_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** EXT_METISsplitting error:\n"
      "    Invalid meshing: 2D or 3D expected for METIS splitting." ) ;
}

//internal--------------------------------------------------------------
void 
EXT_METISsplitting_ERROR:: n1( std::string const& s1, std::string const& s2 )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** EXT_METISsplitting error:\n"
      "    Invalid meshing for METIS splitting,\n"
      "    all cell polyhedra must be of same type\n"
      "          - \""+s1+"\"\n"
      "          - \""+s2+"\"\n"
      "    at least encountered." ) ;
}

//internal--------------------------------------------------------------
void 
EXT_METISsplitting_ERROR:: n2( int rk, size_t nb_rks )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** EXT_METISsplitting METIS internal error:\n"
       << "    METIS gives an invalid cell owner\n"
       << "       cell owner rank: " << rk << "\n"
       << "       number of ranks: " << nb_rks ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
EXT_METISsplitting_ERROR:: n3( std::string const& cell_poly_name )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** EXT_METISsplitting error:\n"
      "    Invalid meshing: unexpected polyhedron \""+cell_poly_name+"\"" ) ;
}

//internal--------------------------------------------------------------
void 
EXT_METISsplitting_ERROR:: n4( int rk )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** EXT_METISsplitting METIS internal error:\n"
       << "    METIS does not find cells for process of rank: " << rk << "\n" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
