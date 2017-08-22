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

#include <GE_CoordinateSplitting.hh>

#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Data.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_GroupExp.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL.hh>

#include <doubleArray2D.hh>

#include <GE_Meshing.hh>

#include <string>

GE_CoordinateSplitting const* 
GE_CoordinateSplitting::PROTOTYPE = new GE_CoordinateSplitting() ;

//-------------------------------------------------------------------------
GE_CoordinateSplitting:: GE_CoordinateSplitting( void )
//------------------------------------------------------------------------
   : GE_SplittingStrategy( "GE_CoordinateSplitting" )
   , CELL_RANK( 0 )
{
}

//-------------------------------------------------------------------------
GE_CoordinateSplitting*
GE_CoordinateSplitting:: create_replica( PEL_Object* a_owner,
			                 PEL_ModuleExplorer const* exp,
			                 GE_Meshing* meshing,
                                         size_t nb_rks, size_t rk ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CoordinateSplitting:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, meshing, nb_rks, rk ) ) ;

   GE_CoordinateSplitting* result =
      new GE_CoordinateSplitting( a_owner, exp, meshing, nb_rks, rk ) ;
   
   PEL_CHECK( create_replica_POST( a_owner, exp, meshing, nb_rks, rk, 
                                   result ) ) ;
   return ( result ) ;
}

//-------------------------------------------------------------------------
GE_CoordinateSplitting:: GE_CoordinateSplitting( PEL_Object* a_owner,
			                         PEL_ModuleExplorer const* exp,
						 GE_Meshing* meshing,
                                                 size_t nb_rks, size_t rk )
//------------------------------------------------------------------------
   : GE_SplittingStrategy( a_owner, exp, meshing, nb_rks, rk )
   , CELL_RANK( meshing->nb_cells() )
{
   PEL_CHECK_INV( invariant() ) ;

   CELL_RANK.set( PEL::bad_index() ) ;

   search_owners_from_coords( meshing, exp ) ;
}

//-------------------------------------------------------------------------
GE_CoordinateSplitting*
GE_CoordinateSplitting:: create_replica( PEL_Object* a_owner,
			                 PEL_ModuleExplorer const* exp,
			                 GE_Meshing* meshing,
                                         PEL_Communicator const* com ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CoordinateSplitting:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp, meshing, com ) ) ;

   GE_CoordinateSplitting* result =
      new GE_CoordinateSplitting( a_owner, exp, meshing, com ) ;
   
   PEL_CHECK( create_replica_POST( a_owner, exp, meshing, com, 
                                   result ) ) ;
   return ( result ) ;
}

//-------------------------------------------------------------------------
GE_CoordinateSplitting:: GE_CoordinateSplitting( PEL_Object* a_owner,
			                         PEL_ModuleExplorer const* exp,
						 GE_Meshing* meshing,
                                                 PEL_Communicator const* com )
//------------------------------------------------------------------------
   : GE_SplittingStrategy( a_owner, exp, meshing, com )
   , CELL_RANK( meshing->nb_cells() )
{
   PEL_CHECK_INV( invariant() ) ;

   CELL_RANK.set( PEL::bad_index() ) ;

   search_owners_from_coords( meshing, exp ) ;
}

//-------------------------------------------------------------------------
GE_CoordinateSplitting:: ~GE_CoordinateSplitting( void )
//------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
size_t 
GE_CoordinateSplitting:: cell_rank( size_t mesh_id ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CoordinateSplitting:: cell_rank" ) ;
   PEL_CHECK_PRE( cell_rank_PRE( mesh_id ) ) ;

   size_t result = CELL_RANK( mesh_id ) ;
 
   PEL_CHECK_POST( cell_rank_POST( mesh_id, result ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
GE_CoordinateSplitting:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( GE_SplittingStrategy::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
void
GE_CoordinateSplitting:: search_owners_from_coords(
                    GE_Meshing* meshing, PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_CoordinateSplitting:: search_owners_from_coords" ) ;
   PEL_CHECK( meshing != 0 ) ;
   PEL_CHECK( meshing->nb_cells() == nb_cells() ) ;
   PEL_CHECK( exp != 0 ) ;
   
   size_t const sp_dim = meshing->nb_space_dimensions() ;

   // Set coordinates of vertices:
   doubleArray2D v_coords( meshing->nb_vertices(), sp_dim ) ;
   {
      size_t iv = 0 ;
      for( meshing->start_vertex_iterator() ;
           meshing->valid_vertex() ;
           meshing->go_next_vertex(), ++iv )
      {
         doubleVector const& coords =  meshing->vertex_coordinates() ;
         for( size_t ic=0 ; ic<sp_dim ; ++ic )
            v_coords( iv, ic ) = coords(ic) ;
      }
   }
   
 
   // Context :
   PEL_ContextSimple* ctx = PEL_ContextSimple::create( 0 ) ;
   PEL_DoubleVector* coords = PEL_DoubleVector::create( ctx, sp_dim ) ;
   ctx->extend( PEL_Variable::object( "DV_X" ), coords ) ;

   // Check formula expression:
   PEL_Data* formula =
          exp->abstract_data( 0, "coordinate_splitting_formula", ctx ) ;
   if( !formula->value_can_be_evaluated(0) )
   {
      PEL_Error::object()->raise_not_evaluable(
         exp, "coordinate_splitting_formula",
         formula->undefined_variables(0) ) ;
   }
   if( formula->data_type() != PEL_Data::Int )
   {
      PEL_Error::object()->raise_bad_data_type(
         exp, "coordinate_splitting_formula", PEL_Data::Int ) ;
   }

   // Optimized evaluation of PEL_Group expressions:
   PEL_GroupExp::set_optimized_evaluation() ;

   // Loop on cells:
   doubleVector center_coord(sp_dim) ;
   size_t idx=0 ;
   for( meshing->start_cell_iterator() ;
        meshing->valid_cell() ;
        meshing->go_next_cell(), ++idx )
   {
      size_t_vector const& vert = meshing->cell_vertices() ;
      center_coord.set(0.0) ;
      size_t const nbv = vert.size() ;
      double const inv_nbv = 1./nbv ;
      for( size_t i=0 ; i<nbv ; i++ )
      {
         for( size_t j=0 ; j<sp_dim ; j++ )
         {
            center_coord(j) += v_coords(vert(i),j)*inv_nbv ;
         }
      }
      coords->set( center_coord ) ;
      CELL_RANK( idx ) = point_owner( formula ) ;
   }
   
   PEL_GroupExp::unset_optimized_evaluation() ;
   
   ctx->destroy() ; ctx = 0 ;
   formula->destroy() ; formula = 0 ;
   coords = 0 ;
}

//----------------------------------------------------------------------
size_t
GE_CoordinateSplitting:: point_owner( PEL_Data const* formula ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_CoordinateSplitting:: point_owner" ) ;
   PEL_CHECK( formula != 0 ) ;
   PEL_CHECK( formula->data_type() == PEL_Data::Int ) ;
   PEL_CHECK( formula->value_can_be_evaluated(0) ) ;
   
   int const idx = formula->to_int() ;
   if( idx < 0 )
   {
      PEL_Error::object()->raise_plain(
         "*** GE_CoordinateSplitting error:\n"
         "    the expression of keyword \"coordinate_splitting_formula\"\n"
         "    has negative values" ) ;
   }
   if( idx >= (int)nb_ranks() )
   {
      PEL_Error::object()->raise_plain(
         "*** GE_CoordinateSplitting error:\n"
         "    the expression of keyword \"coordinate_splitting_formula\"\n"
         "    has some values equal or greater than the number of ranks" ) ;
   }
   size_t const result = idx ;

   PEL_CHECK( result<nb_ranks() ) ;
   return( result ) ;
}



