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

#include <FE_BCupdate.hh>

#include <PEL_ContextSimple.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <FE_TimeIterator.hh>

#include <iostream>
#include <cmath>

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

FE_BCupdate const* FE_BCupdate::PROTOTYPE = new FE_BCupdate() ;

//---------------------------------------------------------------------------
FE_BCupdate:: FE_BCupdate( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "FE_BCupdate" )
   , L_FFs( 0 )
{
}

//---------------------------------------------------------------------------
FE_BCupdate*
FE_BCupdate:: create_replica( PEL_Object* a_owner,
                              PDE_DomainAndFields const* dom,
                              FE_SetOfParameters const* prms,
                              PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_BCupdate:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_BCupdate* result = new FE_BCupdate( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_BCupdate:: FE_BCupdate( PEL_Object* a_owner,
                           PDE_DomainAndFields const* dom,
                           FE_SetOfParameters const* prms,
                           PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , FF( dom->set_of_discrete_fields()->item( exp->string_data( "field" ) ) )
   , L_FFs( exp->intVector_data( "levels_of_field" ) )
   , BCs( dom->set_of_boundary_conditions() )
   , bFE( dom->create_LocalFEbound( this ) )
   , CONTEXT( PEL_ContextSimple::create( this ) )
   , XX( 0 )
   , TT( 0 )
{
   PEL_LABEL( "FE_BCupdate:: FE_BCupdate" ) ;

   for( size_t l=0 ; l<L_FFs.size() ; ++l )
      check_field_storage_depth( FF, L_FFs( l ) ) ;
            
   bFE->require_field_calculation( FF, PDE_LocalFE::N ) ;
 
   XX = PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) ;
   CONTEXT->extend( PEL_Variable::object("DV_X"), XX ) ;

   TT = PEL_Double::create( CONTEXT, 0.0 ) ;
   CONTEXT->extend( PEL_Variable::object("DS_T"), TT ) ;

}

//---------------------------------------------------------------------------
FE_BCupdate:: ~FE_BCupdate( void )
//---------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void
FE_BCupdate:: do_before_inner_iterations_stage( FE_TimeIterator const* t_it )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_BCupdate:: do_before_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "FE_BCupdate:: do_before_inner_iterations_stage" ) ;
   
   TT->set( t_it->time() ) ;

   size_t glob_n = 0 ;

   boolVector done( FF->nb_nodes() ) ;
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   { 
      GE_Color const* color = bFE->color() ;
      if( BCs->has_BC( color, FF ) )
      {
         PEL_ModuleExplorer const* ee =  BCs->BC_explorer( color, FF ) ;
         if( ee->string_data( "type" ) == "Dirichlet_to_update" )
         {
            for( size_t loc_n=0 ; loc_n<bFE->nb_local_nodes( FF ) ; ++loc_n )
            {
               if( bFE->local_node_is_in_mesh( FF, loc_n ) )
               {
                  GE_Point const* pt = bFE->local_node_location( FF, loc_n ) ;
                  glob_n = bFE->global_node( FF, loc_n ) ;
                  if( ! done( glob_n ) )
                  {
                     done( glob_n ) = true ;
                     XX->set( pt->coordinate_vector() ) ;
                     doubleVector const& val = 
                                  ee->doubleVector_data( "value", CONTEXT ) ;
                     for( size_t i=0 ; i<FF->nb_components() ; ++i )
                     {
                        if( FF->DOF_has_imposed_value( glob_n, i ) )
                        {
                           for( size_t l=0 ; l<L_FFs.size() ; ++l )
                           {
                              size_t level = L_FFs( l ) ;
                              FF->set_DOF_value( level, 
                                                 glob_n, val( i ), i ) ;
                           }
                           FF->set_DOF_imposed_value( glob_n, val(i), i ) ;
                        }
                     }
                  }
               }
            }
         }	 
      }
   }

   stop_total_timer() ;

   PEL_CHECK_POST( do_before_inner_iterations_stage_POST( t_it ) ) ;
}

//---------------------------------------------------------------------------
void
FE_BCupdate:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_BCupdate:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
}
