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

#include <doubleVector.hh>

#include <CH_Fsaver.hh>

#include <iostream>

#include <GE_Point.hh>
#include <GE_Mpolyhedron.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <CH_BulkChemicalPotential.hh>
#include <CH_BulkEnergy.hh>

using std::cout ;
using std::endl ;

CH_Fsaver const* CH_Fsaver::PROTOTYPE = new CH_Fsaver() ;

//----------------------------------------------------------------------------
CH_Fsaver:: CH_Fsaver( void )
//----------------------------------------------------------------------------
   : PEL_Application( "CH_Fsaver" )
{
}

//----------------------------------------------------------------------------
CH_Fsaver*
CH_Fsaver:: create_replica( PEL_Object* a_owner, 
                            PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_Fsaver:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   CH_Fsaver* result = new CH_Fsaver( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
CH_Fsaver:: CH_Fsaver( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , EXP( exp->create_clone( this ) )
   , DOM( 0 )
   , cFE( 0 )
{
   PEL_ModuleExplorer* ee =
                       EXP->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   DOM = PDE_DomainAndFields::create( this, ee ) ;   
   ee->destroy() ; ee = 0 ;
   cFE = DOM->create_LocalFEcell( this ) ;
}

//----------------------------------------------------------------------------
CH_Fsaver:: ~CH_Fsaver( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
CH_Fsaver:: parse_arguments( stringVector& args ) 
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
CH_Fsaver:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_Fsaver:: run" ) ;
   
   PDE_ResultSaver* rsaver = DOM->result_saver() ;
   rsaver->start_cycle() ;
   rsaver->save_grid() ;
   rsaver->save_fields( 0 ) ;

   if( EXP->has_module( "plots" ) )
   {
      PEL_ModuleExplorer* ee = EXP->create_subexplorer( 0, "plots" ) ;
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer const* se = ee->create_subexplorer( ee ) ;
         save_associated_field( se, rsaver ) ;
      }
      ee->destroy() ; ee = 0 ;
   }
   rsaver->terminate_cycle() ;
}

//----------------------------------------------------------------------------
void
CH_Fsaver:: save_associated_field( PEL_ModuleExplorer const* exp,
                                   PDE_ResultSaver* rsaver )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_Fsaver:: save_associated_field" ) ;

   PEL_ASSERT( rsaver->grid_is_saved() ) ;

   std::string ptype = exp->string_data( "type" ) ;
   
   std::string name = exp->string_data( "name" ) ;
   PEL::out() << "*** Saving field " << name << std::endl << std::endl ;
   
   double S1 = exp->double_data( "coef_Sigma_1" ) ;
   double S2 = exp->double_data( "coef_Sigma_2" ) ;
   double S3 = exp->double_data( "coef_Sigma_3" ) ;

   CH_BulkChemicalPotential* mu = 0 ;
   CH_BulkEnergy* pp = 0 ;
   if( ptype == "CH_BulkChemicalPotential" )
   {      
      mu = CH_BulkChemicalPotential::create( 0, exp ) ;
   } 
   else if( ptype == "CH_BulkEnergy" ) 
   {
      PEL_ModuleExplorer* ee = exp->create_subexplorer( 0, "CH_BulkEnergy" ) ;
      pp = CH_BulkEnergy::make( 0, S1, S2, S3, ee ) ;
      ee->destroy() ; ee = 0 ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value(
                           exp, "type",
                           "  - \"CH_BulkChemicalPotential\"\n"
                           "  - \"CH_BulkEnergy\"" ) ;
   }

   doubleArray2D values( 0, 0 ) ;
   doubleVector default_values( 0 ) ;
   PDE_ResultSaver::SavingLocation loc = PDE_ResultSaver::AtVertices ;
   rsaver->prepare_for_field_saving( loc, name, 1, values, default_values ) ;
   
   std::string const banner ="*** CH_Fsaver : error saving \""+name+"\"" ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      GE_Mpolyhedron const* poly = cFE->polyhedron() ;
      for( size_t v=0 ; v<poly->nb_vertices() ; ++v )
      {
         GE_Point const* pt = poly->vertex( v ) ;
         size_t const iv = rsaver->vertex_index( pt ) ; 
         double c1 = pt->coordinate( 0 ) 
                   - pt->coordinate( 1 ) / PEL::sqrt(3.) ;
         double c2 = 2.* pt->coordinate( 1 ) / PEL::sqrt(3.) ;
         double c3 = 1. - c1 - c2 ;
         double val= PEL::bad_double() ;
         if( mu != 0 )
         {
            PEL_ASSERT( pp == 0 ) ;
            val = mu->F( c1, c2, c3 ) ;
         }
         else
         {
            PEL_ASSERT( pp != 0 ) ;
            val = pp->F( c1, c2, c3 ) ;
         }
         // the value at a vertex obtained from different cells
         // should be the same
         PDE_ResultSaver::check_value_consistency_at_vertex(
                               banner, pt, values( 0, iv ), val, 1.E-4, 1.E-8) ;
         values( 0, iv ) = val ;
         
         if( ( default_values( 0 ) == 
                       PDE_ResultSaver::undefined_value() ) ||
             ( val > default_values( 0 ) ) )
         {
            default_values( 0 ) = val ;
         }
         
      }
   }
   if( mu != 0 ) { mu->destroy() ; mu = 0 ; }
   if( pp != 0 ) { pp->destroy() ; pp = 0 ; }
   
   rsaver->save_field( values, default_values ) ;
}
