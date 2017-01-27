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

#include <PDE_PointInGridFE_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_PointInGridFE.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_SetOfPoints.hh>

#include <doubleVector.hh>

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

//---------------------------------------------------------------------------
PDE_PointInGridFE_TEST*
PDE_PointInGridFE_TEST:: registered_test = new PDE_PointInGridFE_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_PointInGridFE_TEST:: PDE_PointInGridFE_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_PointInGridFE", "PDE_PointInGridFE_TEST" )
{
}

//---------------------------------------------------------------------------
PDE_PointInGridFE_TEST:: ~PDE_PointInGridFE_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_PointInGridFE_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PointInGridFE_TEST:: process_one_test" ) ;

   std::string const& tname = texp->name() ;

   bool verbose = texp->bool_data( "verbose" ) ;

   // Target point
   PEL_ModuleExplorer* data_exp = texp->create_subexplorer( 0, "DATA" ) ;
   GE_Point const* target_pt = 
      GE_Point::create( 0, data_exp->doubleVector_data( "target_point" ) ) ;
   if( verbose )
   {
      std::cout << " Target point : " ; 
      target_pt->print( std::cout, 1 ) ;
   }

   // Solution :
   PEL_ModuleExplorer* result_exp = texp->create_subexplorer( 0, "RESULT" ) ;
   bool in_grid_sol = result_exp->bool_data( "is_in_grid" ) ;
      
   // Domain and fields :
   PEL_ModuleExplorer* dom_exp = data_exp->create_subexplorer( 0, 
                                                    "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields const* dom = PDE_DomainAndFields::create( 0, dom_exp ) ;
   PDE_LocalFEcell* cFE = dom->create_LocalFEcell( 0 ) ;
   if( cFE->is_excluded( GE_Color::halo_color() ) )
   {
      cFE->include_color( GE_Color::halo_color() ) ;
   }
   PEL_ModuleExplorer* se = 0 ;
   if( data_exp->has_module( "PDE_PointInGridFE" ) )
   {
      se = data_exp->create_subexplorer( data_exp, "PDE_PointInGridFE" ) ;
   }
   PDE_PointInGridFE* pig = PDE_PointInGridFE::create( 0, dom, se ) ;
   
   // Direct search by loop on cells -> solution
   GE_Mpolyhedron const* sol_poly = 0 ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      if( cFE->polyhedron()->contains( target_pt ) )
      {
         sol_poly = cFE->polyhedron() ;
         break ;
      }
   }

   if( verbose )
   {
      if( sol_poly!=0 )
      {
	 std::cout << " lies in polyhedron " << cFE->mesh_id() ;
         sol_poly->print( std::cout, 1 ) ; 
      }
      else 
         std::cout << " is outside the grid " ;
      std::cout << std::endl ;
   }
   PEL_ASSERT( IMPLIES( in_grid_sol, sol_poly!=0 ) ) ;

   // Point location search :
   cFE->start() ;
   bool in_grid = pig->is_in_grid( target_pt, cFE ) ;
   if( verbose )
   {
      std::cout << " PDE_PointInGridFE found it : " ;
      if( sol_poly!=0 )
      {
	 std::cout << " in polyhedron " << cFE->mesh_id() ;
         cFE->polyhedron()->print( std::cout, 1 ) ; 
      }
      else 
         std::cout << " outside the grid " ;
      std::cout << std::endl ;
   }

   bool ok = ( in_grid && in_grid_sol ) || ( !in_grid && !in_grid_sol );
   if( ok && in_grid )
   {
      ok = ok &&  sol_poly->is_equal( cFE->polyhedron() ) ;
   }
   notify_one_test_result( tname+"::point location in grid", ok ) ;

   // Objects destruction :
   pig->destroy() ;
   cFE->destroy() ;
   dom->destroy() ;
   dom_exp->destroy() ;
   target_pt->destroy() ;
   data_exp->destroy() ;
   result_exp->destroy() ;
}
