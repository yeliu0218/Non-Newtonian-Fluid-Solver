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

#include <PDE_CFootFinder_TEST.hh>

#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_CFootFinder.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_PointInGridFE.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <iostream>

using std::cout ;
using std::endl ;

//---------------------------------------------------------------------------
PDE_CFootFinder_TEST*
PDE_CFootFinder_TEST:: REGISTRATOR = new PDE_CFootFinder_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_CFootFinder_TEST:: PDE_CFootFinder_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_CFootFinder", "PDE_CFootFinder_TEST" )
{
}

//---------------------------------------------------------------------------
PDE_CFootFinder_TEST:: ~PDE_CFootFinder_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_CFootFinder_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder_TEST:: process_one_test" ) ;
   
   out() << "| ... "  <<  exp->name() << endl ;
   
   PEL_ModuleExplorer* ee = 
                       exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom = PDE_DomainAndFields::create( 0, ee ) ;
   ee->destroy() ; ee = 0 ;

   PIG = PDE_PointInGridFE::create( dom, dom ) ;
   AA = dom->set_of_discrete_fields()->item( "velocity" ) ;

   cFE = dom->create_LocalFEcell( dom ) ;
   cFE_TEST = dom->create_LocalFEcell( dom ) ;
   cFE->include_color( GE_Color::halo_color() ) ; // for PDE_PointInGridFE

   ee = exp->create_subexplorer( 0, "PDE_CFootFinder" ) ;
   cFE->set_foot_finder( ee ) ;
   FINDER = cFE->foot_finder() ;
   ee->destroy() ; ee = 0 ;

   GE_Point* head = GE_Point::create( 0, dom->nb_space_dimensions() ) ;
   GE_Point* foot = GE_Point::create( 0, dom->nb_space_dimensions() ) ;
   ee = exp->create_subexplorer( 0, "searches" ) ;
   ee->start_module_iterator() ;
   for( ; ee->is_valid_module() ; ee->go_next_module() )
   {
      PEL_ModuleExplorer const* eee = ee->create_subexplorer( 0 ) ;
      head->set_coordinates( eee->doubleVector_data( "head_point" ) ) ;
      foot->set_coordinates( eee->doubleVector_data( "foot_point" ) ) ;
      double dt = eee->double_data( "time_step" ) ;
      bool foot_in = eee->bool_data( "foot_is_interior" ) ;
      
      check_search( eee->name(), dt, head, foot, foot_in ) ;
      
      eee->destroy() ; eee = 0 ;
   }
   ee->destroy() ;
   head->destroy() ;
   foot->destroy() ;
   dom->destroy() ;
}

//--------------------------------------------------------------------------
void
PDE_CFootFinder_TEST:: check_search( std::string const& test_name,
                                     double time_step,
                                     GE_Point const* head,
                                     GE_Point const* foot,
                                     bool interior_foot )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_CFootFinder_TEST:: check_search" ) ;
   
   cFE->start() ;
   PIG->is_in_grid( head, cFE ) ;
   
   FINDER->search_foot( AA, 0, time_step, head ) ;
   
   bool ok = true ;
   ok &= FINDER->foot_has_been_found() ;
   if( interior_foot )
   {
      ok &=  FINDER->foot_is_interior() ;
      ok &= !FINDER->foot_is_on_boundary() ;
   }
   else
   {
      ok &= !FINDER->foot_is_interior() ;
      ok &=  FINDER->foot_is_on_boundary() ;
   }
         
   bool ok_dist ;
   ok_dist = ( foot->distance( FINDER->foot() ) < 1.E-8 ) ;
   if( !ok_dist )
   {
      out() << "foot from jdd : " ; foot->print( out(), 0 ) ; out() << endl ;
      out() << "foot computed : " ; 
      FINDER->foot()->print( out(), 0 ) ; out() << endl ;
   }
   ok &= ok_dist ;
   
   
   cFE_TEST->go_i_th( FINDER->mesh_id_of_foot_cell() ) ;
   ok &= cFE_TEST->polyhedron()->contains( FINDER->foot() ) ;

   notify_one_test_result( test_name, ok ) ;
}

