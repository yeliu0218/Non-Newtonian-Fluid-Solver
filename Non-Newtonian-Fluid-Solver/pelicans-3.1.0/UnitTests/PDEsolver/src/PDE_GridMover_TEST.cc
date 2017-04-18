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

#include <PDE_GridMover_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_GridMover.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <GE_Point.hh>
#include <GE_SetOfPoints.hh>

#include <doubleVector.hh>

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

//---------------------------------------------------------------------------
PDE_GridMover_TEST*
PDE_GridMover_TEST:: registered_test = new PDE_GridMover_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_GridMover_TEST:: PDE_GridMover_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_GridMover", "PDE_GridMover_TEST" )
{
}

//---------------------------------------------------------------------------
PDE_GridMover_TEST:: ~PDE_GridMover_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_GridMover_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridMover_TEST:: process_one_test" ) ;
   
   std::string const& tname = texp->name() ;

   PEL_ModuleExplorer* data_exp = texp->create_subexplorer( 0, "DATA" ) ;
   PEL_ModuleExplorer* result_exp = texp->create_subexplorer( 0, "RESULT" ) ;
   
   // Domain and fields :
   PDE_DomainAndFields const* dom = PDE_DomainAndFields::create( 0, data_exp );

   // Fields :
   PDE_SetOfDiscreteFields const* fields = dom->set_of_discrete_fields() ;

   // Grid motion field :
   std::string const grid_motion_field_name = "GridDeformation" ;
   if( !fields->has( grid_motion_field_name ) )
   {
      PEL_Error::object()->raise_plain(
         "The grid deformation field \""+grid_motion_field_name+"\" is expected." ) ;
   }
   PDE_DiscreteField const* grid_motion_field =
                                   fields->item( grid_motion_field_name ) ;
   if( grid_motion_field->storage_depth()!=1 )
   {
      PEL_Error::object()->raise_plain(
         "The grid deformation field \""+grid_motion_field_name+"\" should have only one storage level." ) ;
   }
   if( grid_motion_field->nb_components()!=dom->nb_space_dimensions() )
   {
      PEL_Error::object()->raise_plain(
         "Bad dimension for the grid deformation field \""+grid_motion_field_name+"\"." ) ;
   }

   // Grid mover :
   PDE_GridMover* grid_mover = PDE_GridMover::create( 0, dom ) ;
   grid_mover->move_grid( grid_motion_field, 0 ) ;
   grid_mover->destroy() ; grid_mover = 0 ;

   // Check result :
   GE_SetOfPoints const* vertices = dom->set_of_vertices() ;
   doubleVector const& verts =
                result_exp->doubleVector_data( "vertices_coordinates" ) ;

   bool ok = verts.size()==vertices->nb_points()*vertices->dimension() ;
   if( !ok )
   {
      std::ostringstream mess ;
      mess << std::setprecision(9)
           << std::setiosflags( std::ios::scientific ) ;
      mess << "vertices_coordinates : " << std::endl ;
      for( size_t i=0 ; i<verts.size() ; ++i )
      {
         mess << "   " << verts(i) << std::endl ;
      }
      mess << "moved vertices : " << std::endl ;
      vertices->print( mess, 3 ) ;
      mess << std::endl << std::endl ;
      mess << "Bad dimension for table \"vertices_coordinates\"." ;
      PEL_Error::object()->raise_plain( mess.str() ) ;
   }

   if( ok )
   {
      for( size_t i=0 ; ok && i<vertices->nb_points() ; ++i )
      {
         GE_Point const* pt = vertices->point(i) ;
         for( size_t j=0 ; ok && j<vertices->dimension() ; ++j )
         {
            double const x = verts(vertices->dimension()*i+j) ;
            ok = PEL::toler( pt->coordinate(j)-x, 1.E-8 ) ;
         }
      }
      if( !ok )
      {
         std::ostringstream mess ;
         mess << "vertices_coordinates : " << std::endl ;
         mess << std::setprecision(9)
              << std::setiosflags( std::ios::scientific ) ;
         for( size_t i=0 ; i<vertices->nb_points() ; ++i )
         {
            for( size_t j=0 ; j<vertices->dimension() ; ++j )
            {
               double const x = verts(vertices->dimension()*i+j) ;
               mess << "   " << x ;
            }
            mess << std::endl ;
         }
         mess << "moved vertices : " << std::endl ;
         vertices->print( mess, 3 ) ;
         PEL_Error::object()->raise_plain( mess.str() ) ;
      }
   }
   
   notify_one_test_result( tname+"::vertices_coordinates", ok ) ;

   // Objects destruction :
   dom->destroy() ;
   data_exp->destroy() ;
   result_exp->destroy() ;
}
