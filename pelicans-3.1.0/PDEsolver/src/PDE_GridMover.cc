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

#include <PDE_GridMover.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>
#include <doubleVector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>

#include <iostream>
#include <iomanip>
#include <sstream>

using std::endl ;

//----------------------------------------------------------------------
PDE_GridMover*
PDE_GridMover:: create( PEL_Object* a_owner,
                        PDE_DomainAndFields const* dom )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridMover:: create" ) ;
   PEL_CHECK_PRE( dom != 0 );

   PDE_GridMover* result = new PDE_GridMover( a_owner, dom ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_space_dimensions() == 
                   dom->nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_GridMover:: PDE_GridMover( PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NB_DIMS( dom->nb_space_dimensions() )
   , VERTICES( dom->set_of_vertices() )
   , cFE( dom->create_LocalFEcell( this ) )
   , DEFO( PEL_Vector::create( this, dom->set_of_vertices()->nb_points() ) )
{
   PEL_LABEL( "PDE_GridMover:: PDE_GridMover" ) ;

   cFE->include_color( GE_Color::halo_color() ) ;
   
   for( size_t i=0 ; i<VERTICES->nb_points() ; ++i )
   {
      DEFO->set_at( i, GE_Vector::create( DEFO, VERTICES->dimension() ) ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_GridMover:: ~PDE_GridMover( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridMover:: ~PDE_GridMover" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_GridMover:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridMover:: nb_space_dimensions" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( NB_DIMS ) ;
}

//----------------------------------------------------------------------
void
PDE_GridMover:: move_grid( PDE_DiscreteField const* strain,
                           size_t level,
                           double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GridMover:: move_grid" ) ;
   PEL_CHECK_PRE( strain->nb_components() == nb_space_dimensions() ) ;
   PEL_CHECK_INV( invariant() ) ;

   std::string const banner =
            "*** PDE_GridMover : error field \""+strain->name()+"\"" ;
   // Le choix des deux constantes suivantes est sujet à caution
   double const dbl_eps = 1.E-4 ;
   double const dbl_min = 1.E-8 ;
   
   size_t const nbvs = VERTICES->nb_points() ;
   for( size_t iv=0 ; iv<nbvs ; ++iv )
   {
      GE_Vector* dM = static_cast<GE_Vector*>( DEFO->at( iv ) ) ;
      for( size_t ic=0 ; ic<NB_DIMS ; ++ic )
      {
         dM->set_component( ic, PDE_ResultSaver::undefined_value() ) ;
      }
   }

   cFE->require_field_calculation( strain, PDE_LocalFE::N ) ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      GE_Mpolyhedron const* poly = cFE->polyhedron() ;

      for( size_t i=0 ; i<poly->nb_vertices() ; ++i )
      {
         GE_Point const* pt = poly->vertex( i ) ;
         PEL_CHECK( VERTICES->has( pt ) ) ;
         size_t iv = VERTICES->index( pt ) ;
         GE_Vector* dM = static_cast<GE_Vector*>( DEFO->at( iv ) ) ;

         cFE->set_calculation_point( pt ) ;
         for( size_t ic=0 ; ic<NB_DIMS ; ++ic )
         {
            double val = coef * cFE->value_at_pt( strain, level, ic ) ;

            // the value at a vertex obtained from different meshes should
            // be the same
            PDE_ResultSaver::check_value_consistency_at_vertex(
               banner, pt, dM->component( ic ), val, dbl_eps, dbl_min ) ;
            dM->set_component( ic, val ) ;
         }
      }
   }
   
   VERTICES->move_and_update( DEFO ) ;

   PEL_CHECK_INV( invariant() ) ;
}
