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

#include <AP_CheckDiscretizationCFV.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <stringVector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_CursorFEside.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

using std::cout ;
using std::string ;
using std::endl ;

AP_CheckDiscretizationCFV const* 
AP_CheckDiscretizationCFV::PROTOTYPE = new AP_CheckDiscretizationCFV() ;

struct AP_CheckDiscretizationCFV_ERROR
{
   static void n0( GE_Mpolyhedron const* poly ) ;
} ;

//---------------------------------------------------------------------------
AP_CheckDiscretizationCFV:: AP_CheckDiscretizationCFV( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "AP_CheckDiscretizationCFV" )
{
}

//---------------------------------------------------------------------------
AP_CheckDiscretizationCFV*
AP_CheckDiscretizationCFV:: create_replica( PEL_Object* a_owner,
  			         	  PDE_DomainAndFields const* dom,
				 	  FE_SetOfParameters const* prms,
				          PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_CheckDiscretizationCFV:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   AP_CheckDiscretizationCFV* result = 
                      new AP_CheckDiscretizationCFV( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_CheckDiscretizationCFV:: AP_CheckDiscretizationCFV( 
                                               PEL_Object* a_owner,
					       PDE_DomainAndFields const* dom,
					       FE_SetOfParameters const* prms,
					       PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , BCs( dom->set_of_boundary_conditions() )
   , MAX_NEG_DISTANCE_TO_SIDE( 
         exp->double_data( "max_allowed_negative_distance_to_face" ) )
   , MAX_SCALAR_PRODUCT( 
         exp->double_data( "max_allowed_normal_scalar_VtoFVcenter" ) )
   , MIN_DISTANCE_CENTERS( 
         exp->double_data( "min_distance_between_centers" ) )
{
   PEL_LABEL( "AP_CheckDiscretizationCFV:: AP_CheckDiscretizationCFV" ) ;

   PEL_ASSERT( MAX_NEG_DISTANCE_TO_SIDE < 0.0 ) ;
   PEL_ASSERT( MAX_SCALAR_PRODUCT > 0.0 ) ;

   {
      stringVector const& nns = exp->stringVector_data( "fields" ) ;
      for( size_t i=0 ; i<nns.size() ; ++i )
      {
         PDE_DiscreteField const* ff = 
                           dom->set_of_discrete_fields()->item( nns( i ) ) ;
         FIELDS.push_back( ff ) ;
         cFE->require_field_calculation( ff, PDE_LocalFE::node ) ;
      }
   }
   {
      stringVector const& nns = exp->stringVector_data( "fields_with_BCs" ) ;
      for( size_t i=0 ; i<nns.size() ; ++i )
      {
         PDE_DiscreteField const* ff = 
                           dom->set_of_discrete_fields()->item( nns( i ) ) ;
         BCFIELDS.push_back( ff ) ;
      }
   }
}

//---------------------------------------------------------------------------
AP_CheckDiscretizationCFV:: ~AP_CheckDiscretizationCFV( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
AP_CheckDiscretizationCFV:: do_before_time_stepping( 
                                                FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_CheckDiscretizationCFV:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;

   start_total_timer( "AP_CheckDiscretizationCFV:: do_before_time_stepping" ) ;
   // -------------

   bool ok ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      GE_Mpolyhedron const* poly = cFE->polyhedron() ;

      GE_Point const* cc = poly->finite_volume_center() ;
      if( cc == 0 ) AP_CheckDiscretizationCFV_ERROR::n0( poly ) ;
      ok = poly->contains( cc ) ;
      if( !ok ) display_not_in( poly ) ;

      std::vector< PDE_DiscreteField const* >::const_iterator it = 
                                                              FIELDS.begin() ;
      for( ; it != FIELDS.end() ; ++it )
      {
         check_nb_local_nodes( *it, cFE, 1 ) ;
      }
   }

   PDE_LocalFEcell const* cFE0 = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* cFE1 = sFE->adjacent_localFEcell( 1 ) ;
   
   GE_Vector* cc = GE_Vector::create( 0, cFE->nb_space_dimensions() ) ;
   GE_Vector* vv = GE_Vector::create( 0, cFE->nb_space_dimensions() ) ;

   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      GE_Point const* c0 = cFE0->polyhedron()->finite_volume_center() ;
      GE_Point const* c1 = cFE1->polyhedron()->finite_volume_center() ;
      PEL_ASSERT( c0!=0 && c1!=0 ) ;

      double dd0 = sFE->distance_to_adjacent_finite_volume_center( 0 ) ;
      ok = ( dd0 > MAX_NEG_DISTANCE_TO_SIDE ) ;

      double dd1 = sFE->distance_to_adjacent_finite_volume_center( 1 ) ;
      ok = ok && ( dd1 > MAX_NEG_DISTANCE_TO_SIDE ) ;

      double dd = dd0 + dd1 ;
      ok = ok && ( dd > MIN_DISTANCE_CENTERS ) ;
      
      bool ok_n = true ;
      if( !sFE->is_periodic() )
      {
         cc->re_initialize( c0, c1 ) ;
         GE_Mpolyhedron const* s_poly = sFE->polyhedron() ;
         for( size_t iv=1 ; iv<s_poly->nb_vertices() ; ++iv )
         {
            vv->re_initialize( s_poly->vertex( 0 ), s_poly->vertex( iv ) ) ;
            double xx = cc->dot_product( vv ) ;
            ok_n = ok_n && ( PEL::abs( xx ) < MAX_SCALAR_PRODUCT ) ;
         }
      }

      if( !ok || !ok_n ) display_side_pb( sFE, dd0, dd1, dd, ok_n ) ;
   }

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      std::vector< PDE_DiscreteField const* >::const_iterator it = 
                                                            BCFIELDS.begin() ;
      for( ; it != BCFIELDS.end() ; ++it )
      {
         check_boundary_condition( BCs, bFE->color(), *it ) ;
      }

      double dd = bFE->distance_to_adjacent_finite_volume_center() ;
      ok = ( dd > MAX_NEG_DISTANCE_TO_SIDE ) ;
      if( !ok ) display_bound_pb( bFE, dd ) ;
   }

   cc->destroy() ;
   vv->destroy() ;

   stop_total_timer() ;
   // -------------
}

//---------------------------------------------------------------------------
void 
AP_CheckDiscretizationCFV:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_CheckDiscretizationCFV:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
}

//---------------------------------------------------------------------------
void
AP_CheckDiscretizationCFV:: display_side_pb( PDE_CursorFEside const* fe,
                                           double dd0, 
                                           double dd1,
                                           double dd,
                                           bool ok_n ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_CheckDiscretizationCFV:: display_side_pb" ) ;

   std::ostringstream msg ;
   msg << "--------------------------" << endl ;
   msg << "the side : " << endl ;
   fe->polyhedron()->print( msg , 6 ) ;
   msg << "has improper finite volume surroundings:" << endl ;
   msg << "   distance to center 0 : " << dd0 << endl ;
   msg << "   distance to center 1 : " << dd1 << endl ;
   msg << "   distance between centers : " << dd << endl ;
   if( !ok_n )
      msg << "   the orthogonality property is not satisfied" << endl ;
   msg << "--------------------------" << endl ;
   PEL_Error::object()->display_info( msg.str() ) ;
}

//---------------------------------------------------------------------------
void
AP_CheckDiscretizationCFV:: display_bound_pb( PDE_LocalFEbound const* fe,
                                            double dd ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_CheckDiscretizationCFV:: display_bound_pb" ) ;

   std::ostringstream msg ;
   msg << "--------------------------" << endl ;
   msg << "the bound : " << endl ;
   fe->polyhedron()->print( msg , 6 ) ;
   msg << "has improper finite volume surroundings:" << endl ;
   msg << "   distance to center 0 : " << dd << endl ;
   PEL_Error::object()->display_info( msg.str() ) ;
}

//---------------------------------------------------------------------------
void
AP_CheckDiscretizationCFV:: display_not_in( GE_Mpolyhedron const* poly ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_CheckDiscretizationCFV:: display_not_in" ) ;

   std::ostringstream msg ;
   msg << "--------------------------" << endl ;
   msg << "the cell : " << endl ;
   poly->print( msg , 6 ) ;
   msg << "does not contain its finite volume center:" << endl ;
   poly->finite_volume_center()->print( msg , 6 ) ;
   msg << endl << "--------------------------" << endl ;
   PEL_Error::object()->display_info( msg.str() ) ;
}

//internal------------------------------------------------------------------
void
AP_CheckDiscretizationCFV_ERROR:: n0( GE_Mpolyhedron const* poly )
//internal------------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** AP_CheckDiscretizationCFV:" << endl << endl ;
   msg << "    the \"finite volume center\" has not been defined" << endl ;
   msg << "    for the following polyhedron:" << endl << endl ;
   poly->print( msg, 4 ) ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
