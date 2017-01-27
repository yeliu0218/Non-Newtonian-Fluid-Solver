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

#include <CH_InterfaceIndicator.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <iostream>
#include <string>

using std::string ;
using std::cout ; using std::endl ;

CH_InterfaceIndicator const*
CH_InterfaceIndicator:: PROTOTYPE = new CH_InterfaceIndicator() ;

//---------------------------------------------------------------------------
CH_InterfaceIndicator:: CH_InterfaceIndicator( void )
//---------------------------------------------------------------------------
   : PDE_AdaptationIndicator( "CH_InterfaceIndicator" )
   , CELL_INTERF( 0 )
{
}

//---------------------------------------------------------------------------
CH_InterfaceIndicator*
CH_InterfaceIndicator:: create_replica( 
                                  PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  PEL_ModuleExplorer const* exp,
                                  size_t a_verbose_level ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_InterfaceIndicator:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp, a_verbose_level ) ) ;

   CH_InterfaceIndicator* result = 
       new CH_InterfaceIndicator( a_owner, dom, exp, a_verbose_level ) ;

   PEL_CHECK( create_replica_POST( result,
                                   a_owner, dom, exp, a_verbose_level ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
CH_InterfaceIndicator:: CH_InterfaceIndicator( 
                                           PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           PEL_ModuleExplorer const* exp,
                                           size_t a_verbose_level  )
//---------------------------------------------------------------------------
   : PDE_AdaptationIndicator( a_owner, a_verbose_level )
   , C1( dom->set_of_discrete_fields()->item( 
                              exp->string_data( "phase_field_1" ) ) )
   , C2( 0 )
   , LL( exp->int_data( "level" ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , QRP( GE_QRprovider::object( 
                           exp->string_data( "quadrature_rule_provider" ) ) )
   , CELL_INTERF( 0 )
   , H_INTERF( exp->double_data( "h_for_interface" ) )
   , VAL_REFI( exp->double_data( "refinement_limit" ) )
   , VAL_UNREFI( exp->double_data( "unrefinement_limit" ) )
{
   PEL_ASSERT( C1->nb_components() == 1 ) ; //??????????

   cFE->include_color( GE_Color::halo_color() ) ;

   cFE->require_field_calculation( C1, PDE_LocalFE::N ) ;

   if( exp->has_entry( "phase_field_2" ) )
   {
      C2 = dom->set_of_discrete_fields()->item( 
                                exp->string_data( "phase_field_2" ) ) ;
      PEL_ASSERT( C2->nb_components() == 1 ) ; //??????????
      cFE->require_field_calculation( C2, PDE_LocalFE::N ) ;
   }

   CELL_INTERF.re_initialize( cFE->nb_meshes() ) ;
}

//---------------------------------------------------------------------------
CH_InterfaceIndicator:: ~CH_InterfaceIndicator( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
CH_InterfaceIndicator:: reset( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
CH_InterfaceIndicator:: build( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_InterfaceIndicator:: build" ) ;

   CELL_INTERF.set( 0.0 ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   { 
      size_t id = cFE->mesh_id() ;
      if( id >= CELL_INTERF.size() )
      {
         CELL_INTERF.resize( id+1 ) ;
      }
      cFE->start_IP_iterator( QRP ) ;
      double ss_1 = 0.0 ;
      double ss_2 = 0.0 ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         ss_1 += cFE->value_at_IP( C1, LL, 0 ) * cFE->weight_of_IP() ;
         if( C2 != 0 )
            ss_2 += cFE->value_at_IP( C2, LL, 0 ) * cFE->weight_of_IP() ;
      }
      double meas = cFE->polyhedron()->measure() ;
      CELL_INTERF( id ) = 
           PEL::max( meas-ss_1-ss_2, PEL::max( ss_1, ss_2 ) ) / meas ;
   }
}

//---------------------------------------------------------------------------
bool
CH_InterfaceIndicator:: to_be_refined( double bf_indicator,
                                       GE_Mpolyhedron const* poly,
                                       PDE_ReferenceElement const* elm,
                                       size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_InterfaceIndicator:: to_be_refined" ) ;
   PEL_CHECK_PRE( to_be_refined_PRE( bf_indicator, poly, elm, local_node ) ) ;

   bool result = false ;

   if( bf_indicator < VAL_REFI )
   {
      double h = poly->inter_vertices_maximum_distance() ;
      if( h > H_INTERF ) result = true ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
bool
CH_InterfaceIndicator:: to_be_unrefined( double bf_indicator,
                                         GE_Mpolyhedron const* poly,
                                         PDE_ReferenceElement const* elm,
                                         size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_InterfaceIndicator:: to_be_unrefined" ) ;
   PEL_CHECK_PRE( to_be_unrefined_PRE( bf_indicator, poly, elm, local_node ) );

   bool result = false ;

   if( bf_indicator > VAL_UNREFI )
   {
      result = true ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------------
double
CH_InterfaceIndicator:: cell_indicator( size_t cell_id ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_InterfaceIndicator:: cell_indicator" ) ;

   cFE->go_i_th( cell_id ) ;
   PEL_ASSERT( cFE->is_valid() ) ;

   double result = CELL_INTERF( cell_id ) ;
   return( result ) ;
}
