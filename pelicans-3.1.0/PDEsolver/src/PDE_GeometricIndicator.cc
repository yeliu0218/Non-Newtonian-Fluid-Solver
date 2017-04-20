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

#include <PDE_GeometricIndicator.hh>

#include <PEL.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>
#include <GE_Vector.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_SetOfBCs.hh>

#include <iostream>
#include <string>

using std::string ;
using std::cout ; using std::endl ;

PDE_GeometricIndicator const*
PDE_GeometricIndicator:: PROTOTYPE = new PDE_GeometricIndicator() ;

//---------------------------------------------------------------------------
PDE_GeometricIndicator:: PDE_GeometricIndicator( void )
//---------------------------------------------------------------------------
   : PDE_AdaptationIndicator( "PDE_GeometricIndicator" )
{
}

//---------------------------------------------------------------------------
PDE_GeometricIndicator*
PDE_GeometricIndicator:: create_replica( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         PEL_ModuleExplorer const* exp,
                                         size_t a_verbose_level  ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricIndicator:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp, a_verbose_level ) ) ;

   PDE_GeometricIndicator* result = 
       new PDE_GeometricIndicator( a_owner, dom, exp, a_verbose_level ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, 
                                        dom, exp, a_verbose_level ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_GeometricIndicator:: PDE_GeometricIndicator( 
                                           PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           PEL_ModuleExplorer const* exp,
                                           size_t a_verbose_level  )
//---------------------------------------------------------------------------
   : PDE_AdaptationIndicator( a_owner, a_verbose_level )
   , CTX( 0 )
   , COORDS( 0 )
   , ITER( 0 )
   , R_INDIC( 0 )
   , U_INDIC( 0 )
   , NB_REFS( exp->int_data( "nb_steps" ) )
   , PT( GE_Point::create( this, dom->nb_space_dimensions() ) )
   , ICALL( 0 )
   , ICALL_MAX( 0 )
{
   ICALL_MAX = NB_REFS ;
   
   PEL_ContextSimple* c = PEL_ContextSimple::create( this ) ;
   COORDS = PEL_DoubleVector::create( c, doubleVector(0) ) ;
   c->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
   ITER = PEL_Int::create( c, 0 ) ;
   c->extend( PEL_Variable::object( "IS_ITER" ), ITER ) ;
   CTX = c ;
   R_INDIC = exp->abstract_data( this, "refinement_indicator", CTX ) ;
   if( !R_INDIC->value_can_be_evaluated() )
   {
      PEL_Error::object()->raise_not_evaluable(
         exp, "refinement_indicator", R_INDIC->undefined_variables() ) ;
   }
   if( R_INDIC->data_type()!=PEL_Data::Bool )
   {
      PEL_Error::object()->raise_bad_data_type(
            exp, "refinement_indicator", PEL_Data::Bool ) ;
   }
   if( exp->has_entry( "unrefinement_indicator" ) )
   {
      U_INDIC = exp->abstract_data( this, "unrefinement_indicator", CTX ) ;
      if( !U_INDIC->value_can_be_evaluated() )
      {
         PEL_Error::object()->raise_not_evaluable(
            exp, "unrefinement_indicator", U_INDIC->undefined_variables() ) ;
      }
      if( U_INDIC->data_type()!=PEL_Data::Bool )
      {
         PEL_Error::object()->raise_bad_data_type(
            exp, "unrefinement_indicator", PEL_Data::Bool ) ;
      }
   }
}

//---------------------------------------------------------------------------
PDE_GeometricIndicator:: ~PDE_GeometricIndicator( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_GeometricIndicator:: reset( void )
//---------------------------------------------------------------------------
{
   ICALL_MAX = ICALL + NB_REFS ;
}

//---------------------------------------------------------------------------
void
PDE_GeometricIndicator:: build( void )
//---------------------------------------------------------------------------
{
   ++ICALL ;
}

//---------------------------------------------------------------------------
bool
PDE_GeometricIndicator:: to_be_refined( double bf_indicator,
                                         GE_Mpolyhedron const* poly,
                                         PDE_ReferenceElement const* elm,
                                         size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricIndicator:: to_be_refined" ) ;
   PEL_CHECK_PRE( to_be_refined_PRE( bf_indicator, poly, elm, local_node ) ) ;

   bool result = false ;
   if( ICALL <= ICALL_MAX )
   {
      poly->apply_mapping( elm->node_location( local_node ), PT ) ;
      COORDS->set( PT->coordinate_vector() ) ;
      ITER->set( ICALL-1 ) ;
      result = ( R_INDIC->to_bool() ) ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
bool
PDE_GeometricIndicator:: to_be_unrefined( double bf_indicator,
                                           GE_Mpolyhedron const* poly,
                                           PDE_ReferenceElement const* elm,
                                           size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricIndicator:: to_be_unrefined" ) ;
   PEL_CHECK_PRE( to_be_unrefined_PRE( bf_indicator, poly, elm, local_node ) );

   bool result = false ;
   if( U_INDIC != 0 ) 
   {
      if( ICALL <= ICALL_MAX )
      {
         poly->apply_mapping( elm->node_location( local_node ), PT ) ;
         COORDS->set( PT->coordinate_vector() ) ;
         ITER->set( ICALL-1 ) ;
         result = ( U_INDIC->to_bool() ) ;
      }
   }

   return( result ) ;
}

//----------------------------------------------------------------------------
double
PDE_GeometricIndicator:: cell_indicator( size_t cell_id ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_GeometricIndicator:: cell_indicator" ) ;

   double result = 0.0 ;
   return( result ) ;
}
