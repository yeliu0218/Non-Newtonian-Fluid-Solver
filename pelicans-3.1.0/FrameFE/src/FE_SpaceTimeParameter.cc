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

#include <FE_SpaceTimeParameter.hh>

#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <FE_TimeIterator.hh>

#include <FE.hh>
#include <FE_SetOfParameters.hh>

#include <iostream>
#include <sstream>

using std::ostringstream ;
using std::endl ;

FE_SpaceTimeParameter const* 
FE_SpaceTimeParameter:: PROTOTYPE = new FE_SpaceTimeParameter() ;

struct FE_SpaceTimeParameter_ERROR
{
      static void n0( FE_SpaceTimeParameter const* prm ) ;
} ;

//-------------------------------------------------------------------------
FE_SpaceTimeParameter:: FE_SpaceTimeParameter( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "FE_SpaceTimeParameter" )
{
}

//-------------------------------------------------------------------------
FE_SpaceTimeParameter*
FE_SpaceTimeParameter:: create_replica( PEL_Object* a_owner,
        	                        PDE_DomainAndFields const* dom,
                                        PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: create_replica"  ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FE_SpaceTimeParameter* result = 
                              new FE_SpaceTimeParameter( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_SpaceTimeParameter:: FE_SpaceTimeParameter( PEL_Object* a_owner,
        		                       PDE_DomainAndFields const* dom,
                                               PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , VAL( 0 )
   , GRAD( 0 )
   , COORDS( 0 )
   , TT( 0 )
   , NBCS( exp->int_data( "nb_components" ) )
{
   PEL_LABEL( "FE_SpaceTimeParameter:: FE_SpaceTimeParameter"  ) ;

   PEL_ContextSimple* ct = PEL_ContextSimple::create( this ) ;
   COORDS = PEL_DoubleVector::create( ct, doubleVector( 0 ) ) ;
   ct->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;   
   TT = PEL_Double::create( ct, 0.0 ) ;
   ct->extend( PEL_Variable::object( "DS_T" ), TT ) ;   

   VAL = exp->abstract_data( this, "value", ct ) ;
   if( !VAL->value_can_be_evaluated() )
   {
      PEL_Error::object()->raise_not_evaluable(
                    exp, "value", VAL->undefined_variables() ) ;
   }
   if( VAL->data_type() != PEL_Data::DoubleVector )
   {
      PEL_Error::object()->raise_bad_data_type(
                        exp, "value", PEL_Data::DoubleVector ) ;
   }
   
   if( exp->has_entry( "gradient" ) )
   {
      GRAD = exp->abstract_data( this, "gradient", ct ) ;
      if( !GRAD->value_can_be_evaluated() )
      {
         PEL_Error::object()->raise_not_evaluable(
                   exp, "gradient" , GRAD->undefined_variables() ) ;
      }
      if( GRAD->data_type() != PEL_Data::DoubleArray2D )
      {
         PEL_Error::object()->raise_bad_data_type(
                        exp, "gradient", PEL_Data::DoubleArray2D ) ;
      }
   }
}

//-------------------------------------------------------------------------
FE_SpaceTimeParameter:: ~FE_SpaceTimeParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
FE_SpaceTimeParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( NBCS ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   set_context( fe->calculation_point(), t_it->time() ) ;
   double result = VAL->to_double_vector()(ic) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: cell_gradient_at_pt( FE_TimeIterator const* t_it,
                                             PDE_LocalFEcell const* fe,
                                             size_t a,
                                             size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: cell_gradient_at_pt" ) ;
   PEL_CHECK_PRE( cell_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   if( GRAD == 0 ) FE_SpaceTimeParameter_ERROR::n0( this ) ;

   set_context( fe->calculation_point(), t_it->time() ) ;
   double result = GRAD->to_double_array2D()( ic, a ) ;
   
   return( result ) ;
}


//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   set_context( fe->coordinates_of_IP(), t_it->time() ) ;
   double result = VAL->to_double_vector()(ic) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                             PDE_LocalFEcell const* fe,
                                             size_t a,
                                             size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: cell_gradient_at_IP" ) ;
   PEL_CHECK_PRE( cell_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   if( GRAD == 0 ) FE_SpaceTimeParameter_ERROR::n0( this ) ;

   set_context( fe->coordinates_of_IP(), t_it->time() ) ;
   double result = GRAD->to_double_array2D()( ic, a ) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   set_context( fe->calculation_point(), t_it->time() ) ;
   double result = VAL->to_double_vector()(ic) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: bound_gradient_at_pt( FE_TimeIterator const* t_it,
                                              PDE_LocalFEbound const* fe,
                                              size_t a,
                                              size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: bound_gradient_at_pt" ) ;
   PEL_CHECK_PRE( bound_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   if( GRAD == 0 ) FE_SpaceTimeParameter_ERROR::n0( this ) ;

   set_context( fe->calculation_point(), t_it->time() ) ;
   double result = GRAD->to_double_array2D()( ic, a ) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   set_context( fe->coordinates_of_IP(), t_it->time() ) ;
   double result = VAL->to_double_vector()(ic) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: bound_gradient_at_IP( FE_TimeIterator const* t_it,
                                              PDE_LocalFEbound const* fe,
                                              size_t a,
                                              size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: bound_gradient_at_IP" ) ;
   PEL_CHECK_PRE( bound_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   if( GRAD == 0 ) FE_SpaceTimeParameter_ERROR::n0( this ) ;

   set_context( fe->coordinates_of_IP(), t_it->time() ) ;
   double result = GRAD->to_double_array2D()( ic, a ) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: side_value_at_pt( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: side_value_at_pt" ) ;
   PEL_CHECK_PRE( side_value_at_pt_PRE( t_it, fe, ic ) ) ;

   set_context( fe->calculation_point(), t_it->time() ) ;
   double result = VAL->to_double_vector()(ic) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: side_gradient_at_pt( FE_TimeIterator const* t_it,
                                             PDE_CursorFEside const* fe,
                                             size_t a,
                                             size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: side_gradient_at_pt" ) ;
   PEL_CHECK_PRE( side_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   if( GRAD == 0 ) FE_SpaceTimeParameter_ERROR::n0( this ) ;

   set_context( fe->calculation_point(), t_it->time() ) ;
   double result = GRAD->to_double_array2D()( ic, a ) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: side_value_at_IP( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: side_value_at_IP" ) ;
   PEL_CHECK_PRE( side_value_at_IP_PRE( t_it, fe, ic ) ) ;

   set_context( fe->coordinates_of_IP(), t_it->time() ) ;
   double result = VAL->to_double_vector()(ic) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SpaceTimeParameter:: side_gradient_at_IP( FE_TimeIterator const* t_it,
                                             PDE_CursorFEside const* fe,
                                             size_t a,
                                             size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: side_gradient_at_IP" ) ;
   PEL_CHECK_PRE( side_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   if( GRAD == 0 ) FE_SpaceTimeParameter_ERROR::n0( this ) ;

   set_context( fe->coordinates_of_IP(), t_it->time() ) ;
   double result = GRAD->to_double_array2D()( ic, a ) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_SpaceTimeParameter:: set_context( GE_Point const* pt, double time ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SpaceTimeParameter:: set_context" ) ;
   PEL_CHECK( pt != 0 ) ;

   TT->set( time ) ;
   COORDS->set( pt->coordinate_vector() ) ;
}

//internal------------------------------------------------------------------
void
FE_SpaceTimeParameter_ERROR:: n0( FE_SpaceTimeParameter const* prm )
//internal------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_Parameter instance of name : \"" << prm->name() << "\"" 
        << endl ;
   mesg << "   the definition MODULE has no entry of keyword \"gradient\"" 
        << endl ;
   mesg << "   (requested to compute gradients)" << endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
