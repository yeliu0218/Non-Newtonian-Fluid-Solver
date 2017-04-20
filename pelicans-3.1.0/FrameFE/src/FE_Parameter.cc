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

#include <FE_Parameter.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <GE_Mpolyhedron.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDomains.hh>

#include <iostream>
#include <sstream>

using std::string ;
using std::endl ;
using std::ostringstream ;
using std::map ;

//-------------------------------------------------------------------------
FE_Parameter*
FE_Parameter:: make( PEL_Object* a_owner,
                     PDE_DomainAndFields const* dom,
		     PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: make" ) ;
   PEL_CHECK_PRE( dom != 0 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   string name = exp->string_data( "concrete_name" ) ;
   FE_Parameter const* proto =
      static_cast<FE_Parameter const*>( plugins_map()->item( name ) ) ;
      
   FE_Parameter* result = proto->create_replica( a_owner, dom, exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_Parameter*
FE_Parameter:: make( PEL_Object* a_owner,
                     PDE_SetOfDomains const* sdoms,
		     PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: make" ) ;
   PEL_CHECK_PRE( sdoms != 0 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   string name = exp->string_data( "concrete_name" ) ;
   FE_Parameter const* proto =
      static_cast<FE_Parameter const*>( plugins_map()->item( name ) ) ;
      
   FE_Parameter* result = proto->create_replica( a_owner, sdoms, exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_Parameter:: FE_Parameter( PEL_Object* a_owner,
                             std::string const& a_name )
//-------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , NAME( a_name )
{
   PEL_CHECK_POST( !is_a_prototype() ) ;
}

//-------------------------------------------------------------------------
FE_Parameter*
FE_Parameter:: create_replica( PEL_Object* a_owner,
                               PDE_SetOfDomains const* sdoms,
                               PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( a_owner, sdoms, exp ) ) ;

   FE_Parameter* result = 0 ;

   if( exp->has_entry( "domain" ) )
   {
      string const& dom_name = exp->string_data( "domain" ) ;
      PDE_DomainAndFields const* dom = sdoms->domain( dom_name ) ;
      result = create_replica( a_owner, dom, exp ) ;
   }
   else
   {
      ostringstream mesg ;
      mesg << type_name() << " : " << endl ;
      mesg << "   an entry of keyword \"domain\" is required" << endl ;
      mesg << "   if there is more than one domain" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( create_replica_POST( result, a_owner, sdoms, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_Parameter:: FE_Parameter( std::string const& a_name )
//-------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
{
   PEL_LABEL( "FE_Parameter:: FE_Parameter" ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//-------------------------------------------------------------------------
FE_Parameter:: ~FE_Parameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
FE_Parameter:: reset( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: reset" ) ;
   PEL_CHECK_PRE( reset_PRE( exp ) ) ;
   
   if( exp != 0 )
   {
      std::string msg =
         "*** FE_Parameter<" + name() + "> error\n"
         "    reset with explorer not implemented" ;
      PEL_Error::object()->raise_plain( msg ) ;
   }
}

//-------------------------------------------------------------------------
std::string const&
FE_Parameter:: name( void ) const
//-------------------------------------------------------------------------
{
   return( NAME ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: transfer_cell_calculation_requirements( 
                                                 PDE_LocalFEcell* fe,
                                                 int requisite )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: transfer_cell_calculation_requirements" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( !fe->is_valid() ) ;
   PEL_CHECK_PRE( requisite == FE_Parameter::Val || 
                  requisite == FE_Parameter::Grad ) ;

   if( requisite & FE_Parameter::Val )
   {
      prepare_for_value_on_cells( fe ) ;
      VAL_ON_CELLS[ fe ]  = true ;
   }
   if( requisite & FE_Parameter::Grad )
   {
      prepare_for_gradient_on_cells( fe ) ;
      GRAD_ON_CELLS[ fe ] = true ;
   }
	 
   PEL_CHECK_POST( ok_for_cell_calculations( fe, requisite ) ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: ok_for_cell_calculations( PDE_LocalFEcell const* fe,
                                         int requisite ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: ok_for_cell_calculations" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( requisite == FE_Parameter::Val || 
                  requisite == FE_Parameter::Grad ) ;
   
   bool result = false ;
   if( requisite & FE_Parameter::Val )
   {
      map< PDE_LocalFEcell const*, bool >::const_iterator it = 
                                                     VAL_ON_CELLS.find( fe ) ;
      if( it != VAL_ON_CELLS.end() ) result = it->second ;
   }
   if( requisite & FE_Parameter::Grad )
   {
      map< PDE_LocalFEcell const*, bool >::const_iterator it = 
                                                     GRAD_ON_CELLS.find( fe ) ;
      if( it != GRAD_ON_CELLS.end() ) result = it->second ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: transfer_side_calculation_requirements( 
                                                 PDE_CursorFEside* fe,
                                                 int requisite )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: transfer_side_calculation_requirements" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( !fe->is_valid() ) ;
   PEL_CHECK_PRE( requisite == FE_Parameter::Val || 
                  requisite == FE_Parameter::Grad ) ;

   if( requisite & FE_Parameter::Val )
   {
      prepare_for_value_on_sides( fe ) ;
      VAL_ON_CELLS[ fe->adjacent_localFEcell(0) ]  = true ;
      VAL_ON_CELLS[ fe->adjacent_localFEcell(1) ]  = true ;
      
   }
   if( requisite & FE_Parameter::Grad )
   {
      prepare_for_gradient_on_sides( fe ) ;
      GRAD_ON_CELLS[ fe->adjacent_localFEcell(0) ]  = true ;
      GRAD_ON_CELLS[ fe->adjacent_localFEcell(1) ]  = true ;
   }
	 
   PEL_CHECK_POST( ok_for_side_calculations( fe, requisite ) ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: ok_for_side_calculations( PDE_CursorFEside const* fe,
                                         int requisite ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: ok_for_side_calculations" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( requisite == FE_Parameter::Val || 
                  requisite == FE_Parameter::Grad ) ;
   
   bool result = ok_for_cell_calculations( fe->adjacent_localFEcell(0), requisite ) ;
   PEL_CHECK(
      EQUIVALENT(
         result,
         ok_for_cell_calculations( fe->adjacent_localFEcell(1), requisite ) ) ) ;

   PEL_CHECK_POST( 
      EQUIVALENT( result,
                  ok_for_cell_calculations( fe->adjacent_localFEcell(0),
                                            requisite ) ) ) ;
   PEL_CHECK_POST( 
      EQUIVALENT( result,
                  ok_for_cell_calculations( fe->adjacent_localFEcell(1),
                                            requisite ) ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: transfer_bound_calculation_requirements( 
                                                 PDE_LocalFEbound* fe,
                                                 int requisite )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: transfer_bound_calculation_requirements" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( !fe->is_valid() ) ;
   PEL_CHECK_PRE( requisite == FE_Parameter::Val || 
                  requisite == FE_Parameter::Grad ) ;

   if( requisite & FE_Parameter::Val )
   {
      prepare_for_value_on_bounds( fe ) ;
      VAL_ON_BOUNDS[ fe ]  = true ;
   }
   if( requisite & FE_Parameter::Grad )
   {
      prepare_for_gradient_on_bounds( fe ) ;
      GRAD_ON_BOUNDS[ fe ] = true ;
   }
	 
   PEL_CHECK_POST( ok_for_bound_calculations( fe, requisite ) ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: ok_for_bound_calculations( PDE_LocalFEbound const* fe,
                                          int requisite ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: ok_for_bound_calculations" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( requisite == FE_Parameter::Val || 
                  requisite == FE_Parameter::Grad ) ;
   
   bool result = false ;
   if( requisite & FE_Parameter::Val )
   {
      map< PDE_LocalFEbound const*, bool >::const_iterator it = 
                                                   VAL_ON_BOUNDS.find( fe ) ;
      if( it != VAL_ON_BOUNDS.end() ) result = it->second ;
   }
   if( requisite & FE_Parameter::Grad )
   {
      map< PDE_LocalFEbound const*, bool >::const_iterator it =
                                                    GRAD_ON_BOUNDS.find( fe ) ;
      if( it != GRAD_ON_BOUNDS.end() ) result = it->second ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK_PRE( prepare_for_value_on_cells_PRE( fe ) ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: prepare_for_gradient_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: prepare_for_gradient_on_cells" ) ;
   PEL_CHECK_PRE( prepare_for_gradient_on_cells_PRE( fe ) ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: prepare_for_value_on_sides( PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: prepare_for_value_on_sides" ) ;
   PEL_CHECK_PRE( prepare_for_value_on_sides_PRE( fe ) ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: prepare_for_gradient_on_sides( PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: prepare_for_gradient_on_sides" ) ;
   PEL_CHECK_PRE( prepare_for_gradient_on_sides_PRE( fe ) ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK_PRE( prepare_for_value_on_bounds_PRE( fe ) ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: prepare_for_gradient_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: prepare_for_gradient_on_bounds" ) ;
   PEL_CHECK_PRE( prepare_for_gradient_on_bounds_PRE( fe ) ) ;
}

//-------------------------------------------------------------------------
void
FE_Parameter:: do_the_links( FE_SetOfParameters const* prms )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: do_the_links" ) ;
   PEL_CHECK_PRE( do_the_links_PRE( prms ) ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: cell_value( FE_TimeIterator const* t_it,
                           PDE_LocalFEcell const* fe,
                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: cell_value" ) ;
   PEL_CHECK_PRE( cell_value_PRE( t_it, fe, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "cell_value" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: cell_gradient( FE_TimeIterator const* t_it,
                              PDE_LocalFEcell const* fe,
                              size_t a,
                              size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: cell_gradient" ) ;
   PEL_CHECK_PRE( cell_gradient_PRE( t_it, fe, a, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "cell_gradient" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                 PDE_LocalFEcell const* fe,
                                 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "cell_value_at_IP" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                    PDE_LocalFEcell const* fe,
                                    size_t a,
                                    size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: cell_gradient_at_IP" ) ;
   PEL_CHECK_PRE( cell_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "cell_gradient_at_IP" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                 PDE_LocalFEcell const* fe,
                                 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "cell_value_at_pt" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: cell_gradient_at_pt( FE_TimeIterator const* t_it,
                                    PDE_LocalFEcell const* fe,
                                    size_t a,
                                    size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: cell_gradient_at_pt" ) ;
   PEL_CHECK_PRE( cell_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "cell_gradient_at_pp" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: bound_value( FE_TimeIterator const* t_it,
                            PDE_LocalFEbound const* fe,
                            size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: bound_value" ) ;
   PEL_CHECK_PRE( bound_value_PRE( t_it, fe, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "bound_value" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: bound_gradient( FE_TimeIterator const* t_it,
                               PDE_LocalFEbound const* fe,
                               size_t a,
                               size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: bound_gradient" ) ;
   PEL_CHECK_PRE( bound_gradient_PRE( t_it, fe, a, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "bound_gradient" ) ;
   return( PEL::bad_double() ) ;
}


//-------------------------------------------------------------------------
double
FE_Parameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                  PDE_LocalFEbound const* fe,
                                  size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "bound_value_at_IP" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: bound_gradient_at_IP( FE_TimeIterator const* t_it,
                                     PDE_LocalFEbound const* fe,
                                     size_t a,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: bound_gradient_at_IP" ) ;
   PEL_CHECK_PRE( bound_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "bound_gradient_at_IP" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                  PDE_LocalFEbound const* fe,
                                  size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "bound_value_at_pt" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: bound_gradient_at_pt( FE_TimeIterator const* t_it,
                                     PDE_LocalFEbound const* fe,
                                     size_t a,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: bound_gradient_at_pt" ) ;
   PEL_CHECK_PRE( bound_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "bound_gradient_at_pt" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: side_value( FE_TimeIterator const* t_it,
                           PDE_CursorFEside const* fe,
                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: side_value" ) ;
   PEL_CHECK_PRE( side_value_PRE( t_it, fe, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "side_value" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: side_gradient( FE_TimeIterator const* t_it,
                              PDE_CursorFEside const* fe,
                              size_t a,
                              size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: side_gradient" ) ;
   PEL_CHECK_PRE( side_gradient_PRE( t_it, fe, a, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "side_gradient" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: side_value_at_pt( FE_TimeIterator const* t_it,
                                 PDE_CursorFEside const* fe,
                                 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: side_value_at_pt" ) ;
   PEL_CHECK_PRE( side_value_at_pt_PRE( t_it, fe, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "side_value_at_pt" ) ;
   return( 0.0 ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: side_gradient_at_pt( FE_TimeIterator const* t_it,
                                    PDE_CursorFEside const* fe,
                                    size_t a,
                                    size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: side_gradient_at_pt" ) ;
   PEL_CHECK_PRE( side_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "side_gradient_at_pt" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: side_value_at_IP( FE_TimeIterator const* t_it,
                                 PDE_CursorFEside const* fe,
                                 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: side_value_at_IP" ) ;
   PEL_CHECK_PRE( side_value_at_IP_PRE( t_it, fe, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "side_value_at_IP" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
double
FE_Parameter:: side_gradient_at_IP( FE_TimeIterator const* t_it,
                                    PDE_CursorFEside const* fe,
                                    size_t a,
                                    size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_Parameter:: side_gradient_at_IP" ) ;
   PEL_CHECK_PRE( side_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "side_gradient_at_IP" ) ;
   return( PEL::bad_double() ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: is_a_prototype( void ) const
//-------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: do_the_links_PRE( FE_SetOfParameters const* prms ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( prms != 0 ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: reset_PRE( PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: prepare_for_value_on_cells_PRE(
                                          PDE_LocalFEcell const* fe ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( !fe->is_valid() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: prepare_for_gradient_on_cells_PRE(
                                          PDE_LocalFEcell const* fe ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( !fe->is_valid() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: prepare_for_value_on_sides_PRE(
                                        PDE_CursorFEside const * fe ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( !fe->is_valid() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: prepare_for_gradient_on_sides_PRE(
                                        PDE_CursorFEside const * fe ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( !fe->is_valid() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: prepare_for_value_on_bounds_PRE(
                                         PDE_LocalFEbound const* fe ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( !fe->is_valid() ) ;
   return( true ) ;
}

//-------------------------------------------------------------------------
bool
FE_Parameter:: prepare_for_gradient_on_bounds_PRE(
                                         PDE_LocalFEbound const* fe ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( !fe->is_valid() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: cell_value_PRE( FE_TimeIterator const* t_it,
                               PDE_LocalFEcell const* fe,
                               size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( ok_for_cell_calculations( fe, FE_Parameter::Val ) ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: cell_gradient_PRE( FE_TimeIterator const* t_it,
                                  PDE_LocalFEcell const* fe,
                                  size_t a,
                                  size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( ok_for_cell_calculations( fe, FE_Parameter::Grad ) ) ;
   PEL_ASSERT( a < fe->nb_space_dimensions() ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: cell_value_at_pt_PRE( FE_TimeIterator const* t_it,
                                     PDE_LocalFEcell const* fe,
                                     size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->calculation_point() != 0 ) ;
   PEL_ASSERT( ok_for_cell_calculations( fe, FE_Parameter::Val ) ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: cell_gradient_at_pt_PRE( FE_TimeIterator const* t_it,
                                        PDE_LocalFEcell const* fe,
                                        size_t a,
                                        size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->calculation_point() != 0 ) ;
   PEL_ASSERT( ok_for_cell_calculations( fe, FE_Parameter::Grad ) ) ;
   PEL_ASSERT( a < fe->nb_space_dimensions() ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: cell_value_at_IP_PRE( FE_TimeIterator const* t_it,
                                     PDE_LocalFEcell const* fe,
                                     size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->valid_IP() ) ;
   PEL_ASSERT( ok_for_cell_calculations( fe, FE_Parameter::Val ) ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: cell_gradient_at_IP_PRE( FE_TimeIterator const* t_it,
                                        PDE_LocalFEcell const* fe,
                                        size_t a,
                                        size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->valid_IP() ) ;
   PEL_ASSERT( ok_for_cell_calculations( fe, FE_Parameter::Grad ) ) ;
   PEL_ASSERT( a < fe->nb_space_dimensions() ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: bound_value_PRE( FE_TimeIterator const* t_it,
                                PDE_LocalFEbound const* fe,
                                size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( ok_for_bound_calculations( fe, FE_Parameter::Val ) ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: bound_gradient_PRE( FE_TimeIterator const* t_it,
                                   PDE_LocalFEbound const* fe,
                                   size_t a,
                                   size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( ok_for_bound_calculations( fe, FE_Parameter::Grad ) ) ;
   PEL_ASSERT( a < fe->nb_space_dimensions() ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: bound_value_at_pt_PRE( FE_TimeIterator const* t_it,
                                      PDE_LocalFEbound const* fe,
                                      size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->calculation_point() != 0 ) ;
   PEL_ASSERT( ok_for_bound_calculations( fe, FE_Parameter::Val ) ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: bound_gradient_at_pt_PRE( FE_TimeIterator const* t_it,
                                         PDE_LocalFEbound const* fe,
                                         size_t a,
                                         size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->calculation_point() != 0 ) ;
   PEL_ASSERT( ok_for_bound_calculations( fe, FE_Parameter::Grad ) ) ;
   PEL_ASSERT( a < fe->nb_space_dimensions() ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: bound_value_at_IP_PRE( FE_TimeIterator const* t_it,
                                      PDE_LocalFEbound const* fe,
                                      size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->valid_IP() ) ;
   PEL_ASSERT( ok_for_bound_calculations( fe, FE_Parameter::Val ) ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: bound_gradient_at_IP_PRE( FE_TimeIterator const* t_it,
                                         PDE_LocalFEbound const* fe,
                                         size_t a,
                                         size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->valid_IP() ) ;
   PEL_ASSERT( ok_for_bound_calculations( fe, FE_Parameter::Grad ) ) ;
   PEL_ASSERT( a < fe->nb_space_dimensions() ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: side_value_PRE( FE_TimeIterator const* t_it,
                               PDE_CursorFEside const* fe,
                               size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( ok_for_side_calculations( fe, FE_Parameter::Val ) ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: side_gradient_PRE( FE_TimeIterator const* t_it,
                                  PDE_CursorFEside const* fe,
                                  size_t a,
                                  size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( ok_for_side_calculations( fe, FE_Parameter::Grad ) ) ;
   PEL_ASSERT( a < fe->nb_space_dimensions() ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: side_value_at_pt_PRE( FE_TimeIterator const* t_it,
                                     PDE_CursorFEside const* fe,
                                     size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->calculation_point() != 0 ) ;
   PEL_ASSERT( ok_for_side_calculations( fe, FE_Parameter::Val ) ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: side_gradient_at_pt_PRE( FE_TimeIterator const* t_it,
                                        PDE_CursorFEside const* fe,
                                        size_t a,
                                        size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->calculation_point() != 0 ) ;
   PEL_ASSERT( ok_for_side_calculations( fe, FE_Parameter::Grad ) ) ;
   PEL_ASSERT( a < fe->nb_space_dimensions() ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: side_value_at_IP_PRE( FE_TimeIterator const* t_it,
                                     PDE_CursorFEside const* fe,
                                     size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->valid_IP() ) ;
   PEL_ASSERT( ok_for_side_calculations( fe, FE_Parameter::Val ) ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: side_gradient_at_IP_PRE( FE_TimeIterator const* t_it,
                                        PDE_CursorFEside const* fe,
                                        size_t a,
                                        size_t ic ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( fe->valid_IP() ) ;
   PEL_ASSERT( ok_for_side_calculations( fe, FE_Parameter::Grad ) ) ;
   PEL_ASSERT( a < fe->nb_space_dimensions() ) ;
   PEL_ASSERT( ic < nb_components() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: create_replica_PRE( PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   PEL_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( dom != 0 ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: create_replica_POST( FE_Parameter const* result,
                                    PEL_Object* a_owner,
                                    PDE_DomainAndFields const* dom,
                                    PEL_ModuleExplorer const* exp) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: create_replica_PRE( PEL_Object* a_owner,
                                   PDE_SetOfDomains const* sdoms,
                                   PEL_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( sdoms != 0 ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
FE_Parameter:: create_replica_POST( FE_Parameter const* result,
                                    PEL_Object* a_owner,
                                    PDE_SetOfDomains const* sdoms,
                                    PEL_ModuleExplorer const* exp) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
FE_Parameter:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
               PEL_ObjectRegister::create( PEL_Root::object(),
                                           "FE_Parameter descendant" ) ;
   return( result ) ;
}
