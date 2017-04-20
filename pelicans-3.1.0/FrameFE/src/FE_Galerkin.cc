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

#include <FE_Galerkin.hh>

#include <GE_QRprovider.hh>

#include <PDE_DiscreteField.hh>
#include <PEL_Error.hh>
#include <PDE_LocalFEcell.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

#include <FE_TimeIterator.hh>

//----------------------------------------------------------------------
FE_Galerkin*
FE_Galerkin:: make( PEL_Object* a_owner,
		    PDE_DomainAndFields const* dom,
                    FE_SetOfParameters const* prms,
                    PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_Galerkin:: make" ) ;
   PEL_CHECK_PRE( dom != 0 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   std::string name = exp->string_data( "concrete_name" ) ;
   FE_Galerkin const* proto =
      static_cast<FE_Galerkin const*>( galerkin_map()->item( name ) ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;
      
   FE_Galerkin* result = proto->create_replica( a_owner, dom, prms, exp ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_Galerkin:: FE_Galerkin( PEL_Object* a_owner,
                           PDE_DomainAndFields const* dom,
                           PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , CV_FIELDS( PEL_Vector::create( this, 0 ) )
   , L_MASKED( PEL::bad_index() )
{
   PEL_LABEL( "FE_Galerkin:: FE_Galerkin" ) ;
}

//----------------------------------------------------------------------
FE_Galerkin:: FE_Galerkin( std::string const& name )
//----------------------------------------------------------------------
   : FE_OneStepIteration( name )
   , CV_FIELDS( 0 )
   , L_MASKED( 0 )
{
   PEL_LABEL( "FE_Galerkin:: FE_Galerkin" ) ;

   galerkin_map()->register_item( name, this ) ;
   
   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
FE_Galerkin:: ~FE_Galerkin( void ) 
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_Galerkin:: add_one_convected_field( PDE_DiscreteField* ff,
                                        size_t a_masked_level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_Galerkin:: add_one_convected_field" ) ;
   PEL_CHECK_PRE( !is_a_prototype() ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( masked_level()==PEL::bad_index() ||
                  masked_level()==a_masked_level ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( L_MASKED != PEL::bad_index() )
   {
      PEL_ASSERT( L_MASKED == a_masked_level ) ;
   }
   else
   {
      L_MASKED = a_masked_level ;
   }

   CV_FIELDS->extend( ff ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( masked_level()==a_masked_level ) ;
   PEL_CHECK_POST( convected_fields()->has( ff ) ) ;
}

//----------------------------------------------------------------------
PEL_Sequence const*
FE_Galerkin:: convected_fields( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_Galerkin:: convected_fields" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( CV_FIELDS ) ;
}

//----------------------------------------------------------------------
size_t
FE_Galerkin:: masked_level( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_Galerkin:: masked_level" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( L_MASKED ) ;
}

//----------------------------------------------------------------------
bool
FE_Galerkin:: QRprovider_for_material_derivative_POST( 
                                         GE_QRprovider const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->is_under_ownership_of( PEL_Root::object() ) );
   return( true ) ;
}

//----------------------------------------------------------------------
bool
FE_Galerkin:: transfer_calculation_requirements_for_material_derivative_PRE( 
                                         PDE_LocalFEcell const* fe ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( !fe->is_valid() );
   return( true ) ;
}

//-----------------------------------------------------------------------
bool
FE_Galerkin:: reset_discrete_problem_PRE( FE_TimeIterator const* t_it ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->is_finished() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------
bool
FE_Galerkin:: build_cell_contribution_to_material_derivative_PRE(
                                                FE_TimeIterator const* t_it,
                                                PDE_LocalFEcell* fe ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->is_finished() ) ;
   PEL_ASSERT( fe != 0 ) ;
   PEL_ASSERT( fe->is_valid() ) ;
   PEL_ASSERT( !fe->valid_IP() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------
bool
FE_Galerkin:: terminate_discrete_problem_PRE( 
                                         FE_TimeIterator const* t_it ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->is_finished() ) ;
   return( true ) ;
}

//-----------------------------------------------------------------------
bool
FE_Galerkin:: invariant( void ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( FE_OneStepIteration::invariant() ) ;

   PEL_ASSERT( is_a_prototype() ||
               ( masked_level() != PEL::bad_index()
                 || convected_fields()->index_limit() == 0 ) ) ;

   PEL_ASSERT(  is_a_prototype() ||
               (masked_level() == PEL::bad_index()
               || convected_fields()->index_limit() > 0 ) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
FE_Galerkin:: galerkin_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "FE_Galerkin descendant" ) ;
   return( result ) ;
}
