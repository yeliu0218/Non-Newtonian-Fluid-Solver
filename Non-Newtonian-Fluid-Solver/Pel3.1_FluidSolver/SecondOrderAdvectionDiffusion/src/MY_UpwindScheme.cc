#include <MY_UpwindScheme.hh>

// PELICANS tools :

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalEquation.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <iostream>

MY_UpwindScheme const* 
MY_UpwindScheme::PROTOTYPE = new MY_UpwindScheme( "MY_UpwindScheme" ) ;

//---------------------------------------------------------------------------
MY_UpwindScheme:: MY_UpwindScheme( std::string const& a_name )
//---------------------------------------------------------------------------
      : MY_AdvectiveScheme( a_name )
{
}

//---------------------------------------------------------------------------
MY_AdvectiveScheme*
MY_UpwindScheme:: create_replica( 
                           PEL_Object* a_owner,
                           PDE_DomainAndFields const* dom,
                           FE_SetOfParameters const* prms,
                           PDE_DiscreteField const* field,
                           PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_UpwindScheme:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, field, exp ) ) ;
  
   MY_UpwindScheme* result = new MY_UpwindScheme( a_owner, dom, prms, field, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, field, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
MY_UpwindScheme:: MY_UpwindScheme( 
                           PEL_Object* a_owner,
                           PDE_DomainAndFields const* dom,
                           FE_SetOfParameters const* prms,
                           PDE_DiscreteField const* field,
                           PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
  : MY_AdvectiveScheme( a_owner, dom, prms, field, exp )
{
}

//---------------------------------------------------------------------------
MY_UpwindScheme:: ~MY_UpwindScheme( void )
//---------------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
MY_UpwindScheme:: assemble_convection_on_side( FE_TimeIterator const* t_it,
					       PDE_CursorFEside* fe,
					       PDE_LocalEquation* leq )
//----------------------------------------------------------------------
{
    PEL_LABEL( "MY_UpwindScheme:: assemble_convection_on_side" ) ;
    PEL_CHECK_PRE( assemble_convection_on_side_PRE( t_it,  fe, leq ) ) ;

    double const mass_flux = pour_rate( t_it, fe ) ;
    if( mass_flux>0. )
    {
       leq->add_to_matrix(  mass_flux, 0, 0 ) ;
       leq->add_to_matrix( -mass_flux, 1, 0 ) ;
    }
    else if( mass_flux<0. )
    {
       leq->add_to_matrix( -mass_flux, 1, 1 ) ;
       leq->add_to_matrix(  mass_flux, 0, 1 ) ;
    }
}

//----------------------------------------------------------------------
void
MY_UpwindScheme:: assemble_convection_on_bound( FE_TimeIterator const* t_it,
						PDE_LocalFEbound* fe,
						double const bound_value,
						PDE_LocalEquation* leq )
//----------------------------------------------------------------------
{
    PEL_LABEL( "MY_UpwindScheme:: assemble_convection_on_bound" ) ;
    PEL_CHECK_PRE( assemble_convection_on_bound_PRE( t_it, fe, bound_value, leq ) ) ;
    
    double const mass_flux = pour_rate( t_it, fe ) ;
    if( mass_flux>0. )
    {
       leq->add_to_matrix( mass_flux, 0, 0 ) ;
    }
    else if( mass_flux<0. )
    {
       leq->add_to_vector( -mass_flux*bound_value, 0 ) ;
    }
}

//----------------------------------------------------------------------
void
MY_UpwindScheme:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "MY_UpwindScheme:: print" ) ;
   
   std::string const s( indent_width, ' ' ) ;
   os << s << "convective scheme : \"MY_UpwindScheme\"" << std::endl ;
}

//----------------------------------------------------------------------
bool
MY_UpwindScheme:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( MY_AdvectiveScheme::invariant() ) ;
   return( true ) ;
}
