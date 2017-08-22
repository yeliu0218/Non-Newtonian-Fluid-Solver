#include <MY_MUSCL_Scheme.hh>

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

#include <MY_MUSCL_DataStructure.hh>

// C++ tools :

#include <iostream>

MY_MUSCL_Scheme const* 
MY_MUSCL_Scheme::PROTOTYPE = new MY_MUSCL_Scheme( "MY_MUSCL_Scheme" ) ;

//---------------------------------------------------------------------------
MY_MUSCL_Scheme:: MY_MUSCL_Scheme( std::string const& a_name )
//---------------------------------------------------------------------------
   : MY_AdvectiveScheme( a_name )
   , MDS( 0 )
{
}

//---------------------------------------------------------------------------
MY_AdvectiveScheme*
MY_MUSCL_Scheme:: create_replica( PEL_Object* a_owner,
				  PDE_DomainAndFields const* dom,
				  FE_SetOfParameters const* prms,
				  PDE_DiscreteField const* field,
				  PEL_ModuleExplorer const* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_Scheme:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, field, exp ) ) ;
  
   MY_MUSCL_Scheme* result 
                     = new MY_MUSCL_Scheme( a_owner, dom, prms, field, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, field, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
MY_MUSCL_Scheme:: MY_MUSCL_Scheme( PEL_Object* a_owner,
				   PDE_DomainAndFields const* dom,
				   FE_SetOfParameters const* prms,
				   PDE_DiscreteField const* field,
                   PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : MY_AdvectiveScheme( a_owner, dom, prms, field, exp )
   , MDS( MY_MUSCL_DataStructure::object( dom ) )
{
   PEL_LABEL( "MY_MUSCL_Scheme:: MY_MUSCL_Scheme" ) ;
}

//---------------------------------------------------------------------------
MY_MUSCL_Scheme:: ~MY_MUSCL_Scheme( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_Scheme:: ~MY_MUSCL_Scheme" ) ;
   
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
MY_MUSCL_Scheme:: assemble_convection_on_side( FE_TimeIterator const* t_it,
					       PDE_CursorFEside* fe,
					       PDE_LocalEquation* leq )
//----------------------------------------------------------------------
{
    PEL_LABEL( "MY_MUSCL_Scheme:: assemble_convection_on_side" ) ;
    PEL_CHECK_PRE(assemble_convection_on_side_PRE(t_it,fe,leq)) ;

    PDE_DiscreteField const* f = discrete_field() ;

    size_t un = MDS->upstream_node( fe, f ) ;
    size_t dn = MDS->downstream_node( fe, f ) ;
    double ud = MDS->upstream_node_dist( fe, f ) ;
    double dd = MDS->downstream_node_dist( fe, f ) ;
    double crg = (f->DOF_value(0,dn)-f->DOF_value(0,un))/( ud + dd ) ; 
    double const mass_flux = pour_rate( t_it, fe ) ;
    if( mass_flux>0. )
    {
       // Implicit value, upwind first order flux
       leq->add_to_matrix(  mass_flux, 0, 0 ) ;
       leq->add_to_matrix( -mass_flux, 1, 0 ) ;

       // Explicit value, upwind second order correction
       if( MDS->has_up_upstream_node( fe, f ) )
       {
          size_t uun = MDS->up_upstream_node( fe, f ) ;
 	  double lhg = ( f->DOF_value( 0, un ) - f->DOF_value( 0, uun ) )/
             MDS->up_upstream_node_dist( fe, f )  ;
	  double limit = limit_slope_ratio( crg, lhg ) ;
	  double contrib = -ud*mass_flux*limit*lhg ;
	  leq->add_to_vector(  contrib, 0 ) ;
	  leq->add_to_vector( -contrib, 1 ) ;
       }
    }
    else if( mass_flux<0. )
    {
       // Implicit value, upwind first order flux
       leq->add_to_matrix(  mass_flux, 0, 1 ) ;
       leq->add_to_matrix( -mass_flux, 1, 1 ) ;

       // Explicit value, upwind second order correction
       if( MDS->has_down_downstream_node( fe, f ) )
       {
          size_t ddn = MDS->down_downstream_node( fe, f ) ;
	  double rhg = ( f->DOF_value( 0, ddn ) - f->DOF_value( 0, dn ) )/
               MDS->down_downstream_node_dist( fe, f )   ;
	  double limit = limit_slope_ratio( crg, rhg ) ;
	  double contrib = dd*mass_flux*limit*rhg ;
	  leq->add_to_vector(  contrib, 0 ) ;
	  leq->add_to_vector( -contrib, 1 ) ;
       }
    }
}

//----------------------------------------------------------------------
void
MY_MUSCL_Scheme:: assemble_convection_on_bound( FE_TimeIterator const* t_it,
						PDE_LocalFEbound* fe,
						double const bound_value,
						PDE_LocalEquation* leq )
//----------------------------------------------------------------------
{
    PEL_LABEL( "MY_MUSCL_Scheme:: assemble_convection_on_bound" ) ;
    PEL_CHECK_PRE(assemble_convection_on_bound_PRE(t_it,fe,bound_value,leq)) ;
    
    double const mass_flux = pour_rate( t_it, fe ) ;
    if( mass_flux>0. )
    {
       leq->add_to_matrix( mass_flux, 0, 0 ) ;

       // Explicit value, upwind second order correction
       PDE_DiscreteField const* f = discrete_field() ;
       if( MDS->has_up_upstream_node( fe, f ) )
       {
       	  size_t un = MDS->upstream_node( fe, f ) ;
	  double ud = MDS->upstream_node_dist( fe, f ) ;
	 double crg = (bound_value-f->DOF_value(0,un))/ud ; 
         size_t uun = MDS->up_upstream_node( fe, f ) ;
 	 double lhg = ( f->DOF_value( 0, un ) - f->DOF_value( 0, uun ) ) / MDS->up_upstream_node_dist( fe, f )  ;
	  double limit = limit_slope_ratio( crg, lhg ) ;
	  leq->add_to_vector( -ud*mass_flux*limit*lhg, 0 ) ;

          // if Neumann homogeneous bound_value=ud and the MUSCL correction
          // vanishes -> upwind 1er ordre 
      }
    }
    else if( mass_flux<0. )
    {
       leq->add_to_vector( -mass_flux*bound_value, 0 ) ;
    }
}

//----------------------------------------------------------------------
double
MY_MUSCL_Scheme:: limit_slope_ratio( double numerator_slope, 
				     double denominator_slope,
                                     double smallest_nonzero_den ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_Scheme:: limit_slope_ratio" ) ;

   double result = 0. ;
   if( PEL::abs( denominator_slope )>smallest_nonzero_den )
   {
      double const ratio = numerator_slope/denominator_slope ;
      if( ratio>0. )
      {
         result = 2.*ratio/(1.+ratio) ;
      }
      //result = 1.5*(ratio)*(ratio+1.0)/(ratio*ratio+ratio+1.0);
   }
   return( result ) ;
}


//----------------------------------------------------------------------
void
MY_MUSCL_Scheme:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_Scheme:: print" ) ;
   
   std::string const s( indent_width, ' ' ) ;
   os << s << "convective scheme: \"MY_MUSCL_Scheme\"" << std::endl ;
}

//----------------------------------------------------------------------
bool
MY_MUSCL_Scheme:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( MY_AdvectiveScheme::invariant() ) ;
   return( true ) ;
}
