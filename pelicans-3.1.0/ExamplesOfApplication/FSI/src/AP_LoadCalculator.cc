#include <AP_LoadCalculator.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <LA_PelMatrix.hh>
#include <LA_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfDomains.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>

using std::cout ; using std::endl ;
using std::set ;
using std::string ; using std::ostringstream ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

struct AP_LoadCalculator_ERROR {
   static void n0( std::string const& a_name ) ;
   static void n1( std::string const& a_name ) ;
   static void n2( void ) ;
   static void n3( void ) ;
   static void n4( void ) ;
   static void n5( double xx1, double xx2, double dbl_eps, double dbl_min ) ;
} ;

std::map< std::string, AP_LoadCalculator* > AP_LoadCalculator::OBJS ;

AP_LoadCalculator const* 
AP_LoadCalculator::PROTOTYPE = new AP_LoadCalculator() ;

//---------------------------------------------------------------------------
AP_LoadCalculator:: AP_LoadCalculator( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "AP_LoadCalculator" )
{
}

//---------------------------------------------------------------------------
AP_LoadCalculator*
AP_LoadCalculator:: create_replica( PEL_Object* a_owner,
			            PDE_DomainAndFields const* dom,
			            FE_SetOfParameters const* prms,
			            PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   PEL_Error::object()->raise_plain( "only multidomains" ) ;

   AP_LoadCalculator* result = 0 ;
   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_LoadCalculator*
AP_LoadCalculator:: create_replica( PEL_Object* a_owner,
			         PDE_SetOfDomains const* sdoms,
			         FE_SetOfParameters const* prms,
			         PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, sdoms, prms, exp ) ) ;

   AP_LoadCalculator* result = 
                    new AP_LoadCalculator( a_owner, sdoms, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, sdoms, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
AP_LoadCalculator:: AP_LoadCalculator( 
                             PEL_Object* a_owner,
			     PDE_SetOfDomains const* sdoms,
			     FE_SetOfParameters const* prms,
			     PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, sdoms, exp )
   , NAME( exp->has_entry( "name") ? exp->string_data( "name" ) : "unknown" )
   , SD( 0 )
   , FV( 0 )
   , L_FV( exp->int_data( "level_of_fluid_velocity" ) )
   , FP( 0 )
   , L_FP( exp->int_data( "level_of_fluid_pressure" ) )
   , MU( prms->item( exp->string_data( "param_viscous" ) ) )
   , NB_COMPS( PEL::bad_index() )
   , NB_DIMS( PEL::bad_index() )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( GE_QRprovider::object( 
                         exp->string_data( "quadrature_rule_provider" ) ) )
   , bFE( 0 )
   , SURFCOL( GE_Color::object( 
        exp->string_data( "color_of_fluid_interface_with_structure" ) ) )
   , S2F( 0 )
   , MY_DBL_EPS( exp->has_entry( "dbl_epsilon" ) ?
                 exp->double_data( "dbl_epsilon" ) : 1.E-8 )
   , MY_DBL_MIN( exp->has_entry( "dbl_minimum" ) ?
                 exp->double_data( "dbl_minimum" ) : 1.E-8 )
   , F_LOAD( 0 )
   , S_LOAD( 0 )
{
   PEL_LABEL( "AP_LoadCalculator:: AP_LoadCalculator" ) ;

   if( OBJS.count( NAME ) != 0 ) AP_LoadCalculator_ERROR::n0( NAME ) ;
   OBJS[ NAME ] = this ;
   
   check_param_nb_components( MU, "param_viscous", 1 ) ;

   PDE_DomainAndFields const* s_dom = 
                  sdoms->domain( exp->string_data( "domain_of_structure" ) ) ;
   SD = s_dom->set_of_discrete_fields()->item( 
                  exp->string_data( "structure_displacement" ) ) ;
   NB_COMPS = SD->nb_components() ;

 
   PDE_DomainAndFields const* f_dom = 
                  sdoms->domain( exp->string_data( "domain_of_fluid" ) ) ;
   FV = f_dom->set_of_discrete_fields()->item( 
                  exp->string_data( "fluid_velocity" ) ) ;
   FP = f_dom->set_of_discrete_fields()->item( 
                  exp->string_data( "fluid_pressure" ) ) ;
   NB_DIMS = f_dom->nb_space_dimensions() ;

   check_field_storage_depth( FV, L_FV ) ;
   check_field_storage_depth( FP, L_FP ) ;
   if( FV->nb_components() != NB_COMPS ) AP_LoadCalculator_ERROR::n2() ;
   if( NB_COMPS != NB_DIMS ) AP_LoadCalculator_ERROR::n3() ;

   S2F = LA_PelMatrix::create( this, 
                               NB_COMPS * SD->nb_nodes(),
                               NB_COMPS * FV->nb_nodes() ) ;
   // not distributed
   F_LOAD = LA_SeqVector::create( this, NB_COMPS * FV->nb_nodes() ) ;
   S_LOAD = LA_SeqVector::create( this, NB_COMPS * SD->nb_nodes() ) ;

   bFE = f_dom->create_LocalFEbound( this ) ;
   bFE->require_field_calculation( FV, PDE_LocalFE::N ) ;
   bFE->require_field_calculation( FV, PDE_LocalFE::dN ) ;
   bFE->require_field_calculation( FP, PDE_LocalFE::N ) ;
   bFE->require_field_calculation( SD, PDE_LocalFE::N ) ;
   MU->transfer_bound_calculation_requirements( bFE, FE_Parameter::Val ) ;
}

//---------------------------------------------------------------------------
AP_LoadCalculator:: ~AP_LoadCalculator( void )
//---------------------------------------------------------------------------
{
}

//--------------------------------------------------------------------------
AP_LoadCalculator*
AP_LoadCalculator:: object( std::string const& a_name )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: object" ) ;
   std::map<std::string,AP_LoadCalculator*>::const_iterator it = 
                                                         OBJS.find( a_name ) ;
   if( it == OBJS.end() ) AP_LoadCalculator_ERROR::n1( a_name ) ;
   
   AP_LoadCalculator* result = (*it).second ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
std::string const&
AP_LoadCalculator:: name( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: name" ) ;

   return( NAME ) ;
}

//---------------------------------------------------------------------------
PDE_DiscreteField const* 
AP_LoadCalculator:: fluid_velocity( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: fluid_velocity" ) ;

   PDE_DiscreteField const* result = FV ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_DiscreteField const* 
AP_LoadCalculator:: structure_displacement( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: structure_displacement" ) ;

   PDE_DiscreteField const* result = SD ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
AP_LoadCalculator:: do_before_time_stepping( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;

   start_total_timer( "AP_LoadCalculator:: do_before_time_stepping" ) ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      if( SURFCOL->is_overlapping( bFE->color() ) )
      {
         if( bFE->nb_local_nodes( SD ) == 0 ) AP_LoadCalculator_ERROR::n4() ;
         bFE->set_row_and_col_fields( SD, SD ) ;

         ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(), 
                                 NB_COMPS,
                                 bFE->col_field_node_connectivity(),
                                 NB_COMPS ) ;

         for( size_t j=0 ; j<bFE->nb_local_nodes( FV ) ; ++j )
         {
            GE_Point const* f_pt =  bFE->local_node_location( FV, j ) ;
            if( bFE->polyhedron()->contains( f_pt ) )
            {
               size_t f_n = bFE->global_node( FV, j ) ;
               bFE->set_calculation_point( f_pt ) ;
               for( size_t i=0 ; i<bFE->nb_local_nodes( SD ) ; ++i )
               {
                  GE_Point const* s_pt = bFE->local_node_location( SD, i ) ;
                  if( bFE->polyhedron()->contains( s_pt ) )
                  {
                     size_t s_n = bFE->global_node( SD, i ) ;
                     double xx = bFE->N_at_pt( SD, i ) ;
                     for( size_t ic=0 ; ic<NB_COMPS ; ++ic )
                     {
                        size_t I = idx( SD, s_n, ic );
                        size_t J = idx( FV, f_n, ic );
                        if( S2F->item( I, J ) != 0.0 )
                        {
                           bool ok = PEL::double_equality( S2F->item( I, J ),
                                                           xx,
                                                           MY_DBL_EPS,
                                                           MY_DBL_MIN ) ;

                           if( !ok ) 
                              AP_LoadCalculator_ERROR::n5( S2F->item( I, J ),
                                                           xx, MY_DBL_EPS,
                                                           MY_DBL_MIN ) ;
                        }
                        S2F->set_item( I, J, xx ) ;
                     }
                  }
               }
            }
         }
      }
   }

   stop_total_timer() ;
}

//---------------------------------------------------------------------------
void 
AP_LoadCalculator:: do_one_inner_iteration( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "AP_LoadCalculator:: do_one_inner_iteration" ) ;

   doubleVector coef( NB_COMPS ) ;

   size_t nb_bds = 0 ;
   F_LOAD->nullify() ;
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      if( SURFCOL->is_overlapping( bFE->color() ) )
      {
         nb_bds++ ;
         bFE->set_row_and_col_fields( FV, FV ) ;
         ELEMENT_EQ->initialize( bFE->row_field_node_connectivity(),
                                 FV->nb_components(),
                                 bFE->col_field_node_connectivity(),
                                 FV->nb_components() ) ;

         bFE->start_IP_iterator( QRP ) ;
         for( ; bFE->valid_IP() ; bFE->go_next_IP() )
         { 
            GE_Vector const* nor = bFE->outward_normal() ;
            double mu = MU->bound_value_at_IP( t_it, bFE ) ;
            coef.set( 0.0 ) ;
            for( size_t d=0 ; d<NB_DIMS ; ++d )
            {
               coef( d ) = - bFE->value_at_IP( FP, L_FP ) * 
                             nor->component( d ) ;
               for( size_t i=0 ; i<NB_DIMS ; ++i )
               {
                  coef( d ) += ( bFE->gradient_at_IP( FV, L_FV, d, i ) 
 			       + bFE->gradient_at_IP( FV, L_FV, i, d ) )
                             *  mu * nor->component( i ) ;
               }
            }
            FE::add_row( ELEMENT_EQ, bFE, coef ) ;
         }

         for( size_t il=0 ; il<ELEMENT_EQ->nb_rows() ; ++il )
         {
            for( size_t ic=0 ; ic<ELEMENT_EQ->nb_row_sub_indices() ; ++ic )
            {
               size_t n = ELEMENT_EQ->row_node( il ) ;
               size_t J = idx( FV, n, ic ) ;
               F_LOAD->add_to_item( J, - ELEMENT_EQ->vector_item( il, ic ) ) ;
            }     
         }
      }
   }
   
   F_LOAD->synchronize() ;
   S2F->synchronize() ;
   S2F->multiply_vec_then_add( F_LOAD, S_LOAD, 1.0 , 0.0 ) ;

   if( verbose_level() >= 2 )
   {
      PEL::out() << indent() << "   found " << nb_bds << " bounds" << endl ;
   }
   stop_total_timer() ;
}

//---------------------------------------------------------------------------
double
AP_LoadCalculator:: fluid_load( size_t n, size_t ic ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: fluid_load" ) ;
   PEL_CHECK_PRE( n < fluid_velocity()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < fluid_velocity()->nb_components() ) ;

   double result = F_LOAD->item( idx( FV, n, ic ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
double
AP_LoadCalculator:: solid_load( size_t n, size_t ic ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: solid_load" ) ;
   PEL_CHECK_PRE( n < structure_displacement()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < structure_displacement()->nb_components() ) ;

   double result = S_LOAD->item( idx( SD, n, ic ) ) ;
   return( result ) ;
}
//---------------------------------------------------------------------------
void
AP_LoadCalculator:: update_load( PDE_LinkDOF2Unknown const* link,
                                 LA_SeqVector* load ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: update_load" ) ;
   PEL_CHECK_PRE( link != 0 ) ;
   PEL_CHECK_PRE( link->field() == structure_displacement() ) ;
   PEL_CHECK_PRE( load != 0 ) ;
   PEL_CHECK_PRE( load->nb_rows() == link->unknown_vector_size() ) ;

   size_t ifin = 0 ;
   for( size_t n=0 ; n<SD->nb_nodes() ; ++n )
   {
      for ( size_t ic=0 ; ic<NB_COMPS ; ++ic )
      {
         size_t ii = idx( SD, n, ic ) ;
         if( link->DOF_is_unknown( n, ic ) )
         {
            ++ifin ;
            size_t i_unk = link->unknown_linked_to_DOF( n, ic ) ;
            load->set_item( i_unk, S_LOAD->item( ii ) ) ;
         }
      }
   }
   PEL_ASSERT( ifin == link->unknown_vector_size() ) ;
}

//---------------------------------------------------------------------------
size_t
AP_LoadCalculator:: idx( PDE_DiscreteField const* ff,
                         size_t n, size_t ic ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LoadCalculator:: idx" ) ;
   PEL_CHECK( n < ff->nb_nodes() ) ;
   PEL_CHECK( ic < ff->nb_components() ) ;

   return( n*NB_COMPS + ic ) ;
}

//internal-------------------------------------------------------------------
void
AP_LoadCalculator_ERROR:: n0( std::string const& a_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "AP_LoadCalculator :" << endl ;
   mesg << "   attempt to create two instances with the" << endl ;
   mesg << "   same name :\"" << a_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
AP_LoadCalculator_ERROR:: n1( std::string const& a_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "AP_LoadCalculator :" << endl ;
   mesg << "   there is no recorded instance of name :\"" << a_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
AP_LoadCalculator_ERROR:: n2( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "AP_LoadCalculator :" << endl ;
   mesg << "   the fields \"structure_displacement\" and \"fluid_velocity\" " 
        << endl ;
   mesg << "   should have the same number of components" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
AP_LoadCalculator_ERROR:: n3( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "AP_LoadCalculator :" << endl ;
   mesg << "   the number of components of the field \"fluid_velocity\"" 
        << endl ;
   mesg << "   should be equal to the number of space dimensions of the"
        << endl ;
   mesg << "   fluid domain" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
AP_LoadCalculator_ERROR:: n4( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "AP_LoadCalculator :" << endl ;
   mesg << "   there are bounds of the declared fluid-structure interface with"
        << endl
        << "   no local discretization of the structure displacement" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
AP_LoadCalculator_ERROR:: n5( double xx1, double xx2, 
                              double dbl_eps, double dbl_min )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << std::setprecision( 10 )
        << std::setiosflags( std::ios::scientific ) ;
   mesg << "AP_LoadCalculator :" << endl ;
   mesg << "   the two values : " << xx1 << endl ;
   mesg << "              and : " << xx2 << endl ;
   mesg << "   should be equal." << endl ;
   mesg << "If it is not an error, do adjust the values of the keywords " 
        << endl ;
   mesg << "\"dbl_epsilon\" (currently : " << dbl_eps << ") and" << endl ;
   mesg << "\"dbl_minimum\" (currently : " << dbl_min << ")" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

